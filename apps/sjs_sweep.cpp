// apps/sjs_sweep.cpp
//
// Parameter sweep harness (SIGMOD-style experiments).
//
// This app repeatedly calls the single-run protocol (same as sjs_run),
// but over a grid/list of parameters, and writes:
//   (1) sweep_raw.csv     : one row per run (repeat)
//   (2) sweep_summary.csv : aggregated stats per (dataset params, method, variant, t)
//
// Main input is a JSON sweep file (no external JSON dependency; we parse a
// small JSON subset). You can keep your sweep files under config/sweeps/*.json.
//
// Usage:
//   ./sjs_sweep --config=config/sweeps/alpha.json
//
// Minimal JSON schema (all fields optional; unspecified values fall back to
// Config defaults):
// {
//   "base": {
//     "dataset": {
//       "source": "synthetic",          // synthetic|binary|csv
//       "name": "alpha_sweep",
//       "dim": 2,
//       "path_r": "", "path_s": "",
//       "synthetic": {
//         "generator": "stripe",        // stripe|uniform|clustered|hetero_sizes
//         "n_r": 200000,
//         "n_s": 200000,
//         "alpha": 1e-6,
//         "seed": 1,                    // gen_seed
//         "params": { "gap_factor": 0.1 }
//       }
//     },
//     "run": {
//       "method": "ours",
//       "variant": "sampling",
//       "t": 10000,
//       "seed": 1,
//       "repeats": 3,
//       "enum_cap": 0,
//       "j_star": 0,
//       "write_samples": false,
//       "extra": { "csv_sep": "tab" }
//     },
//     "output": { "out_dir": "out/alpha_sweep" },
//     "logging": { "level": "info" },
//     "sys": { "threads": 1 }
//   },
//
//   "sweep": {
//     "alpha":   [1e-8, 3e-8, 1e-7, 3e-7, 1e-6],
//     "t":       [10000, 100000],
//     "method":  ["ours", "aabb"],
//     "variant": ["sampling", "enum_sampling", "adaptive"],
//     "seed":    [1,2,3]              // optional; overrides repeats if provided
//   },
//
//   "files": {
//     "raw": "sweep_raw.csv",
//     "summary": "sweep_summary.csv"
//   }
// }
//
// Notes:
//   - For synthetic datasets, dataset is re-generated for each
//     (n_r,n_s,alpha,gen_seed) combination.
//   - For binary/csv datasets, dataset is loaded once; only method/variant/t/seed vary.

#include "baselines/baseline_factory_2d.h"

#include "sjs/baselines/runners/adaptive_runner.h"
#include "sjs/baselines/runners/enum_sampling_runner.h"
#include "sjs/baselines/runners/sampling_runner.h"

#include "sjs/core/config.h"
#include "sjs/core/logging.h"
#include "sjs/core/stats.h"
#include "sjs/core/timer.h"
#include "sjs/core/types.h"

#include "sjs/data/synthetic/generator.h"
#include "sjs/io/binary_io.h"
#include "sjs/io/csv_io.h"
#include "sjs/io/dataset.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace sjs {
namespace apps {

namespace {

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool IsHelpRequested(const sjs::ArgMap& args) {
  return args.Has("help") || args.Has("h") || args.Has("-h") || args.Has("--help");
}

inline bool EnsureDir(const fs::path& p, std::string* err) {
  try {
    if (p.empty()) return true;
    fs::create_directories(p);
    return true;
  } catch (const std::exception& e) {
    SetErr(err, std::string("create_directories failed: ") + e.what());
    return false;
  }
}

inline std::string SanitizeFilename(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    const bool ok = (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_' || c == '-' || c == '.';
    out.push_back(ok ? c : '_');
  }
  if (out.empty()) out = "x";
  return out;
}

inline std::string FormatSci(double x, int prec = 3) {
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(prec) << x;
  std::string s = oss.str();
  // Make it filename-friendly.
  for (char& c : s) {
    if (c == '+') c = 'p';
    if (c == '-') c = 'm';
    if (c == '.') c = '_';
  }
  return s;
}

// --------------------------
// Dataset load/generate (same logic as sjs_run)
// --------------------------

template <int Dim>
bool LoadOrGenerateDataset(const sjs::Config& cfg,
                           sjs::Dataset<Dim, sjs::Scalar>* out,
                           sjs::synthetic::Report* gen_report,
                           std::string* err) {
  if (!out) {
    SetErr(err, "LoadOrGenerateDataset: out is null");
    return false;
  }

  if (cfg.dataset.dim != Dim) {
    SetErr(err, "Config dim mismatch: cfg.dataset.dim != compiled Dim");
    return false;
  }

  if (cfg.dataset.source == sjs::DataSource::Synthetic) {
    sjs::synthetic::DatasetSpec spec;
    spec.name = cfg.dataset.name;
    spec.n_r = cfg.dataset.synthetic.n_r;
    spec.n_s = cfg.dataset.synthetic.n_s;
    spec.alpha = cfg.dataset.synthetic.alpha;
    spec.seed = cfg.dataset.synthetic.seed;
    spec.params = cfg.dataset.synthetic.extra;

    sjs::synthetic::Report rep;
    std::string local_err;
    if (!sjs::synthetic::GenerateDataset<Dim, sjs::Scalar>(
            cfg.dataset.synthetic.generator, spec, out, &rep, &local_err)) {
      SetErr(err, "Synthetic generation failed: " + local_err);
      return false;
    }
    if (gen_report) *gen_report = rep;
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::Binary) {
    if (cfg.dataset.path_r.empty() || cfg.dataset.path_s.empty()) {
      SetErr(err, "Binary dataset requires --path_r and --path_s");
      return false;
    }
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;
    sjs::binary::RelationFileInfo infoR, infoS;
    sjs::binary::BinaryReadOptions opt;
    opt.generate_ids_if_missing = true;
    opt.drop_empty = false;
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, &infoR, opt, &local_err)) {
      SetErr(err, "Failed reading binary R: " + local_err);
      return false;
    }
    if (!sjs::binary::ReadRelationBinary<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, &infoS, opt, &local_err)) {
      SetErr(err, "Failed reading binary S: " + local_err);
      return false;
    }
    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "binary" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "Binary dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  if (cfg.dataset.source == sjs::DataSource::CSV) {
    if (cfg.dataset.path_r.empty() || cfg.dataset.path_s.empty()) {
      SetErr(err, "CSV dataset requires --path_r and --path_s");
      return false;
    }
    sjs::Relation<Dim, sjs::Scalar> R, S;
    std::string local_err;

    // Simple reader uses sep=',' by default. You can change to TSV by passing --csv_sep=tab.
    char sep = ',';
    if (auto v = cfg.run.extra.find("csv_sep"); v != cfg.run.extra.end()) {
      if (!v->second.empty()) {
        if (v->second == "tab" || v->second == "\\t") sep = '\t';
        else sep = v->second[0];
      }
    }

    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_r, &R, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV R: " + local_err);
      return false;
    }
    if (!sjs::csv::ReadBoxesSimple<Dim, sjs::Scalar>(cfg.dataset.path_s, &S, sep, /*has_header=*/true, &local_err)) {
      SetErr(err, "Failed reading CSV S: " + local_err);
      return false;
    }

    sjs::Dataset<Dim, sjs::Scalar> ds;
    ds.name = cfg.dataset.name.empty() ? "csv" : cfg.dataset.name;
    ds.half_open = true;
    ds.R = std::move(R);
    ds.S = std::move(S);

    std::string v;
    if (!ds.Validate(/*require_proper=*/true, &v)) {
      SetErr(err, "CSV dataset validation failed: " + v);
      return false;
    }
    *out = std::move(ds);
    return true;
  }

  SetErr(err, "Unsupported dataset_source");
  return false;
}

// --------------------------
// CSV headers
// --------------------------

inline std::vector<std::string> RawHeader() {
  return {
      "dataset",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "rep",
      "seed",
      "ok",
      "wall_ms",
      "count_value",
      "count_exact",
      "used_enumeration",
      "enum_truncated",
      "enum_cap",
      "adaptive_branch",
      "adaptive_pilot_pairs",
      "phases_json",
      "note",
      "config_json",
      "gen_report_json",
  };
}

inline std::vector<std::string> SummaryHeader() {
  return {
      "dataset",
      "generator",
      "alpha",
      "n_r",
      "n_s",
      "method",
      "variant",
      "t",
      "repeats",
      "ok_rate",
      "wall_mean_ms",
      "wall_stdev_ms",
      "wall_median_ms",
      "wall_p95_ms",
      "count_mean",
      "count_stdev",
      "exact_frac",
      "note",
  };
}

// --------------------------
// Minimal JSON parser (subset)
// --------------------------

namespace json {

enum class Type : u8 { Null, Bool, Number, String, Array, Object };

struct Value {
  Type type{Type::Null};
  bool b{false};
  double num{0.0};
  std::string str;
  std::vector<Value> arr;
  std::unordered_map<std::string, Value> obj;

  bool IsNull() const { return type == Type::Null; }
  bool IsBool() const { return type == Type::Bool; }
  bool IsNumber() const { return type == Type::Number; }
  bool IsString() const { return type == Type::String; }
  bool IsArray() const { return type == Type::Array; }
  bool IsObject() const { return type == Type::Object; }
};

class Parser {
 public:
  explicit Parser(std::string_view s) : s_(s) {}

  bool Parse(Value* out, std::string* err) {
    if (!out) {
      SetErr(err, "json::Parser: out is null");
      return false;
    }
    SkipWs();
    if (!ParseValue(out, err)) return false;
    SkipWs();
    if (i_ != s_.size()) {
      SetErr(err, "json parse: trailing characters");
      return false;
    }
    return true;
  }

 private:
  void SkipWs() {
    while (i_ < s_.size()) {
      char c = s_[i_];
      if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
        ++i_;
      } else {
        break;
      }
    }
  }

  bool Match(std::string_view kw) {
    if (s_.substr(i_, kw.size()) == kw) {
      i_ += kw.size();
      return true;
    }
    return false;
  }

  bool ParseValue(Value* out, std::string* err) {
    SkipWs();
    if (i_ >= s_.size()) {
      SetErr(err, "json parse: unexpected end");
      return false;
    }

    const char c = s_[i_];
    if (c == '{') return ParseObject(out, err);
    if (c == '[') return ParseArray(out, err);
    if (c == '"') return ParseString(out, err);
    if (c == 't') {
      if (!Match("true")) { SetErr(err, "json parse: expected true"); return false; }
      out->type = Type::Bool;
      out->b = true;
      return true;
    }
    if (c == 'f') {
      if (!Match("false")) { SetErr(err, "json parse: expected false"); return false; }
      out->type = Type::Bool;
      out->b = false;
      return true;
    }
    if (c == 'n') {
      if (!Match("null")) { SetErr(err, "json parse: expected null"); return false; }
      out->type = Type::Null;
      return true;
    }
    // number
    if (c == '-' || (c >= '0' && c <= '9')) {
      return ParseNumber(out, err);
    }

    SetErr(err, std::string("json parse: unexpected char '") + c + "'");
    return false;
  }

  bool ParseString(Value* out, std::string* err) {
    if (s_[i_] != '"') {
      SetErr(err, "json parse: expected string");
      return false;
    }
    ++i_;
    std::string res;
    while (i_ < s_.size()) {
      char c = s_[i_++];
      if (c == '"') {
        out->type = Type::String;
        out->str = std::move(res);
        return true;
      }
      if (c == '\\') {
        if (i_ >= s_.size()) { SetErr(err, "json parse: bad escape"); return false; }
        char e = s_[i_++];
        switch (e) {
          case '"': res.push_back('"'); break;
          case '\\': res.push_back('\\'); break;
          case '/': res.push_back('/'); break;
          case 'b': res.push_back('\b'); break;
          case 'f': res.push_back('\f'); break;
          case 'n': res.push_back('\n'); break;
          case 'r': res.push_back('\r'); break;
          case 't': res.push_back('\t'); break;
          case 'u': {
            // Minimal handling: skip 4 hex digits and emit '?' (good enough for configs).
            if (i_ + 4 > s_.size()) { SetErr(err, "json parse: bad \\u escape"); return false; }
            i_ += 4;
            res.push_back('?');
            break;
          }
          default:
            SetErr(err, "json parse: unsupported escape");
            return false;
        }
      } else {
        res.push_back(c);
      }
    }
    SetErr(err, "json parse: unterminated string");
    return false;
  }

  bool ParseNumber(Value* out, std::string* err) {
    const usize start = i_;
    if (s_[i_] == '-') ++i_;
    while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    if (i_ < s_.size() && s_[i_] == '.') {
      ++i_;
      while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    }
    if (i_ < s_.size() && (s_[i_] == 'e' || s_[i_] == 'E')) {
      ++i_;
      if (i_ < s_.size() && (s_[i_] == '+' || s_[i_] == '-')) ++i_;
      while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
    }

    const std::string token(s_.substr(start, i_ - start));
    try {
      out->type = Type::Number;
      out->num = std::stod(token);
      return true;
    } catch (...) {
      SetErr(err, "json parse: bad number");
      return false;
    }
  }

  bool ParseArray(Value* out, std::string* err) {
    if (s_[i_] != '[') { SetErr(err, "json parse: expected ["); return false; }
    ++i_;
    SkipWs();
    out->type = Type::Array;
    out->arr.clear();

    if (i_ < s_.size() && s_[i_] == ']') {
      ++i_;
      return true;
    }

    while (true) {
      Value v;
      if (!ParseValue(&v, err)) return false;
      out->arr.push_back(std::move(v));
      SkipWs();
      if (i_ >= s_.size()) { SetErr(err, "json parse: unterminated array"); return false; }
      if (s_[i_] == ',') {
        ++i_;
        SkipWs();
        continue;
      }
      if (s_[i_] == ']') {
        ++i_;
        return true;
      }
      SetErr(err, "json parse: expected , or ] in array");
      return false;
    }
  }

  bool ParseObject(Value* out, std::string* err) {
    if (s_[i_] != '{') { SetErr(err, "json parse: expected {"); return false; }
    ++i_;
    SkipWs();
    out->type = Type::Object;
    out->obj.clear();

    if (i_ < s_.size() && s_[i_] == '}') {
      ++i_;
      return true;
    }

    while (true) {
      Value key;
      if (!ParseString(&key, err)) return false;
      SkipWs();
      if (i_ >= s_.size() || s_[i_] != ':') { SetErr(err, "json parse: expected :"); return false; }
      ++i_;
      SkipWs();
      Value val;
      if (!ParseValue(&val, err)) return false;
      out->obj.emplace(std::move(key.str), std::move(val));
      SkipWs();
      if (i_ >= s_.size()) { SetErr(err, "json parse: unterminated object"); return false; }
      if (s_[i_] == ',') {
        ++i_;
        SkipWs();
        continue;
      }
      if (s_[i_] == '}') {
        ++i_;
        return true;
      }
      SetErr(err, "json parse: expected , or } in object");
      return false;
    }
  }

  std::string_view s_;
  usize i_{0};
};

inline const Value* Get(const Value& obj, std::string_view key) {
  if (!obj.IsObject()) return nullptr;
  auto it = obj.obj.find(std::string(key));
  if (it == obj.obj.end()) return nullptr;
  return &it->second;
}

inline bool GetString(const Value& v, std::string* out) {
  if (!out) return false;
  if (v.IsString()) { *out = v.str; return true; }
  return false;
}

inline bool GetBool(const Value& v, bool* out) {
  if (!out) return false;
  if (v.IsBool()) { *out = v.b; return true; }
  return false;
}

inline bool GetNumber(const Value& v, double* out) {
  if (!out) return false;
  if (v.IsNumber()) { *out = v.num; return true; }
  if (v.IsString()) {
    try {
      *out = std::stod(v.str);
      return true;
    } catch (...) {
      return false;
    }
  }
  return false;
}

inline bool GetU64(const Value& v, u64* out) {
  if (!out) return false;
  double x;
  if (GetNumber(v, &x)) {
    if (x < 0) return false;
    *out = static_cast<u64>(x);
    return true;
  }
  return false;
}

inline std::vector<double> GetDoubleList(const Value* v, const std::vector<double>& fallback) {
  if (!v) return fallback;
  if (v->IsArray()) {
    std::vector<double> out;
    out.reserve(v->arr.size());
    for (const auto& e : v->arr) {
      double x;
      if (GetNumber(e, &x)) out.push_back(x);
    }
    return out.empty() ? fallback : out;
  }
  double x;
  if (GetNumber(*v, &x)) return {x};
  return fallback;
}

inline std::vector<u64> GetU64List(const Value* v, const std::vector<u64>& fallback) {
  if (!v) return fallback;
  if (v->IsArray()) {
    std::vector<u64> out;
    out.reserve(v->arr.size());
    for (const auto& e : v->arr) {
      u64 x;
      if (GetU64(e, &x)) out.push_back(x);
    }
    return out.empty() ? fallback : out;
  }
  u64 x;
  if (GetU64(*v, &x)) return {x};
  return fallback;
}

inline std::vector<std::string> GetStringList(const Value* v, const std::vector<std::string>& fallback) {
  if (!v) return fallback;
  if (v->IsArray()) {
    std::vector<std::string> out;
    out.reserve(v->arr.size());
    for (const auto& e : v->arr) {
      if (e.IsString()) out.push_back(e.str);
    }
    return out.empty() ? fallback : out;
  }
  if (v->IsString()) return {v->str};
  return fallback;
}

inline std::unordered_map<std::string, std::string> GetStringMap(const Value* v) {
  std::unordered_map<std::string, std::string> out;
  if (!v || !v->IsObject()) return out;
  out.reserve(v->obj.size());
  for (const auto& kv : v->obj) {
    const Value& vv = kv.second;
    if (vv.IsString()) out.emplace(kv.first, vv.str);
    else if (vv.IsNumber()) out.emplace(kv.first, std::to_string(vv.num));
    else if (vv.IsBool()) out.emplace(kv.first, vv.b ? "true" : "false");
    else {
      // ignore arrays/objects/null
    }
  }
  return out;
}

}  // namespace json

// --------------------------
// Sweep spec
// --------------------------

struct SweepSpec {
  sjs::Config base;

  // sweep lists
  std::vector<double> alphas;
  std::vector<u64> n_r_list;
  std::vector<u64> n_s_list;
  std::vector<u64> gen_seeds;
  std::vector<u64> t_list;
  std::vector<sjs::Method> methods;
  std::vector<sjs::Variant> variants;
  std::vector<u64> seeds;  // optional explicit list (overrides repeats)

  // output
  std::string raw_file = "sweep_raw.csv";
  std::string summary_file = "sweep_summary.csv";
};

inline bool ApplyBaseConfigFromJson(const json::Value& root, sjs::Config* cfg, std::string* err) {
  if (!cfg) {
    SetErr(err, "ApplyBaseConfigFromJson: cfg is null");
    return false;
  }

  const json::Value* base = json::Get(root, "base");
  if (!base || !base->IsObject()) {
    // No base section -> keep defaults.
    return true;
  }

  // dataset
  if (const json::Value* ds = json::Get(*base, "dataset")) {
    if (ds->IsObject()) {
      if (const json::Value* v = json::Get(*ds, "source")) {
        std::string s;
        if (json::GetString(*v, &s)) {
          sjs::DataSource src;
          if (sjs::ParseDataSource(s, &src)) cfg->dataset.source = src;
        }
      }
      if (const json::Value* v = json::Get(*ds, "name")) {
        std::string s;
        if (json::GetString(*v, &s)) cfg->dataset.name = s;
      }
      if (const json::Value* v = json::Get(*ds, "dim")) {
        sjs::i32 dim;
        if (json::GetU64(*v, reinterpret_cast<u64*>(&dim))) cfg->dataset.dim = dim;
      }
      if (const json::Value* v = json::Get(*ds, "path_r")) {
        std::string s;
        if (json::GetString(*v, &s)) cfg->dataset.path_r = s;
      }
      if (const json::Value* v = json::Get(*ds, "path_s")) {
        std::string s;
        if (json::GetString(*v, &s)) cfg->dataset.path_s = s;
      }

      if (const json::Value* syn = json::Get(*ds, "synthetic")) {
        if (syn->IsObject()) {
          if (const json::Value* v = json::Get(*syn, "generator")) {
            std::string s;
            if (json::GetString(*v, &s)) cfg->dataset.synthetic.generator = s;
          }
          if (const json::Value* v = json::Get(*syn, "n_r")) {
            u64 x;
            if (json::GetU64(*v, &x)) cfg->dataset.synthetic.n_r = x;
          }
          if (const json::Value* v = json::Get(*syn, "n_s")) {
            u64 x;
            if (json::GetU64(*v, &x)) cfg->dataset.synthetic.n_s = x;
          }
          if (const json::Value* v = json::Get(*syn, "alpha")) {
            double x;
            if (json::GetNumber(*v, &x)) cfg->dataset.synthetic.alpha = x;
          }
          if (const json::Value* v = json::Get(*syn, "seed")) {
            u64 x;
            if (json::GetU64(*v, &x)) cfg->dataset.synthetic.seed = x;
          }
          if (const json::Value* v = json::Get(*syn, "params")) {
            cfg->dataset.synthetic.extra = json::GetStringMap(v);
          }
        }
      }
    }
  }

  // run
  if (const json::Value* run = json::Get(*base, "run")) {
    if (run->IsObject()) {
      if (const json::Value* v = json::Get(*run, "method")) {
        std::string s;
        if (json::GetString(*v, &s)) {
          sjs::Method m;
          if (sjs::ParseMethod(s, &m)) cfg->run.method = m;
        }
      }
      if (const json::Value* v = json::Get(*run, "variant")) {
        std::string s;
        if (json::GetString(*v, &s)) {
          sjs::Variant vv;
          if (sjs::ParseVariant(s, &vv)) cfg->run.variant = vv;
        }
      }
      if (const json::Value* v = json::Get(*run, "t")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->run.t = x;
      }
      if (const json::Value* v = json::Get(*run, "seed")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->run.seed = x;
      }
      if (const json::Value* v = json::Get(*run, "repeats")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->run.repeats = x;
      }
      if (const json::Value* v = json::Get(*run, "enum_cap")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->run.enum_cap = x;
      }
      if (const json::Value* v = json::Get(*run, "j_star")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->run.j_star = x;
      }
      if (const json::Value* v = json::Get(*run, "write_samples")) {
        bool b;
        if (json::GetBool(*v, &b)) cfg->run.write_samples = b;
      }
      if (const json::Value* v = json::Get(*run, "extra")) {
        cfg->run.extra = json::GetStringMap(v);
      }
    }
  }

  // output
  if (const json::Value* out = json::Get(*base, "output")) {
    if (out->IsObject()) {
      if (const json::Value* v = json::Get(*out, "out_dir")) {
        std::string s;
        if (json::GetString(*v, &s)) cfg->output.out_dir = s;
      }
      if (const json::Value* v = json::Get(*out, "run_tag")) {
        std::string s;
        if (json::GetString(*v, &s)) cfg->output.run_tag = s;
      }
    }
  }

  // logging
  if (const json::Value* log = json::Get(*base, "logging")) {
    if (log->IsObject()) {
      if (const json::Value* v = json::Get(*log, "level")) {
        std::string s;
        if (json::GetString(*v, &s)) {
          sjs::LogLevel lvl;
          if (sjs::detail::ParseLogLevel(s, &lvl)) cfg->logging.level = lvl;
        }
      }
      if (const json::Value* v = json::Get(*log, "with_timestamp")) {
        bool b;
        if (json::GetBool(*v, &b)) cfg->logging.with_timestamp = b;
      }
      if (const json::Value* v = json::Get(*log, "with_thread_id")) {
        bool b;
        if (json::GetBool(*v, &b)) cfg->logging.with_thread_id = b;
      }
    }
  }

  // sys
  if (const json::Value* sys = json::Get(*base, "sys")) {
    if (sys->IsObject()) {
      if (const json::Value* v = json::Get(*sys, "threads")) {
        u64 x;
        if (json::GetU64(*v, &x)) cfg->sys.threads = static_cast<sjs::i32>(x);
      }
    }
  }

  return true;
}

inline bool LoadSweepSpecJson(const std::string& path, SweepSpec* spec, std::string* err) {
  if (!spec) {
    SetErr(err, "LoadSweepSpecJson: spec is null");
    return false;
  }

  std::ifstream in(path);
  if (!in) {
    SetErr(err, "Cannot open sweep config: " + path);
    return false;
  }
  std::stringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();

  json::Value root;
  std::string jerr;
  if (!json::Parser(std::string_view(text)).Parse(&root, &jerr)) {
    SetErr(err, "JSON parse failed: " + jerr);
    return false;
  }
  if (!root.IsObject()) {
    SetErr(err, "Sweep JSON must be a top-level object");
    return false;
  }

  // Start from defaults (caller should fill spec->base before calling if desired).
  // Apply base config overrides.
  if (!ApplyBaseConfigFromJson(root, &spec->base, err)) return false;

  // sweep lists
  const json::Value* sweep = json::Get(root, "sweep");
  if (sweep && sweep->IsObject()) {
    spec->alphas = json::GetDoubleList(json::Get(*sweep, "alpha"), spec->alphas);
    spec->n_r_list = json::GetU64List(json::Get(*sweep, "n_r"), spec->n_r_list);
    spec->n_s_list = json::GetU64List(json::Get(*sweep, "n_s"), spec->n_s_list);
    spec->gen_seeds = json::GetU64List(json::Get(*sweep, "gen_seed"), spec->gen_seeds);
    spec->t_list = json::GetU64List(json::Get(*sweep, "t"), spec->t_list);

    // method/variant as strings -> enums
    if (const json::Value* m = json::Get(*sweep, "method")) {
      std::vector<std::string> ms = json::GetStringList(m, {});
      if (!ms.empty()) {
        spec->methods.clear();
        for (const auto& s : ms) {
          sjs::Method mm;
          if (sjs::ParseMethod(s, &mm) && mm != sjs::Method::Unknown) {
            spec->methods.push_back(mm);
          }
        }
      }
    }
    if (const json::Value* v = json::Get(*sweep, "variant")) {
      std::vector<std::string> vs = json::GetStringList(v, {});
      if (!vs.empty()) {
        spec->variants.clear();
        for (const auto& s : vs) {
          sjs::Variant vv;
          if (sjs::ParseVariant(s, &vv)) spec->variants.push_back(vv);
        }
      }
    }

    spec->seeds = json::GetU64List(json::Get(*sweep, "seed"), spec->seeds);

    // Optional override for repeats
    if (const json::Value* r = json::Get(*sweep, "repeats")) {
      u64 x;
      if (json::GetU64(*r, &x)) spec->base.run.repeats = x;
    }
  }

  // files
  if (const json::Value* files = json::Get(root, "files")) {
    if (files->IsObject()) {
      if (const json::Value* v = json::Get(*files, "raw")) {
        std::string s;
        if (json::GetString(*v, &s)) spec->raw_file = s;
      }
      if (const json::Value* v = json::Get(*files, "summary")) {
        std::string s;
        if (json::GetString(*v, &s)) spec->summary_file = s;
      }
    }
  }

  // Fill default lists from base config if still empty.
  if (spec->alphas.empty()) spec->alphas = {spec->base.dataset.synthetic.alpha};
  if (spec->n_r_list.empty()) spec->n_r_list = {spec->base.dataset.synthetic.n_r};
  if (spec->n_s_list.empty()) spec->n_s_list = {spec->base.dataset.synthetic.n_s};
  if (spec->gen_seeds.empty()) spec->gen_seeds = {spec->base.dataset.synthetic.seed};
  if (spec->t_list.empty()) spec->t_list = {spec->base.run.t};
  if (spec->methods.empty()) spec->methods = {spec->base.run.method};
  if (spec->variants.empty()) spec->variants = {spec->base.run.variant};

  return true;
}

inline void PrintUsage() {
  std::cerr
      << "sjs_sweep: run a parameter sweep from a JSON config\n\n"
      << "Required:\n"
      << "  --config=<path/to/sweep.json>    (or provide JSON path as a positional arg)\n\n"
      << "Common overrides (optional):\n"
      << "  --out_dir=<dir>                 (overrides base.output.out_dir)\n"
      << "  --raw_file=<path>               (overrides files.raw)\n"
      << "  --summary_file=<path>           (overrides files.summary)\n\n"
      << "Baselines supported in this build (Dim=2):\n"
      << sjs::baselines::BaselineHelp2D()
      << "\n";
}

}  // namespace

}  // namespace apps
}  // namespace sjs

int main(int argc, char** argv) {
  using sjs::u64;
  using sjs::usize;

  sjs::ArgMap args = sjs::ArgMap::FromArgv(argc, argv);
  if (sjs::apps::IsHelpRequested(args)) {
    sjs::apps::PrintUsage();
    return 0;
  }

  // Locate sweep JSON path.
  std::string sweep_path;
  if (auto v = args.Get("config")) {
    sweep_path = std::string(*v);
  } else if (!args.Positional().empty()) {
    sweep_path = args.Positional().front();
  }
  if (sweep_path.empty()) {
    SJS_LOG_ERROR("Missing --config=<sweep.json>.");
    sjs::apps::PrintUsage();
    return 2;
  }

  // Start from CLI config (so CLI can provide quick overrides), then JSON overrides base.
  sjs::Config cfg0 = sjs::Config::FromArgs(argc, argv);

  // Load sweep spec.
  sjs::apps::SweepSpec spec;
  spec.base = cfg0;

  std::string err;
  if (!sjs::apps::LoadSweepSpecJson(sweep_path, &spec, &err)) {
    SJS_LOG_ERROR(err);
    return 3;
  }

  // Apply a few CLI overrides on top (out_dir/raw_file/summary_file).
  if (auto v = args.Get("out_dir")) spec.base.output.out_dir = std::string(*v);
  if (auto v = args.Get("raw_file")) spec.raw_file = std::string(*v);
  if (auto v = args.Get("summary_file")) spec.summary_file = std::string(*v);

  if (!spec.base.Validate(&err)) {
    SJS_LOG_ERROR("Base config validation failed: ", err);
    return 4;
  }
  sjs::Logger::Instance().SetConfig(spec.base.logging);

  if (spec.base.dataset.dim != 2) {
    SJS_LOG_ERROR("This build currently supports only Dim=2. Got --dim=", spec.base.dataset.dim);
    return 4;
  }

  // Output files.
  const fs::path out_dir(spec.base.output.out_dir);
  if (!sjs::apps::EnsureDir(out_dir, &err)) {
    SJS_LOG_ERROR("Cannot create out_dir: ", err);
    return 5;
  }

  const std::string raw_path = (out_dir / spec.raw_file).string();
  const std::string summary_path = (out_dir / spec.summary_file).string();

  if (!sjs::apps::EnsureDir(fs::path(raw_path).parent_path(), &err) ||
      !sjs::apps::EnsureDir(fs::path(summary_path).parent_path(), &err)) {
    SJS_LOG_ERROR("Cannot create output directories: ", err);
    return 5;
  }

  sjs::csv::Writer raw_writer(raw_path, sjs::csv::Dialect{','}, &err);
  if (!raw_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open raw CSV: ", raw_path, " err=", err);
    return 5;
  }
  sjs::csv::Writer sum_writer(summary_path, sjs::csv::Dialect{','}, &err);
  if (!sum_writer.Ok()) {
    SJS_LOG_ERROR("Cannot open summary CSV: ", summary_path, " err=", err);
    return 5;
  }

  raw_writer.WriteHeader(sjs::apps::RawHeader(), &err);
  sum_writer.WriteHeader(sjs::apps::SummaryHeader(), &err);

  // For binary/csv, load dataset once (no dataset-param sweep).
  sjs::Dataset<2, sjs::Scalar> fixed_ds;
  sjs::synthetic::Report fixed_gen_report;
  bool has_fixed_ds = false;
  if (spec.base.dataset.source != sjs::DataSource::Synthetic) {
    std::string local_err;
    if (!sjs::apps::LoadOrGenerateDataset<2>(spec.base, &fixed_ds, &fixed_gen_report, &local_err)) {
      SJS_LOG_ERROR(local_err);
      return 6;
    }
    has_fixed_ds = true;
    SJS_LOG_INFO("Loaded fixed dataset: ", fixed_ds.name,
                 " R=", static_cast<unsigned long long>(fixed_ds.R.Size()),
                 " S=", static_cast<unsigned long long>(fixed_ds.S.Size()));
  }

  // --------------------------
  // Sweep loops
  // --------------------------
  u64 total_runs = 0;
  u64 ok_runs = 0;

  const bool sweep_dataset = (spec.base.dataset.source == sjs::DataSource::Synthetic);

  // If dataset is not synthetic, we treat n_r/n_s/alpha/gen_seed lists as singletons.
  const std::vector<u64> n_r_list = sweep_dataset ? spec.n_r_list : std::vector<u64>{0};
  const std::vector<u64> n_s_list = sweep_dataset ? spec.n_s_list : std::vector<u64>{0};
  const std::vector<double> alpha_list = sweep_dataset ? spec.alphas : std::vector<double>{0.0};
  const std::vector<u64> gen_seeds = sweep_dataset ? spec.gen_seeds : std::vector<u64>{0};

  for (u64 n_r : n_r_list) {
    for (u64 n_s : n_s_list) {
      for (double alpha : alpha_list) {
        for (u64 gen_seed : gen_seeds) {
          // Generate/load dataset for this combo.
          sjs::Dataset<2, sjs::Scalar> ds;
          sjs::synthetic::Report gen_report;
          const sjs::synthetic::Report* gen_report_ptr = nullptr;

          sjs::Config cfg_ds = spec.base;

          if (sweep_dataset) {
            cfg_ds.dataset.synthetic.n_r = n_r;
            cfg_ds.dataset.synthetic.n_s = n_s;
            cfg_ds.dataset.synthetic.alpha = alpha;
            cfg_ds.dataset.synthetic.seed = gen_seed;

            // Auto-name dataset so rows remain self-describing.
            {
              std::ostringstream nm;
              nm << (spec.base.dataset.name.empty() ? "synthetic" : spec.base.dataset.name)
                 << "__nr" << static_cast<unsigned long long>(n_r)
                 << "__ns" << static_cast<unsigned long long>(n_s)
                 << "__a" << sjs::apps::FormatSci(alpha, 3)
                 << "__g" << static_cast<unsigned long long>(gen_seed);
              cfg_ds.dataset.name = sjs::apps::SanitizeFilename(nm.str());
            }

            std::string local_err;
            if (!sjs::apps::LoadOrGenerateDataset<2>(cfg_ds, &ds, &gen_report, &local_err)) {
              SJS_LOG_ERROR("Dataset generation failed: ", local_err);
              return 6;
            }
            gen_report_ptr = &gen_report;
          } else {
            (void)n_r; (void)n_s; (void)alpha; (void)gen_seed;
            if (!has_fixed_ds) {
              SJS_LOG_ERROR("Internal error: expected fixed dataset to be loaded.");
              return 6;
            }
            ds = fixed_ds;
            gen_report_ptr = nullptr;
          }

          for (u64 t : spec.t_list) {
            for (sjs::Method method : spec.methods) {
              for (sjs::Variant variant : spec.variants) {
                // Prepare cfg for this group.
                sjs::Config cfg = cfg_ds;
                cfg.run.t = t;
                cfg.run.method = method;
                cfg.run.variant = variant;

                // Determine seeds.
                std::vector<u64> seeds;
                if (!spec.seeds.empty()) {
                  seeds = spec.seeds;
                } else {
                  seeds.resize(static_cast<usize>(cfg.run.repeats));
                  for (u64 rep = 0; rep < cfg.run.repeats; ++rep) {
                    seeds[static_cast<usize>(rep)] = cfg.run.seed + rep;
                  }
                }

                // Create baseline.
                std::string b_err;
                auto baseline = sjs::baselines::CreateBaseline2D(method, variant, &b_err);
                if (!baseline) {
                  SJS_LOG_ERROR("CreateBaseline2D failed: ", b_err);
                  return 7;
                }

                // Collect per-repeat stats.
                std::vector<double> wall_ms;
                wall_ms.reserve(seeds.size());
                std::vector<double> count_vals;
                count_vals.reserve(seeds.size());
                u64 exact_count = 0;
                u64 ok_in_group = 0;

                for (usize i = 0; i < seeds.size(); ++i) {
                  const u64 seed = seeds[i];
                  sjs::baselines::RunReport rep_out;
                  std::string local_err;

                  sjs::Stopwatch sw;
                  bool ok = false;
                  switch (variant) {
                    case sjs::Variant::Sampling:
                      ok = sjs::baselines::RunSamplingOnce<2, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                    case sjs::Variant::EnumSampling:
                      ok = sjs::baselines::RunEnumSamplingOnce<2, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                    case sjs::Variant::Adaptive:
                      ok = sjs::baselines::RunAdaptiveOnce<2, sjs::Scalar>(
                          baseline.get(), ds, cfg, seed, &rep_out, &local_err);
                      break;
                  }

                  const double wms = sw.ElapsedMillis();
                  wall_ms.push_back(wms);
                  count_vals.push_back(static_cast<double>(rep_out.count.value));

                  total_runs++;
                  if (!ok) {
                    rep_out.ok = false;
                    rep_out.error = local_err;
                  } else {
                    ok_runs++;
                    ok_in_group++;
                    if (rep_out.count.exact) exact_count++;
                  }

                  // Raw row.
                  raw_writer.WriteRowV(
                      ds.name,
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                      (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha
                                                                        : std::numeric_limits<double>::quiet_NaN()),
                      static_cast<unsigned long long>(ds.R.Size()),
                      static_cast<unsigned long long>(ds.S.Size()),
                      sjs::ToString(method),
                      sjs::ToString(variant),
                      static_cast<unsigned long long>(t),
                      static_cast<unsigned long long>(i),
                      static_cast<unsigned long long>(seed),
                      (rep_out.ok ? 1 : 0),
                      wms,
                      static_cast<double>(rep_out.count.value),
                      (rep_out.count.exact ? 1 : 0),
                      (rep_out.used_enumeration ? 1 : 0),
                      (rep_out.enumeration_truncated ? 1 : 0),
                      static_cast<unsigned long long>(rep_out.enumeration_cap),
                      rep_out.adaptive_branch,
                      static_cast<unsigned long long>(rep_out.adaptive_pilot_pairs),
                      rep_out.phases.ToJsonMillis(),
                      rep_out.note,
                      cfg.ToJsonLite(),
                      (gen_report_ptr ? gen_report_ptr->ToJsonLite() : "{}"));
                }

                // Summary row.
                const sjs::Summary wall_sum = sjs::Summarize(wall_ms);
                const sjs::Summary cnt_sum = sjs::Summarize(count_vals);
                const double ok_rate = (seeds.empty() ? 0.0 : static_cast<double>(ok_in_group) /
                                                          static_cast<double>(seeds.size()));
                const double exact_frac = (seeds.empty() ? 0.0 : static_cast<double>(exact_count) /
                                                          static_cast<double>(seeds.size()));

                sum_writer.WriteRowV(
                    ds.name,
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.generator : ""),
                    (cfg.dataset.source == sjs::DataSource::Synthetic ? cfg.dataset.synthetic.alpha
                                                                      : std::numeric_limits<double>::quiet_NaN()),
                    static_cast<unsigned long long>(ds.R.Size()),
                    static_cast<unsigned long long>(ds.S.Size()),
                    sjs::ToString(method),
                    sjs::ToString(variant),
                    static_cast<unsigned long long>(t),
                    static_cast<unsigned long long>(seeds.size()),
                    ok_rate,
                    wall_sum.mean,
                    wall_sum.stdev,
                    wall_sum.median,
                    wall_sum.p95,
                    cnt_sum.mean,
                    cnt_sum.stdev,
                    exact_frac,
                    "");
              }
            }
          }
        }
      }
    }
  }

  SJS_LOG_INFO("Sweep finished. total_runs=", static_cast<unsigned long long>(total_runs),
               " ok_runs=", static_cast<unsigned long long>(ok_runs));
  SJS_LOG_INFO("Raw CSV: ", raw_path);
  SJS_LOG_INFO("Summary CSV: ", summary_path);

  // Copy sweep config next to outputs for reproducibility.
  try {
    const fs::path dst = out_dir / "sweep_config.json";
    std::error_code ec;
    fs::copy_file(fs::path(sweep_path), dst, fs::copy_options::overwrite_existing, ec);
    if (!ec) {
      SJS_LOG_INFO("Saved sweep config to: ", dst.string());
    }
  } catch (...) {
    // ignore
  }

  return 0;
}

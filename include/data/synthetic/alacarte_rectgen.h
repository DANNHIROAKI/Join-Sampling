#pragma once
// data/synthetic/alacarte_rectgen.h
//
// Synthetic dataset generator backed by the local Alacarte source tree
// (third_party/Alacarte/alacarte_rectgen.py) in this repository.
//
// Motivation
// ----------
// Your project needs a dimension-agnostic synthetic generator with a knob for
// output density:
//   alpha_out = |J(R,S)| / (|R| + |S|)
//
// The local Alacarte module provides:
//   R, S, info = alacarte_rectgen.make_rectangles_R_S(...)
// and supports tuning expected alpha_out by adjusting a coverage parameter.
//
// Design
// ------
// We do not embed Python. Instead, we spawn the helper script:
//   tools/alacarte_rectgen_generate.py
// which writes SJSBOX-v1 binary relations for R and S, and a JSON report.
// Then we load those binaries back into C++ using io/binary_io.h.
//
// This makes it easy to keep the synthetic method consistent with the upstream
// Python implementation and avoids adding a CPython link dependency.
//
// Configuration (spec.params)
// ---------------------------
// Optional stringly-typed keys:
//   - python          : Python executable (default: env SJS_PYTHON or "python3")
//   - rectgen_script  : Path to tools/alacarte_rectgen_generate.py
//                      (default: env SJS_ALACARTE_RECTGEN_SCRIPT or
//                       "tools/alacarte_rectgen_generate.py")
//   - alacarte_module : Path to third_party/Alacarte/alacarte_rectgen.py
//                      (default: env SJS_ALACARTE_MODULE or
//                       resolved from rectgen_script location)
//   - audit_pairs     : u64 number of pairs for pair-sampling audit (default 2,000,000)
//                       set to 0 to disable audit
//   - audit_seed      : u64 seed for audit sampling (default 0)
//   - temp_dir        : Override base temp directory for intermediate binaries
//                      (default: env SJS_TEMP_DIR or system temp)
//   - keep_files      : bool ("0/1", "true/false") keep temp files (default false)

#include "data/synthetic/generator.h"

#include "io/binary_io.h"

#include <cerrno>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace sjs {
namespace synthetic {

namespace detail {

inline std::string ShellQuote(std::string_view s) {
#if defined(_WIN32)
  // Minimal Windows quoting.
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('"');
  for (char c : s) {
    if (c == '"') out += "\\\"";
    else out.push_back(c);
  }
  out.push_back('"');
  return out;
#else
  // POSIX sh-safe single-quote escaping: ' -> '\''
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('\'');
  for (char c : s) {
    if (c == '\'') out += "'\\''";
    else out.push_back(c);
  }
  out.push_back('\'');
  return out;
#endif
}

inline bool ParseBool(std::string_view s, bool def) {
  using detail::EqualsIgnoreCase;
  if (EqualsIgnoreCase(s, "1") || EqualsIgnoreCase(s, "true") || EqualsIgnoreCase(s, "yes") ||
      EqualsIgnoreCase(s, "y") || EqualsIgnoreCase(s, "on")) {
    return true;
  }
  if (EqualsIgnoreCase(s, "0") || EqualsIgnoreCase(s, "false") || EqualsIgnoreCase(s, "no") ||
      EqualsIgnoreCase(s, "n") || EqualsIgnoreCase(s, "off")) {
    return false;
  }
  return def;
}

inline std::string GetParamOrEnv(const std::unordered_map<std::string, std::string>& params,
                                 std::string_view key,
                                 const char* env_key,
                                 std::string def) {
  if (auto it = params.find(std::string(key)); it != params.end()) return it->second;
  if (env_key) {
    if (const char* v = std::getenv(env_key); v && *v) return std::string(v);
  }
  return def;
}

inline u64 GetParamU64(const std::unordered_map<std::string, std::string>& params,
                       std::string_view key,
                       u64 def) {
  auto it = params.find(std::string(key));
  if (it == params.end()) return def;
  u64 v = 0;
  if (!TryParseU64(it->second, &v)) return def;
  return v;
}

inline double GetParamDouble(const std::unordered_map<std::string, std::string>& params,
                           std::string_view key,
                           double def) {
  auto it = params.find(std::string(key));
  if (it == params.end()) return def;
  double v = def;
  if (!TryParseDouble(it->second, &v)) return def;
  return v;
}

inline bool GetParamBool(const std::unordered_map<std::string, std::string>& params,
                         std::string_view key,
                         bool def) {
  auto it = params.find(std::string(key));
  if (it == params.end()) return def;
  return ParseBool(it->second, def);
}

inline std::string ReadFileToString(const std::filesystem::path& p, std::string* err) {
  std::ifstream in(p, std::ios::binary);
  if (!in) {
    if (err) *err = "Cannot open file: " + p.string();
    return {};
  }
  std::ostringstream oss;
  oss << in.rdbuf();
  return oss.str();
}

// Extremely small helper: extract a top-level JSON number field by key.
// Accepts either:
//   "key": 1.23
//   "key": null
// Returns true if a numeric value was parsed.
inline bool ExtractJsonNumber(const std::string& json,
                              std::string_view key,
                              double* out) {
  if (!out) return false;
  const std::string pat = "\"" + std::string(key) + "\"";
  std::size_t pos = json.find(pat);
  if (pos == std::string::npos) return false;
  pos = json.find(':', pos + pat.size());
  if (pos == std::string::npos) return false;
  ++pos;
  while (pos < json.size() && std::isspace(static_cast<unsigned char>(json[pos]))) ++pos;
  if (pos + 4 <= json.size() && json.compare(pos, 4, "null") == 0) return false;
  const char* begin = json.c_str() + pos;
  char* end = nullptr;
  errno = 0;
  const double v = std::strtod(begin, &end);
  if (end == begin || errno != 0) return false;
  *out = v;
  return true;
}

inline std::filesystem::path MakeTempDir(u64 seed,
                                        const std::string& base_override,
                                        std::string* err) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::path base;
  if (!base_override.empty()) {
    base = fs::path(base_override);
    if (!fs::exists(base, ec)) {
      if (err) *err = "Temp base dir does not exist: " + base.string();
      return {};
    }
  } else {
    base = fs::temp_directory_path(ec);
    if (ec) {
      if (err) *err = "temp_directory_path failed: " + ec.message();
      return {};
    }
  }

  const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::ostringstream tag;
  tag << "sjs_rectgen_" << seed << "_" << now;
  fs::path dir = base / tag.str();

  if (!fs::create_directories(dir, ec) || ec) {
    if (err) *err = "Failed to create temp dir: " + dir.string() + " (" + ec.message() + ")";
    return {};
  }
  return dir;
}

}  // namespace detail

// ----------------------------
// Generator
// ----------------------------

template <int Dim, class T = Scalar>
class AlacarteRectGenGenerator final : public ISyntheticGenerator<Dim, T> {
 public:
  using DatasetT = Dataset<Dim, T>;

  std::string_view Name() const noexcept override { return "alacarte_rectgen"; }

  bool Generate(const DatasetSpec& spec, DatasetT* out_ds, Report* report, std::string* err) override {
    namespace fs = std::filesystem;

    if (!out_ds) {
      detail::SetErr(err, "AlacarteRectGenGenerator: out_ds is null");
      return false;
    }
    if (spec.n_r == 0 || spec.n_s == 0) {
      detail::SetErr(err, "AlacarteRectGenGenerator: n_r and n_s must be > 0");
      return false;
    }
    if (!detail::CheckIdFits(spec.n_r, err, "AlacarteRectGenGenerator(R)")) return false;
    if (!detail::CheckIdFits(spec.n_s, err, "AlacarteRectGenGenerator(S)")) return false;

    // Resolve python + script.
    const std::string python = detail::GetParamOrEnv(spec.params, "python", "SJS_PYTHON", "python3");
    std::string script = detail::GetParamOrEnv(spec.params,
                                               "rectgen_script",
                                               "SJS_ALACARTE_RECTGEN_SCRIPT",
                                               "tools/alacarte_rectgen_generate.py");
    if (!fs::exists(script)) {
      // Common fallback if running from build dir.
      const fs::path p0 = fs::path("tools") / "alacarte_rectgen_generate.py";
      const fs::path p1 = fs::path("../tools") / "alacarte_rectgen_generate.py";
      const fs::path p2 = fs::path("../../tools") / "alacarte_rectgen_generate.py";
      const fs::path p3 = fs::path("../../../tools") / "alacarte_rectgen_generate.py";
      if (fs::exists(p0)) script = p0.string();
      else if (fs::exists(p1)) script = p1.string();
      else if (fs::exists(p2)) script = p2.string();
      else if (fs::exists(p3)) script = p3.string();
    }
    if (!fs::exists(script)) {
      detail::SetErr(err,
                     "AlacarteRectGenGenerator: rectgen_script not found: " + script +
                     " (set spec.params['rectgen_script'] or env SJS_ALACARTE_RECTGEN_SCRIPT)");
      return false;
    }

    const u64 audit_pairs = detail::GetParamU64(spec.params, "audit_pairs", 2'000'000ULL);
    const u64 audit_seed = detail::GetParamU64(spec.params, "audit_seed", 0ULL);
    const bool keep_files = detail::GetParamBool(spec.params, "keep_files", false);
    const std::string alacarte_module = detail::GetParamOrEnv(
        spec.params, "alacarte_module", "SJS_ALACARTE_MODULE", "");
    const std::string volume_dist = detail::GetParamOrEnv(spec.params, "volume_dist", nullptr, "");
    const double volume_cv = detail::GetParamDouble(spec.params, "volume_cv", std::numeric_limits<double>::quiet_NaN());
    const double shape_sigma = detail::GetParamDouble(spec.params, "shape_sigma", std::numeric_limits<double>::quiet_NaN());
    const u64 tune_samples = detail::GetParamU64(spec.params, "tune_samples", 0ULL);
    const double tune_tol_rel = detail::GetParamDouble(spec.params, "tune_tol_rel", std::numeric_limits<double>::quiet_NaN());
    const u64 tune_max_iter = detail::GetParamU64(spec.params, "tune_max_iter", 0ULL);
    const std::string dtype = detail::GetParamOrEnv(spec.params, "dtype", nullptr, "");

    // Temp output.
    // You may override the base temp directory via:
    //   spec.params["temp_dir"] or env SJS_TEMP_DIR.
    const std::string temp_base = detail::GetParamOrEnv(spec.params, "temp_dir", "SJS_TEMP_DIR", "");
    std::string terr;
    fs::path tmp_dir = detail::MakeTempDir(spec.seed, temp_base, &terr);
    if (tmp_dir.empty()) {
      detail::SetErr(err, "AlacarteRectGenGenerator: " + terr);
      return false;
    }

    const fs::path out_r = tmp_dir / "R.bin";
    const fs::path out_s = tmp_dir / "S.bin";
    const fs::path rep_path = tmp_dir / "gen_report.json";

    // Build command.
    std::ostringstream cmd;
    cmd << detail::ShellQuote(python)
        << " " << detail::ShellQuote(script)
        << " --nR=" << spec.n_r
        << " --nS=" << spec.n_s
        << " --d=" << Dim
        << " --alpha_out=" << spec.alpha
        << " --seed=" << spec.seed
        << " --out_r=" << detail::ShellQuote(out_r.string())
        << " --out_s=" << detail::ShellQuote(out_s.string())
        << " --dataset_name=" << detail::ShellQuote(spec.name)
        << " --report_path=" << detail::ShellQuote(rep_path.string())
        << " --audit_pairs=" << audit_pairs
        << " --audit_seed=" << audit_seed
        << " --domain_lo=" << spec.domain_lo
        << " --domain_hi=" << spec.domain_hi;
    if (!alacarte_module.empty()) {
      cmd << " --alacarte_module=" << detail::ShellQuote(alacarte_module);
    }
    if (!volume_dist.empty()) {
      cmd << " --volume_dist=" << detail::ShellQuote(volume_dist);
    }
    if (std::isfinite(volume_cv)) {
      cmd << " --volume_cv=" << volume_cv;
    }
    if (std::isfinite(shape_sigma)) {
      cmd << " --shape_sigma=" << shape_sigma;
    }
    if (tune_samples > 0) {
      cmd << " --tune_samples=" << tune_samples;
    }
    if (std::isfinite(tune_tol_rel)) {
      cmd << " --tune_tol_rel=" << tune_tol_rel;
    }
    if (tune_max_iter > 0) {
      cmd << " --tune_max_iter=" << tune_max_iter;
    }
    if (!dtype.empty()) {
      cmd << " --dtype=" << detail::ShellQuote(dtype);
    }

    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
      detail::SetErr(err,
                     "AlacarteRectGenGenerator: python generator failed rc=" + std::to_string(rc) +
                     ". Ensure local Alacarte source is available (third_party/Alacarte/alacarte_rectgen.py) "
                     "or set spec.params['alacarte_module'] / env SJS_ALACARTE_MODULE.");
      if (!keep_files) {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
      }
      return false;
    }

    if (!fs::exists(out_r) || !fs::exists(out_s)) {
      detail::SetErr(err,
                     "AlacarteRectGenGenerator: expected output files missing under " + tmp_dir.string());
      if (!keep_files) {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
      }
      return false;
    }

    // Load binary into memory.
    sjs::binary::BinaryReadOptions ropt;
    ropt.generate_ids_if_missing = true;
    ropt.drop_empty = false;
    if (!sjs::binary::ReadDatasetBinaryPair<Dim, T>(out_r.string(), out_s.string(), out_ds, ropt, err)) {
      if (!keep_files) {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
      }
      return false;
    }

    // Fill report.
    if (report) {
      report->generator = std::string(Name());
      report->dataset_name = spec.name;
      report->n_r = spec.n_r;
      report->n_s = spec.n_s;
      report->has_exact_k = false;
      report->k_target = 0;
      report->k_achieved = 0;
      report->alpha_target = spec.alpha;
      report->alpha_achieved = 0.0;

      std::string rep_json = detail::ReadFileToString(rep_path, nullptr);
      double alpha_expected_est = 0.0;
      double coverage = 0.0;
      double p_est = 0.0;
      double alpha_hat_est = 0.0;
      double eps_alpha = 0.0;
      const bool has_alpha_expected = detail::ExtractJsonNumber(rep_json, "alpha_expected_est", &alpha_expected_est);
      const bool has_coverage = detail::ExtractJsonNumber(rep_json, "coverage", &coverage);
      const bool has_p_est = detail::ExtractJsonNumber(rep_json, "pair_intersection_prob_est", &p_est);
      const bool has_alpha_hat = detail::ExtractJsonNumber(rep_json, "alpha_hat_est", &alpha_hat_est);
      const bool has_eps = detail::ExtractJsonNumber(rep_json, "epsilon_alpha", &eps_alpha);

      // Best-effort achieved alpha.
      if (has_alpha_hat) report->alpha_achieved = alpha_hat_est;
      else if (has_alpha_expected) report->alpha_achieved = alpha_expected_est;
      else report->alpha_achieved = spec.alpha;

      std::ostringstream notes;
      notes << "alacarte-local";
      if (has_alpha_expected) notes << " alpha_expected_est=" << alpha_expected_est;
      if (has_coverage) notes << " coverage=" << coverage;
      if (has_p_est) notes << " pair_p_est=" << p_est;
      if (has_alpha_hat) notes << " alpha_hat_est=" << alpha_hat_est;
      if (has_eps) notes << " epsilon_alpha=" << eps_alpha;
      notes << " audit_pairs=" << audit_pairs;
      report->notes = notes.str();
    }

    if (!keep_files) {
      std::error_code ec;
      fs::remove_all(tmp_dir, ec);
    }
    return true;
  }
};

}  // namespace synthetic
}  // namespace sjs

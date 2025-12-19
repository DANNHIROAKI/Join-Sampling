#pragma once
// sjs/io/csv_io.h
//
// CSV/TSV utilities for writing experiment outputs and (optionally) reading
// simple box datasets.
//
// Primary use in this project:
//  - Write results CSV/TSV for plotting (runtime, memory, breakdown, etc.)
//  - Optionally write sampled pairs for inspection
//
// We keep the implementation lightweight and dependency-free.
// Writer supports proper escaping for commas/quotes/newlines.
// Reader is intentionally simple and does NOT support quoted separators.
//
// Writer is safe for large outputs (streaming).
// Reader is intended for debugging / small to medium datasets.

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/io/dataset.h"

#include <cerrno>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace csv {

struct Dialect {
  char sep = ',';          // ',' for CSV, '\t' for TSV
  char quote = '"';
  bool always_quote = false;
  bool write_header = true;
};

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

// Escape a cell according to RFC 4180-ish CSV rules:
//  - If cell contains sep/quote/newline, wrap with quotes and double embedded quotes.
//  - If always_quote, always wrap and escape quotes.
inline std::string EscapeCell(std::string_view cell, const Dialect& d) {
  bool need_quote = d.always_quote;
  for (char c : cell) {
    if (c == d.sep || c == d.quote || c == '\n' || c == '\r') {
      need_quote = true;
      break;
    }
  }
  if (!need_quote) return std::string(cell);

  std::string out;
  out.reserve(cell.size() + 2);
  out.push_back(d.quote);
  for (char c : cell) {
    if (c == d.quote) out.push_back(d.quote);  // escape by doubling
    out.push_back(c);
  }
  out.push_back(d.quote);
  return out;
}

class Writer {
 public:
  Writer() = default;

  explicit Writer(const std::string& path, Dialect d = {}, std::string* err = nullptr)
      : d_(d), out_(path) {
    if (!out_) {
      SetErr(err, "Cannot open CSV for writing: " + path);
      ok_ = false;
    } else {
      ok_ = true;
    }
  }

  bool Ok() const noexcept { return ok_ && static_cast<bool>(out_); }

  bool WriteHeader(const std::vector<std::string>& cols, std::string* err = nullptr) {
    if (!Ok()) { SetErr(err, "Writer not ok"); return false; }
    if (!d_.write_header) return true;
    return WriteRow(cols, err);
  }

  bool WriteRow(const std::vector<std::string>& cols, std::string* err = nullptr) {
    if (!Ok()) { SetErr(err, "Writer not ok"); return false; }
    for (usize i = 0; i < cols.size(); ++i) {
      if (i) out_ << d_.sep;
      out_ << EscapeCell(cols[i], d_);
    }
    out_ << "\n";
    if (!out_) {
      SetErr(err, "CSV write failed");
      return false;
    }
    return true;
  }

  // Variadic row writer, converts each arg via operator<< into a string.
  template <class... Args>
  bool WriteRowV(Args&&... args) {
    if (!Ok()) return false;
    bool first = true;
    auto write_one = [&](auto&& x) {
      if (!first) out_ << d_.sep;
      first = false;
      std::ostringstream oss;
      oss << std::forward<decltype(x)>(x);
      out_ << EscapeCell(oss.str(), d_);
    };
    (write_one(std::forward<Args>(args)), ...);
    out_ << "\n";
    return static_cast<bool>(out_);
  }

 private:
  Dialect d_{};
  std::ofstream out_;
  bool ok_{false};
};

// --------------------------
// Convenience writers
// --------------------------

// Write sampled pairs as a 2-column CSV/TSV: r_id, s_id
inline bool WritePairs(const std::string& path,
                       Span<const PairId> pairs,
                       const Dialect& d = Dialect{'\t'},
                       std::string* err = nullptr) {
  Writer w(path, d, err);
  if (!w.Ok()) return false;
  if (!w.WriteHeader({"r_id", "s_id"}, err)) return false;
  for (const auto& p : pairs) {
    if (!w.WriteRowV(p.r, p.s)) {
      SetErr(err, "Failed writing pairs");
      return false;
    }
  }
  return true;
}

// Write boxes for visualization/debug:
// columns: id, lo0..lo(d-1), hi0..hi(d-1)
template <int Dim, class T>
bool WriteBoxes(const std::string& path,
                const Relation<Dim, T>& rel,
                const Dialect& d = Dialect{','},
                std::string* err = nullptr) {
  Writer w(path, d, err);
  if (!w.Ok()) return false;

  std::vector<std::string> header;
  header.reserve(1 + 2 * Dim);
  header.push_back("id");
  for (int i = 0; i < Dim; ++i) header.push_back("lo" + std::to_string(i));
  for (int i = 0; i < Dim; ++i) header.push_back("hi" + std::to_string(i));
  if (!w.WriteHeader(header, err)) return false;

  const bool has_ids = rel.HasExplicitIds();
  std::vector<std::string> row;
  row.reserve(1 + 2 * Dim);

  for (usize i = 0; i < rel.boxes.size(); ++i) {
    row.clear();
    const Id id = has_ids ? rel.ids[i] : static_cast<Id>(i);
    row.push_back(std::to_string(id));
    for (int k = 0; k < Dim; ++k) row.push_back(std::to_string(rel.boxes[i].lo.v[k]));
    for (int k = 0; k < Dim; ++k) row.push_back(std::to_string(rel.boxes[i].hi.v[k]));
    if (!w.WriteRow(row, err)) return false;
  }
  return true;
}

// --------------------------
// Simple CSV/TSV reader for boxes (optional).
// Intended for debugging, not giant production datasets.
//
// Expected format (header optional):
//   id, lo0, lo1, ..., hi0, hi1, ...
// If id column is missing, sequential ids are used.
// Lines starting with '#' are skipped.
// No support for quoted separators.
// --------------------------

inline std::vector<std::string_view> SplitLine(std::string_view line, char sep) {
  std::vector<std::string_view> out;
  usize start = 0;
  while (start <= line.size()) {
    const usize pos = line.find(sep, start);
    if (pos == std::string_view::npos) {
      out.push_back(line.substr(start));
      break;
    }
    out.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
  return out;
}

inline std::string_view Trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return s;
}

inline bool ParseDouble(std::string_view s, double* out) {
  if (!out) return false;
  s = Trim(s);
  if (s.empty()) return false;
  std::string tmp(s);
  char* end = nullptr;
  errno = 0;
  const double v = std::strtod(tmp.c_str(), &end);
  if (errno != 0 || end == tmp.c_str()) return false;
  while (*end != '\0') {
    if (!std::isspace(static_cast<unsigned char>(*end))) return false;
    ++end;
  }
  *out = v;
  return true;
}

inline bool ParseU32(std::string_view s, u32* out) {
  if (!out) return false;
  s = Trim(s);
  if (s.empty()) return false;
  std::string tmp(s);
  char* end = nullptr;
  errno = 0;
  const unsigned long v = std::strtoul(tmp.c_str(), &end, 10);
  if (errno != 0 || end == tmp.c_str()) return false;
  while (*end != '\0') {
    if (!std::isspace(static_cast<unsigned char>(*end))) return false;
    ++end;
  }
  if (v > std::numeric_limits<u32>::max()) return false;
  *out = static_cast<u32>(v);
  return true;
}

template <int Dim, class T>
bool ReadBoxesSimple(const std::string& path,
                     Relation<Dim, T>* out_rel,
                     char sep = ',',
                     bool has_header = true,
                     std::string* err = nullptr) {
  if (!out_rel) { SetErr(err, "out_rel is null"); return false; }

  std::ifstream in(path);
  if (!in) {
    SetErr(err, "Cannot open CSV for reading: " + path);
    return false;
  }

  Relation<Dim, T> rel;
  std::string line;
  bool first_line = true;
  while (std::getline(in, line)) {
    std::string_view sv(line);
    sv = Trim(sv);
    if (sv.empty()) continue;
    if (!sv.empty() && sv.front() == '#') continue;

    if (first_line && has_header) {
      first_line = false;
      continue;
    }
    first_line = false;

    const auto cols = SplitLine(sv, sep);

    // Two supported shapes:
    //  A) with id: 1 + 2*Dim columns
    //  B) without id: 2*Dim columns
    if (cols.size() != static_cast<usize>(1 + 2 * Dim) &&
        cols.size() != static_cast<usize>(2 * Dim)) {
      SetErr(err, "Bad column count in line: expected " + std::to_string(2 * Dim) +
                  " or " + std::to_string(1 + 2 * Dim));
      return false;
    }

    usize offset = 0;
    std::optional<Id> id;
    if (cols.size() == static_cast<usize>(1 + 2 * Dim)) {
      u32 idv;
      if (!ParseU32(cols[0], &idv)) {
        SetErr(err, "Failed to parse id column");
        return false;
      }
      id = static_cast<Id>(idv);
      offset = 1;
    }

    Box<Dim, T> b;
    for (int i = 0; i < Dim; ++i) {
      double x;
      if (!ParseDouble(cols[offset + static_cast<usize>(i)], &x)) {
        SetErr(err, "Failed to parse lo coordinate");
        return false;
      }
      b.lo.v[static_cast<usize>(i)] = static_cast<T>(x);
    }
    for (int i = 0; i < Dim; ++i) {
      double x;
      if (!ParseDouble(cols[offset + static_cast<usize>(Dim + i)], &x)) {
        SetErr(err, "Failed to parse hi coordinate");
        return false;
      }
      b.hi.v[static_cast<usize>(i)] = static_cast<T>(x);
    }
    rel.Add(b, id);
  }

  *out_rel = std::move(rel);
  return true;
}

}  // namespace csv
}  // namespace sjs

#pragma once
// sjs/io/binary_io.h
//
// Custom binary format for fast, robust loading/saving of relations (sets of boxes).
//
// Why a custom format?
//  - CSV is large and slow to parse for big N.
//  - We want deterministic, stable parsing for experiments.
//  - We want to preserve stable ids and dimensionality.
//  - We want forward compatibility (versioned headers).
//
// Format (little-endian):
//   FileHeader (fixed size)
//   name_len: u32
//   name bytes (name_len bytes, not null-terminated)
//   Records (count times):
//     lo[dim] (scalar)
//     hi[dim] (scalar)
//     (optional) id (u32) if flags & kHasIds
//
// Notes:
//  - Scalars can be stored as float32 or float64; we currently write float64 by default.
//  - Half-open semantics are indicated by a header flag (kHalfOpen).
//  - Endianness is assumed little-endian; file header includes a marker to detect mismatch.
//
// This header provides:
//  - BinaryWriteOptions / BinaryReadOptions
//  - WriteRelationBinary / ReadRelationBinary
//  - Convenience: WriteDatasetBinaryPair / ReadDatasetBinaryPair

#include "sjs/core/types.h"
#include "sjs/core/assert.h"
#include "sjs/io/dataset.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sjs {
namespace binary {

inline constexpr u64 kEndianMarker = 0x0102030405060708ULL;
inline constexpr u32 kFormatVersion = 1;

// 8-byte magic: "SJSBOX\0\0"
inline constexpr char kMagic[8] = {'S','J','S','B','O','X','\0','\0'};

enum HeaderFlags : u32 {
  kHalfOpen = 1u << 0,
  kHasIds   = 1u << 1,
};

enum class ScalarEncoding : u32 {
  Float32 = 32,
  Float64 = 64,
};

#pragma pack(push, 1)
struct FileHeader {
  char magic[8];
  u32 version;
  u32 dim;
  u32 scalar_bits;   // 32 or 64
  u32 flags;         // bitset HeaderFlags
  u64 count;         // number of records
  u64 endian;        // kEndianMarker
  u64 reserved[4];   // for forward compatibility (set to 0)
};
#pragma pack(pop)

static_assert(sizeof(FileHeader) == 8 + 4*4 + 8*6, "Unexpected FileHeader size");

struct RelationFileInfo {
  std::string name;
  u32 version = 0;
  u32 dim = 0;
  ScalarEncoding scalar = ScalarEncoding::Float64;
  u32 flags = 0;
  u64 count = 0;
};

struct BinaryWriteOptions {
  bool write_ids = true;
  bool write_name = true;
  ScalarEncoding scalar = ScalarEncoding::Float64;
  bool half_open = true;
};

struct BinaryReadOptions {
  // If true, generate sequential ids when file doesn't contain ids.
  bool generate_ids_if_missing = true;

  // If true, drop empty boxes (lo >= hi) on load.
  bool drop_empty = false;
};

// -------------- helpers --------------
inline bool IsLittleEndianHost() noexcept {
  const u16 x = 1;
  return *reinterpret_cast<const u8*>(&x) == 1;
}

inline void SetErr(std::string* err, const std::string& msg) {
  if (err) *err = msg;
}

inline bool ReadExact(std::ifstream& in, void* dst, std::size_t n, std::string* err) {
  in.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n));
  if (!in) {
    SetErr(err, "Binary read failed (unexpected EOF or IO error)");
    return false;
  }
  return true;
}

inline bool WriteExact(std::ofstream& out, const void* src, std::size_t n, std::string* err) {
  out.write(reinterpret_cast<const char*>(src), static_cast<std::streamsize>(n));
  if (!out) {
    SetErr(err, "Binary write failed (IO error)");
    return false;
  }
  return true;
}

inline bool ReadU32(std::ifstream& in, u32* v, std::string* err) {
  return ReadExact(in, v, sizeof(u32), err);
}
inline bool WriteU32(std::ofstream& out, u32 v, std::string* err) {
  return WriteExact(out, &v, sizeof(u32), err);
}

template <class T>
inline bool ReadScalar(std::ifstream& in, T* v, std::string* err) {
  return ReadExact(in, v, sizeof(T), err);
}
template <class T>
inline bool WriteScalar(std::ofstream& out, const T& v, std::string* err) {
  return WriteExact(out, &v, sizeof(T), err);
}

// -------------- write --------------
template <int Dim, class T>
bool WriteRelationBinary(const std::string& path,
                         const Relation<Dim, T>& rel,
                         const BinaryWriteOptions& opt = {},
                         std::string* err = nullptr) {
  if (!IsLittleEndianHost()) {
    SetErr(err, "Host is not little-endian; binary format assumes little-endian");
    return false;
  }

  // Prepare ids if needed, without mutating input.
  std::vector<Id> tmp_ids;
  const bool has_ids = opt.write_ids;
  if (has_ids && rel.ids.empty()) {
    tmp_ids.resize(rel.boxes.size());
    for (usize i = 0; i < rel.boxes.size(); ++i) tmp_ids[i] = static_cast<Id>(i);
  }

  std::ofstream out(path, std::ios::binary);
  if (!out) {
    SetErr(err, "Cannot open file for writing: " + path);
    return false;
  }

  FileHeader h{};
  std::memcpy(h.magic, kMagic, sizeof(kMagic));
  h.version = kFormatVersion;
  h.dim = static_cast<u32>(Dim);
  h.scalar_bits = static_cast<u32>(opt.scalar == ScalarEncoding::Float32 ? 32 : 64);
  h.flags = 0;
  if (opt.half_open) h.flags |= kHalfOpen;
  if (has_ids) h.flags |= kHasIds;
  h.count = static_cast<u64>(rel.boxes.size());
  h.endian = kEndianMarker;
  for (auto& x : h.reserved) x = 0;

  if (!WriteExact(out, &h, sizeof(h), err)) return false;

  // name block
  if (opt.write_name) {
    const std::string& nm = rel.name;
    if (nm.size() > std::numeric_limits<u32>::max()) {
      SetErr(err, "Relation name too long to serialize");
      return false;
    }
    const u32 len = static_cast<u32>(nm.size());
    if (!WriteU32(out, len, err)) return false;
    if (len > 0) {
      if (!WriteExact(out, nm.data(), len, err)) return false;
    }
  } else {
    const u32 len = 0;
    if (!WriteU32(out, len, err)) return false;
  }

  // records
  if (opt.scalar == ScalarEncoding::Float64) {
    for (usize i = 0; i < rel.boxes.size(); ++i) {
      const auto& b = rel.boxes[i];
      for (int d = 0; d < Dim; ++d) {
        const double x = static_cast<double>(b.lo.v[d]);
        if (!WriteScalar(out, x, err)) return false;
      }
      for (int d = 0; d < Dim; ++d) {
        const double x = static_cast<double>(b.hi.v[d]);
        if (!WriteScalar(out, x, err)) return false;
      }
      if (has_ids) {
        const Id id = rel.ids.empty() ? tmp_ids[i] : rel.ids[i];
        if (!WriteScalar(out, id, err)) return false;
      }
    }
  } else {  // Float32
    for (usize i = 0; i < rel.boxes.size(); ++i) {
      const auto& b = rel.boxes[i];
      for (int d = 0; d < Dim; ++d) {
        const float x = static_cast<float>(b.lo.v[d]);
        if (!WriteScalar(out, x, err)) return false;
      }
      for (int d = 0; d < Dim; ++d) {
        const float x = static_cast<float>(b.hi.v[d]);
        if (!WriteScalar(out, x, err)) return false;
      }
      if (has_ids) {
        const Id id = rel.ids.empty() ? tmp_ids[i] : rel.ids[i];
        if (!WriteScalar(out, id, err)) return false;
      }
    }
  }

  return true;
}

// Convenience: write a Dataset<Dim> into two files (R and S).
template <int Dim, class T>
bool WriteDatasetBinaryPair(const std::string& path_r,
                            const std::string& path_s,
                            const Dataset<Dim, T>& ds,
                            const BinaryWriteOptions& opt = {},
                            std::string* err = nullptr) {
  if (!WriteRelationBinary<Dim, T>(path_r, ds.R, opt, err)) return false;
  if (!WriteRelationBinary<Dim, T>(path_s, ds.S, opt, err)) return false;
  return true;
}

// -------------- read --------------
template <int Dim, class T>
bool ReadRelationBinary(const std::string& path,
                        Relation<Dim, T>* out_rel,
                        RelationFileInfo* out_info = nullptr,
                        const BinaryReadOptions& opt = {},
                        std::string* err = nullptr) {
  if (!out_rel) {
    SetErr(err, "out_rel is null");
    return false;
  }
  if (!IsLittleEndianHost()) {
    SetErr(err, "Host is not little-endian; binary format assumes little-endian");
    return false;
  }

  std::ifstream in(path, std::ios::binary);
  if (!in) {
    SetErr(err, "Cannot open file for reading: " + path);
    return false;
  }

  FileHeader h{};
  if (!ReadExact(in, &h, sizeof(h), err)) return false;

  if (std::memcmp(h.magic, kMagic, sizeof(kMagic)) != 0) {
    SetErr(err, "Bad magic header (not an SJSBOX file)");
    return false;
  }
  if (h.endian != kEndianMarker) {
    SetErr(err, "Endianness marker mismatch; file may be from different endian host");
    return false;
  }
  if (h.version == 0 || h.version > kFormatVersion) {
    SetErr(err, "Unsupported binary format version: " + std::to_string(h.version));
    return false;
  }
  if (h.dim != static_cast<u32>(Dim)) {
    SetErr(err, "Dimension mismatch: file dim=" + std::to_string(h.dim) +
                    " but compiled Dim=" + std::to_string(Dim));
    return false;
  }
  if (!(h.scalar_bits == 32 || h.scalar_bits == 64)) {
    SetErr(err, "Unsupported scalar_bits in file: " + std::to_string(h.scalar_bits));
    return false;
  }
  const bool file_half_open = (h.flags & kHalfOpen) != 0;
  const bool file_has_ids = (h.flags & kHasIds) != 0;

  u32 name_len = 0;
  if (!ReadU32(in, &name_len, err)) return false;

  std::string name;
  name.resize(name_len);
  if (name_len > 0) {
    if (!ReadExact(in, name.data(), name_len, err)) return false;
  }

  if (out_info) {
    out_info->name = name;
    out_info->version = h.version;
    out_info->dim = h.dim;
    out_info->scalar = (h.scalar_bits == 32) ? ScalarEncoding::Float32 : ScalarEncoding::Float64;
    out_info->flags = h.flags;
    out_info->count = h.count;
  }

  Relation<Dim, T> rel;
  rel.name = name;
  rel.boxes.resize(static_cast<usize>(h.count));
  if (file_has_ids) rel.ids.resize(static_cast<usize>(h.count));

  // Read records with conversion if needed.
  if (h.scalar_bits == 64) {
    for (u64 i = 0; i < h.count; ++i) {
      Box<Dim, T> b;
      for (int d = 0; d < Dim; ++d) {
        double x;
        if (!ReadScalar(in, &x, err)) return false;
        b.lo.v[static_cast<usize>(d)] = static_cast<T>(x);
      }
      for (int d = 0; d < Dim; ++d) {
        double x;
        if (!ReadScalar(in, &x, err)) return false;
        b.hi.v[static_cast<usize>(d)] = static_cast<T>(x);
      }
      rel.boxes[static_cast<usize>(i)] = b;
      if (file_has_ids) {
        Id id;
        if (!ReadScalar(in, &id, err)) return false;
        rel.ids[static_cast<usize>(i)] = id;
      }
    }
  } else {  // 32
    for (u64 i = 0; i < h.count; ++i) {
      Box<Dim, T> b;
      for (int d = 0; d < Dim; ++d) {
        float x;
        if (!ReadScalar(in, &x, err)) return false;
        b.lo.v[static_cast<usize>(d)] = static_cast<T>(x);
      }
      for (int d = 0; d < Dim; ++d) {
        float x;
        if (!ReadScalar(in, &x, err)) return false;
        b.hi.v[static_cast<usize>(d)] = static_cast<T>(x);
      }
      rel.boxes[static_cast<usize>(i)] = b;
      if (file_has_ids) {
        Id id;
        if (!ReadScalar(in, &id, err)) return false;
        rel.ids[static_cast<usize>(i)] = id;
      }
    }
  }

  // Post-processing options
  if (!file_has_ids && opt.generate_ids_if_missing) {
    rel.EnsureIds();
  }
  if (opt.drop_empty) {
    rel.RemoveEmptyBoxes();
  }

  *out_rel = std::move(rel);

  // If file wasn't half-open, surface a warning via err string (non-fatal).
  if (!file_half_open && err) {
    *err = "Warning: file is not marked half-open; project assumes half-open semantics";
  }

  return true;
}

template <int Dim, class T>
bool ReadDatasetBinaryPair(const std::string& path_r,
                           const std::string& path_s,
                           Dataset<Dim, T>* out_ds,
                           const BinaryReadOptions& opt = {},
                           std::string* err = nullptr) {
  if (!out_ds) {
    SetErr(err, "out_ds is null");
    return false;
  }
  Dataset<Dim, T> ds;
  ds.half_open = true;

  if (!ReadRelationBinary<Dim, T>(path_r, &ds.R, nullptr, opt, err)) return false;
  if (!ReadRelationBinary<Dim, T>(path_s, &ds.S, nullptr, opt, err)) return false;
  *out_ds = std::move(ds);
  return true;
}

}  // namespace binary
}  // namespace sjs

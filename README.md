# Spatial Join Sampling (C++)

Lightweight C++ framework for running spatial join sampling experiments. The
project builds a static library (`sjs`) plus a small collection of command-line
tools (e.g., `sjs_run`, `sjs_sweep`).

## Quick start

```bash
mkdir -p build
cd build
cmake ..                # defaults to Release if not specified
cmake --build . -j

# List available flags
./sjs_run --help
```

## Notes

- CMake 3.16+ and a C++17-capable compiler are required.
- Tests live under `tests/` and can be enabled with `-DSJS_BUILD_TESTS=ON`.
- Example datasets and configs are stored under `data/` and `config/`.

#pragma once
// sjs/core/logging.h
//
// Lightweight logging for experiments (no fmt dependency).
// Usage:
//   SJS_LOG_INFO("hello", x, y);
//   sjs::Logger::Instance().SetLevel(sjs::LogLevel::Debug);
//
// By default logs go to stderr.

#include <atomic>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

namespace sjs {

enum class LogLevel : std::uint8_t {
  Trace = 0,
  Debug = 1,
  Info  = 2,
  Warn  = 3,
  Error = 4,
  Off   = 5,
};

inline constexpr std::string_view ToString(LogLevel lvl) noexcept {
  switch (lvl) {
    case LogLevel::Trace: return "TRACE";
    case LogLevel::Debug: return "DEBUG";
    case LogLevel::Info:  return "INFO";
    case LogLevel::Warn:  return "WARN";
    case LogLevel::Error: return "ERROR";
    case LogLevel::Off:   return "OFF";
  }
  return "UNKNOWN";
}

struct LoggingConfig {
  LogLevel level = LogLevel::Info;
  bool with_timestamp = true;
  bool with_thread_id = false;
};

namespace detail {

inline std::tm LocalTime(std::time_t t) {
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  return tm;
}

template <class... Args>
inline void AppendSpaceSeparated(std::ostringstream& oss, Args&&... args) {
  bool first = true;
  auto append_one = [&](auto&& v) {
    if (!first) oss << ' ';
    first = false;
    oss << std::forward<decltype(v)>(v);
  };
  (append_one(std::forward<Args>(args)), ...);
}

}  // namespace detail

class Logger {
 public:
  static Logger& Instance() {
    static Logger inst;
    return inst;
  }

  void SetConfig(const LoggingConfig& cfg) {
    std::lock_guard<std::mutex> lk(mu_);
    cfg_ = cfg;
    level_.store(cfg.level, std::memory_order_relaxed);
  }

  void SetLevel(LogLevel lvl) { level_.store(lvl, std::memory_order_relaxed); }

  LogLevel Level() const { return level_.load(std::memory_order_relaxed); }

  // Set output stream (default: &std::cerr). Caller owns the stream.
  void SetOutput(std::ostream* out) {
    std::lock_guard<std::mutex> lk(mu_);
    out_ = out ? out : &std::cerr;
  }

  template <class... Args>
  void Log(LogLevel lvl, Args&&... args) {
    const LogLevel cur = level_.load(std::memory_order_relaxed);
    if (lvl < cur || cur == LogLevel::Off) return;

    std::ostringstream oss;
    {
      // prefix
      std::lock_guard<std::mutex> lk(mu_);  // serialize and also read cfg_/out_ safely

      if (cfg_.with_timestamp) {
        using clock = std::chrono::system_clock;
        const auto now = clock::now();
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

        const std::time_t t = clock::to_time_t(now);
        const std::tm tm = detail::LocalTime(t);

        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S")
            << '.' << std::setfill('0') << std::setw(3) << ms.count()
            << ' ';
      }

      oss << '[' << ToString(lvl) << ']';

      if (cfg_.with_thread_id) {
        oss << "[tid=" << std::this_thread::get_id() << ']';
      }

      oss << ' ';
      detail::AppendSpaceSeparated(oss, std::forward<Args>(args)...);
      oss << '\n';

      (*out_) << oss.str() << std::flush;
    }
  }

 private:
  Logger() : level_(LogLevel::Info), cfg_{}, out_(&std::cerr) {}

  std::atomic<LogLevel> level_;
  LoggingConfig cfg_;
  std::ostream* out_;
  mutable std::mutex mu_;
};

// Convenience macros
#define SJS_LOG_TRACE(...) ::sjs::Logger::Instance().Log(::sjs::LogLevel::Trace, __VA_ARGS__)
#define SJS_LOG_DEBUG(...) ::sjs::Logger::Instance().Log(::sjs::LogLevel::Debug, __VA_ARGS__)
#define SJS_LOG_INFO(...)  ::sjs::Logger::Instance().Log(::sjs::LogLevel::Info,  __VA_ARGS__)
#define SJS_LOG_WARN(...)  ::sjs::Logger::Instance().Log(::sjs::LogLevel::Warn,  __VA_ARGS__)
#define SJS_LOG_ERROR(...) ::sjs::Logger::Instance().Log(::sjs::LogLevel::Error, __VA_ARGS__)

}  // namespace sjs

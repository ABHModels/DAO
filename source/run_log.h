#ifndef DAO_RUN_LOG_H
#define DAO_RUN_LOG_H

#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

// Diagnostic output only: this state never participates in the solver or hash.
// Header-only so standalone numerical tools do not acquire a link dependency.
namespace dao_log {
inline std::ofstream file;
inline std::mutex mutex;
inline bool verbose = false;
inline pid_t owner = 0;
inline bool write_warning = false;

inline void open(const std::string& path, bool echo_details)
{
    std::lock_guard<std::mutex> lock(mutex);
    if(file.is_open()) file.close();
    file.clear();
    file.open(path, std::ios::app);
    if(!file) throw std::runtime_error("Cannot open diagnostic log: " + path);
    verbose = echo_details;
    owner = getpid();
    write_warning = false;
}

inline void write(bool detail, bool error, const char* format, va_list args)
{
    va_list copy;
    va_copy(copy, args);
    const int length = std::vsnprintf(nullptr, 0, format, copy);
    va_end(copy);
    if(length < 0) return;
    std::vector<char> text(static_cast<size_t>(length) + 1);
    std::vsnprintf(text.data(), text.size(), format, args);

    // Cloudy children have their own redirected stdout/stderr. Do not let a
    // forked worker write through the parent's inherited log stream or mutex.
    if(owner != 0 && owner != getpid()) {
        FILE* output = error ? stderr : stdout;
        std::fwrite(text.data(), 1, length, output);
        std::fflush(output);
        return;
    }
    std::lock_guard<std::mutex> lock(mutex);
    if(file.is_open()) {
        file.write(text.data(), length);
        file.flush();
        if(!file && !write_warning) {
            std::fprintf(stderr, "Warning: diagnostic log write failed; check disk space.\n");
            write_warning = true;
        }
    }
    // Standalone tools which do not configure a run log retain their output.
    if(!detail || verbose || owner == 0) {
        FILE* output = error ? stderr : stdout;
        std::fwrite(text.data(), 1, length, output);
        std::fflush(output);
    }
}

inline void detail(const char* format, ...)
{
    va_list args; va_start(args, format);
    write(true, false, format, args);
    va_end(args);
}
inline void info(const char* format, ...)
{
    va_list args; va_start(args, format);
    write(false, false, format, args);
    va_end(args);
}
inline void error(const char* format, ...)
{
    va_list args; va_start(args, format);
    write(false, true, format, args);
    va_end(args);
}
} // namespace dao_log
#endif

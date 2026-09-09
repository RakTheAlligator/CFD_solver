#pragma once

#include <filesystem>
#include <fstream>

namespace cfd
{

struct SimpleIterationInfo;

/// Writes completed SIMPLE-iteration diagnostics to a live-readable CSV file.
///
/// The destination is truncated and initialized with a fixed header at
/// construction. Each completed data row is flushed synchronously so a
/// separate process can observe convergence while the solver runs.
class SimpleConvergenceCsvWriter
{
  public:
    /// Opens `file_path` and writes the fixed convergence header.
    ///
    /// The parent directory must already exist.
    ///
    /// @throws std::runtime_error If the file cannot be opened, written, or
    ///         flushed.
    explicit SimpleConvergenceCsvWriter(const std::filesystem::path &file_path);

    SimpleConvergenceCsvWriter(const SimpleConvergenceCsvWriter &) = delete;
    SimpleConvergenceCsvWriter &operator=(const SimpleConvergenceCsvWriter &) = delete;
    SimpleConvergenceCsvWriter(SimpleConvergenceCsvWriter &&) = delete;
    SimpleConvergenceCsvWriter &operator=(SimpleConvergenceCsvWriter &&) = delete;

    ~SimpleConvergenceCsvWriter() = default;

    /// Appends and flushes one completed SIMPLE iteration.
    ///
    /// @throws std::runtime_error If the row cannot be written or flushed.
    void write(const SimpleIterationInfo &info);

  private:
    std::ofstream output_;
};

} // namespace cfd

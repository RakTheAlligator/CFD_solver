#pragma once

#include <filesystem>
#include <fstream>

namespace cfd
{
struct EulerSteadyObservation;

/// Checked, flushed Euler observation CSV, separate from SIMPLE's schema.
/// Parent must exist. Constructor truncates; valid writes allocate no field buffers.
class EulerSteadyCsvWriter
{
  public:
    /// @throws std::runtime_error On open/header/write/flush failure.
    explicit EulerSteadyCsvWriter(const std::filesystem::path &path);
    EulerSteadyCsvWriter(const EulerSteadyCsvWriter &) = delete;
    EulerSteadyCsvWriter &operator=(const EulerSteadyCsvWriter &) = delete;
    EulerSteadyCsvWriter(EulerSteadyCsvWriter &&) = delete;
    EulerSteadyCsvWriter &operator=(EulerSteadyCsvWriter &&) = delete;
    ~EulerSteadyCsvWriter() = default;
    /// @throws std::runtime_error On write/flush failure.
    void write(const EulerSteadyObservation &observation);

  private:
    std::ofstream output_;
};
} // namespace cfd

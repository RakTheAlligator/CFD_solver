#include "app/LiveConvergenceLauncher.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace cfd::app
{
namespace
{

void warn(const std::string_view message)
{
    std::cerr << "[Warning]\n  Live convergence plotter: " << message << '\n';
}

[[noreturn]]
void report_child_error_and_exit(const int error_pipe, const int error_number) noexcept
{
    ssize_t write_result{::write(error_pipe, &error_number, sizeof(error_number))};
    while (write_result == -1 && errno == EINTR)
    {
        write_result = ::write(error_pipe, &error_number, sizeof(error_number));
    }
    ::_exit(127);
}

void wait_for_intermediate_child(const pid_t child_process)
{
    int child_status{};
    while (::waitpid(child_process, &child_status, 0) == -1)
    {
        if (errno != EINTR)
        {
            throw std::runtime_error(std::string{"unable to reap launcher child: "} + std::strerror(errno));
        }
    }
}

void launch(const std::filesystem::path &convergence_file)
{
    std::error_code filesystem_error;
    const std::filesystem::path executable_path{std::filesystem::canonical("/proc/self/exe", filesystem_error)};
    if (filesystem_error)
    {
        throw std::runtime_error("unable to locate CFD_solver executable: " + filesystem_error.message());
    }

    const std::filesystem::path repository_root{executable_path.parent_path().parent_path()};
    const std::filesystem::path plotter_path{repository_root / "tools" / "live_convergence.py"};
    const bool plotter_exists{std::filesystem::is_regular_file(plotter_path, filesystem_error)};
    if (filesystem_error)
    {
        throw std::runtime_error("unable to inspect " + plotter_path.string() + ": " + filesystem_error.message());
    }
    if (!plotter_exists)
    {
        throw std::runtime_error("script not found at " + plotter_path.string());
    }

    const std::filesystem::path absolute_convergence_file{
        std::filesystem::absolute(convergence_file, filesystem_error)};
    if (filesystem_error)
    {
        throw std::runtime_error("unable to resolve convergence CSV path: " + filesystem_error.message());
    }
    std::string plotter_argument{plotter_path.string()};
    std::string convergence_argument{absolute_convergence_file.string()};

    std::array<int, 2> error_pipe{};
    if (::pipe2(error_pipe.data(), O_CLOEXEC) == -1)
    {
        throw std::runtime_error(std::string{"unable to create launch-status pipe: "} + std::strerror(errno));
    }

    const pid_t intermediate_child{::fork()};
    if (intermediate_child == -1)
    {
        const int error_number{errno};
        static_cast<void>(::close(error_pipe[0]));
        static_cast<void>(::close(error_pipe[1]));
        throw std::runtime_error(std::string{"unable to fork launcher child: "} + std::strerror(error_number));
    }

    if (intermediate_child == 0)
    {
        static_cast<void>(::close(error_pipe[0]));
        if (::setsid() == -1)
        {
            report_child_error_and_exit(error_pipe[1], errno);
        }

        // The parent reaps this intermediate child. The detached grandchild is
        // adopted by init, so CFD_solver neither waits for it nor retains a zombie.
        const pid_t detached_child{::fork()};
        if (detached_child == -1)
        {
            report_child_error_and_exit(error_pipe[1], errno);
        }
        if (detached_child != 0)
        {
            ::_exit(0);
        }

        std::string interpreter{"python3"};
        std::array<char *, 4> arguments{
            interpreter.data(),
            plotter_argument.data(),
            convergence_argument.data(),
            nullptr,
        };
        ::execvp(arguments[0], arguments.data());
        report_child_error_and_exit(error_pipe[1], errno);
    }

    static_cast<void>(::close(error_pipe[1]));
    wait_for_intermediate_child(intermediate_child);

    int launch_error{};
    ssize_t read_result{::read(error_pipe[0], &launch_error, sizeof(launch_error))};
    while (read_result == -1 && errno == EINTR)
    {
        read_result = ::read(error_pipe[0], &launch_error, sizeof(launch_error));
    }
    static_cast<void>(::close(error_pipe[0]));

    if (read_result == -1)
    {
        throw std::runtime_error(std::string{"unable to read launch status: "} + std::strerror(errno));
    }
    if (read_result != 0 && read_result != static_cast<ssize_t>(sizeof(launch_error)))
    {
        throw std::runtime_error("received an incomplete detached-process launch status");
    }
    if (read_result != 0)
    {
        throw std::runtime_error(std::string{"unable to launch detached python3 process: "} +
                                 std::strerror(launch_error));
    }
}

} // namespace

void launch_live_convergence_plotter(const std::filesystem::path &convergence_file) noexcept
{
    try
    {
        launch(convergence_file);
    }
    catch (const std::exception &error)
    {
        warn(error.what());
    }
}

} // namespace cfd::app

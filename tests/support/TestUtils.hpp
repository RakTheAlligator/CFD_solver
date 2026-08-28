#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace cfd::test
{

inline constexpr double test_tolerance{1.0e-12};

[[noreturn]]
inline void fail(const std::string &message)
{
    throw std::runtime_error(message);
}

inline void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        fail(message);
    }
}

inline void require_near(const double actual, const double expected, const double tolerance, const std::string &message)
{
    if (!std::isfinite(actual) || !std::isfinite(expected))
    {
        fail(message + " (comparison requires finite values).");
    }

    if (!std::isfinite(tolerance) || tolerance < 0.0)
    {
        fail(message + " (tolerance must be finite and non-negative).");
    }

    // Combine absolute and relative behavior: values below unit scale use the
    // supplied tolerance directly, while larger values scale it with magnitude.
    const double scale{std::max({1.0, std::abs(actual), std::abs(expected)})};

    if (std::abs(actual - expected) > tolerance * scale)
    {
        fail(message + " (actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected) + ")");
    }
}

inline void require_contains(const std::string &text, const std::string_view expected, const std::string &message)
{
    if (text.find(expected) == std::string::npos)
    {
        fail(message);
    }
}

[[nodiscard]]
inline std::string read_text_file(const std::filesystem::path &file_path)
{
    std::ifstream input{file_path};

    if (!input)
    {
        fail("Unable to open test file: " + file_path.string());
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    if (!input.eof() && input.fail())
    {
        fail("Error while reading test file: " + file_path.string());
    }

    return buffer.str();
}

template <typename Exception = std::exception, typename Function>
void require_throws(Function &&function, const std::string &message)
{
    try
    {
        std::forward<Function>(function)();
    }
    catch (const Exception &)
    {
        return;
    }
    catch (...)
    {
        fail(message + " (unexpected exception type).");
    }

    fail(message + " (no exception was thrown).");
}

template <typename Exception = std::exception, typename Function>
void require_throws_with_message(Function &&function, const std::string_view expected_message,
                                 const std::string &failure_message)
{
    try
    {
        std::forward<Function>(function)();
    }
    catch (const Exception &error)
    {
        const std::string_view actual_message{error.what()};

        if (actual_message.find(expected_message) == std::string_view::npos)
        {
            fail(failure_message + " (unexpected exception message: \"" + error.what() + "\").");
        }

        return;
    }
    catch (...)
    {
        fail(failure_message + " (unexpected exception type).");
    }

    fail(failure_message + " (no exception was thrown).");
}

template <typename TestFunction> int run_test(const std::string_view name, TestFunction test_function)
{
    try
    {
        test_function();

        std::cout << "[PASS] " << name << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] " << name << ": " << error.what() << '\n';
        return 1;
    }
}

inline int finish_tests(const int failure_count, const std::string_view suite_name)
{
    if (failure_count == 0)
    {
        std::cout << "[PASS] All " << suite_name << " tests passed.\n";
        return 0;
    }

    std::cerr << "[FAIL] " << failure_count << ' ' << suite_name << " test(s) failed.\n";

    return 1;
}

} // namespace cfd::test
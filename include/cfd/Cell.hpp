#pragma once

#include <cstdint>

namespace cfd {

enum class CellType : std::uint8_t {
    Triangle,
    Quadrilateral
};

} // namespace cfd

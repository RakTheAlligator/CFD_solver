#pragma once

#include "cfd/Types.hpp"

#include <array>

namespace cfd {

struct Face {
    std::array<Index, 2> node_ids{};
};

} // namespace cfd

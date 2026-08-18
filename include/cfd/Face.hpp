#pragma once

#include "cfd/Types.hpp"

#include <array>

namespace cfd {

struct Face {
    std::array<Index, 2> node_ids{};
};
struct FaceAdjacency {
    Index owner{invalid_index};
    Index neighbor{invalid_index};

    [[nodiscard]]
    bool is_boundary() const noexcept
    {
        return neighbor == invalid_index;
    }
};

} // namespace cfd
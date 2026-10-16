#include "cfd/compressible/ConservativeField.hpp"

namespace cfd
{

ConservativeField::ConservativeField(const Index cell_count, const ConservativeState &initial_value)
    : density_(cell_count, initial_value.density), momentum_x_(cell_count, initial_value.momentum_x),
      momentum_r_(cell_count, initial_value.momentum_r),
      total_energy_density_(cell_count, initial_value.total_energy_density)
{
}

ConservativeState ConservativeField::state(const Index cell) const noexcept
{
    return {density_[cell], momentum_x_[cell], momentum_r_[cell], total_energy_density_[cell]};
}

void ConservativeField::set_state(const Index cell, const ConservativeState &value) noexcept
{
    density_[cell] = value.density;
    momentum_x_[cell] = value.momentum_x;
    momentum_r_[cell] = value.momentum_r;
    total_energy_density_[cell] = value.total_energy_density;
}

} // namespace cfd

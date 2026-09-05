#pragma once

#include <vector>

namespace cfd
{

class CellMomentumPressureResponse;
class CellScalarField;
class CellVectorField;
class CellVelocityField;
class FaceFluxField;
class FacePressureResponseField;
class Mesh;

/// Momentum-weighted mass-flux interpolation for internal mesh faces.
///
/// The resulting flux is integrated and follows the Mesh owner orientation.
/// This component handles internal faces only; boundary entries in both output
/// fields are left exactly unchanged for a future boundary-flux component.
///
/// @note The referenced Mesh is not owned and must outlive this object.
/// @note Construction allocates one fixed geometry cache. Repeated valid
///       updates perform no dynamic allocation and copy no field-sized arrays.
class RhieChowInternalFaceInterpolation
{
  public:
    /// Constructs interpolation data for a fixed Mesh and constant density.
    ///
    /// @throws std::invalid_argument If `density` is not finite and strictly
    ///         positive.
    /// @throws std::runtime_error If an internal face has unusable geometry.
    RhieChowInternalFaceInterpolation(const Mesh &mesh, double density);

    RhieChowInternalFaceInterpolation(const RhieChowInternalFaceInterpolation &) = delete;
    RhieChowInternalFaceInterpolation &operator=(const RhieChowInternalFaceInterpolation &) = delete;

    RhieChowInternalFaceInterpolation(RhieChowInternalFaceInterpolation &&) noexcept = default;
    RhieChowInternalFaceInterpolation &operator=(RhieChowInternalFaceInterpolation &&) noexcept = delete;

    ~RhieChowInternalFaceInterpolation() = default;

    /// Overwrites the mass flux and pressure response on every internal face.
    ///
    /// Boundary entries are left exactly unchanged. All inputs and computed
    /// internal-face results are validated before either output is modified.
    ///
    /// @throws std::invalid_argument If any field cardinality is incompatible
    ///         with the Mesh.
    /// @throws std::runtime_error If a used input or computed face quantity is
    ///         not finite, or a required pressure response is not positive.
    void update_internal_faces(const CellVelocityField &velocity, const CellScalarField &pressure,
                               const CellVectorField &pressure_gradient,
                               const CellMomentumPressureResponse &momentum_response, FaceFluxField &mass_flux,
                               FacePressureResponseField &face_pressure_response) const;

  private:
    struct InternalFaceGeometry
    {
        double interpolation_weight{};
        double inverse_area_dot_displacement{};
    };

    const Mesh *mesh_;
    double density_;
    std::vector<InternalFaceGeometry> internal_face_geometry_;
};

} // namespace cfd

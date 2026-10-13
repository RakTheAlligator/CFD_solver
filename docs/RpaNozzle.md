# RPA nozzle mesh input

This is a **meshing-only** Cartesian meridian `(x,r)` workflow. It does not
implement axisymmetric metrics, compressible flow or Euler. Do not run the
current incompressible SIMPLE application as a nozzle-flow model.

## Dependency and invocation

The CAD adapter uses Python 3 and `ezdxf` 1.4.3 for DXF entities and block
transformations. Gmsh C++ remains the only mesh generator; the existing
RawMeshData, MeshBuilder, statistics and VTU writer are reused.

```sh
python3 -m venv /tmp/cfd-cad
/tmp/cfd-cad/bin/python -m pip install -r tools/requirements-cad.txt
cmake -S . -B build -DPython3_EXECUTABLE=/tmp/cfd-cad/bin/python
cmake --build build --target cfd_mesh_preview
./build/cfd_mesh_preview cases/rpaNozzle
```

There is no manual Python conversion step. C++ invokes the configured
interpreter/script with an argument vector, without a shell, and propagates
failed imports. This implementation uses the project's existing POSIX/Linux
process conventions. The build must retain access to the source-tree adapter
and configured interpreter. RPA tests are enabled when that interpreter has
ezdxf; CMake warns explicitly when they are unavailable. Existing geometry
tests remain enabled regardless.

Output is `cases/rpaNozzle/results/mesh.vtu`. No CFD field files are supplied.

## meshDict

```text
geometry
{
    type rpaNozzle;
    file "constant/geometry/engine.dxf";
}
mesh
{
    cellType quadrilateral; // or triangle
    size 0.01;
    automaticMeshing true;
}
```

The standard FoamFile dictionary header is also required. The DXF path is
relative to the case directory, never the launch directory. Coordinates are
converted from `$INSUNITS` into metres. Missing/ambiguous units are errors;
`dxfUnits mm;` (also `cm`, `m`, `in`) can supply absent units but cannot conflict
with declared units. Drawing dimensions/text are not geometric inputs.

LINE, ARC, arbitrary-degree SPLINE, POLYLINE/LWPOLYLINE (including bulges) and
nested INSERT transforms are supported. Elliptical arcs resulting from
non-uniform block scaling are also flattened. An explicit centreline (CENTER
linetype or axis/centerline layer) defines the local axis; otherwise drawing
x/y are taken as x/r. The positive meridian must be a unique continuous open
chain with strictly increasing x and positive radius, an interior minimum
radius and a downstream exit. Optional mirrored negative meridians are ignored.
If only the negative side is present (e.g. an axis-normal reversal in a reflected
block), its normal is flipped into the positive radial convention.
Unsupported geometric entities are rejected rather than silently discarded;
annotations are ignored before block expansion. Non-planar, branched,
disconnected, backtracking or closed/self-intersecting profiles are rejected.

`geometricTolerance` defaults to `1e-6` metres for CAD curve flattening.
Source-entity junctions remain mandatory Gmsh vertices. Interior curve samples
describe polylines, not an assumed quadratic fit. This CAD approximation is
distinct from the additional chord error of a first-order FV boundary mesh.
Conical and non-parabolic cubic contours are tested. **No real RPA TIC export
has been validated**; support is limited to the stated entity/profile contract.

## Domain conventions

Absent geometry values select:

| Entry | Automatic value |
|---|---|
| `inletPosition` | first usable contour x |
| `outerUpstreamLength` | `min(max(Le, 1.5*De), exit_x-inlet_x)` |
| `wallThickness` | `0.01*De`, a **radial** offset |
| `downstreamLength` | `10*De` |
| `farfieldRadius` | `5*De`, measured from the axis |

`De=2*r_exit`, `Le=exit_x-throat_x`. All overrides are metres.
The inlet must remain upstream of the throat. Exterior length must fit the
retained contour; thickness/lengths/radius must be finite and positive, and
farfield radius must exceed the whole exterior wall. Unrepresentable extents
are rejected. The offset wall is a fictitious CFD obstruction, **not a mechanical
nozzle design**. It contains no solid cells. Interior and exterior fluid connect
around the finite exit lip; the nozzle exit plane is not a CFD boundary.

Physical boundaries: `inlet`, `axis`, `inner_wall`, `outer_wall`, `nozzle_lip`,
`farfield_upstream`, `farfield_radial`, `farfield_downstream`; surface: `fluid`.
Gmsh subdivision completes recombination into an all-QUAD mesh when requested.
Unexpected final element types are still errors; no mixed mesh is accepted.

## Automatic size policy

Common options are reused: `maximumGrowthRate`, `wallRefinementFactor`,
`wallRefinementLayers`. Region factors relative to `size` are:

| Mesh entry | Default |
|---|---:|
| `convergentRefinementFactor` | 0.75 |
| `throatRefinementFactor` | 0.25 |
| `divergentRefinementFactor` | 0.75 |
| `exitRefinementFactor` | 0.25 |
| `plumeRefinementFactor` | 0.5 |

All factors must be in `(0,1]`. Optional positive `throatRefinementLength`,
`exitRefinementLength`, `nearPlumeLength` specify full axial zone extents;
defaults use detected throat/exit dimensions and base size, clipped to the
available domain. In particular near plume defaults to `2*De`.

Distance/threshold and smoothly blended box fields are combined by minimum.
The short exterior-start edge and lip have geometry-driven, smoothly extended
sizes. Far plume returns to the base size. `automaticMeshing false` disables
all regional fields, not mandatory CAD junctions or all-QUAD completion.

`maximumGrowthRate` controls requested size-field ramp slopes. It is **not**
a guaranteed cell-to-cell ratio after unstructured triangulation/recombination.
The test checks actual wall/throat/exit sizes and reports actual neighboring-cell
ratios; the preview reports mesh quality and worst ratios. Thin lips and sharp
CAD junctions may require finer cells than the nominal target.

Adapter protocol `CFD_RPA_PROFILE_V1` is an ordered point count/list of metre
coordinates followed by `CURVE_JUNCTIONS`, count and ordered point indices.
The C++ reader validates cardinalities, finiteness, positive radii and strictly
ordered junctions before building the domain. Nonzero adapter status carries
its diagnostic rather than falling back to a different geometry.

References: [ezdxf units](https://ezdxf.readthedocs.io/en/stable/concepts/units.html),
[ezdxf decomposition](https://ezdxf.readthedocs.io/en/stable/disassemble.html),
[Gmsh quadrilateral meshing](https://gmsh.info/doc/texinfo/gmsh.html).

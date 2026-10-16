# RPA meridian mesh

The DXF is the source of truth for the engine profile. The Python 3 adapter
uses the existing `ezdxf==1.4.3` dependency (`tools/requirements-cad.txt`).
Configure its interpreter with `-DPython3_EXECUTABLE=/path/to/python`.
The C++ importer keeps the flattened curve, tolerance and CAD junctions;
it validates continuity, increasing x, positive radii and an interior throat.
Paths are case-relative; the adapter is invoked without a shell.

This is a Cartesian (x,r) fluid mesh, not an axisymmetric CFD formulation.
The fictitious shell is excluded from the fluid domain.

## Configuration

The demo configuration below matches `cases/rpaNozzle/system/meshDict`.

```foam
FoamFile { format ascii; class dictionary; object meshDict; }
geometry
{
    type rpaNozzle;
    engine
    {
        file "constant/geometry/engine.dxf";
        wallThicknessFactor 0.02;
        // Optional: dxfUnits auto; geometricTolerance 1e-6;
    }
    box
    {
        upstreamFactor auto;
        downstreamFactor 20.0;
        outerRadiusFactor 10.0;
    }
}
regions
{
    bandThicknessFactor 0.5;
    plumeOutletRadiusFactor 6.0;
}
mesh
{
    cellType hybrid;
    baseCellSize 0.01;
    minCellQuality 0.05;
}
```

All factors multiply the **inner exit radius**, not the diameter.
Internal box defaults are `upstreamFactor auto`, `downstreamFactor 80`,
`outerRadiusFactor 65`. The demo deliberately uses the smaller 20/10 box.
Auto upstream uses the complete available retained engine length. Explicit
upstream extents cannot exceed that profile. Shell thickness defaults to
0.02 exit radii. No input factor is silently clipped.

`regions` is optional; defaults are band thickness 0.5 exit radii and plume
outlet radius **16.0** exit radii. The demo explicitly retains its factor **6.0**
and requests HYBRID, 0.01 m and angular threshold 0.05. The inlet plume radius
is exit radius + shell thickness + band thickness. Its upper boundary is the
straight line to the outlet. Require inlet plume radius < outlet plume radius
< box radius, and no BAND/FARFIELD intersection.

Old flat RPA geometry entries and regional refinement settings are errors.
Use `geometry.engine`, `geometry.box` and `regions` instead. In particular
`mesh.size` becomes `baseCellSize`; throat/convergent/divergent/exit factors,
near-plume/corridor lengths, opening slopes, resolution policies/levels and
directional overrides are removed. Rectangle/BFS syntax is unchanged.

## Responsibilities

- `RpaNozzleGeometry` / `RpaNozzleDomain`: imported physical profile, shell,
  box, tolerance and mandatory CAD junctions; no cell counts.
- `RpaRegionLayout`: the five geometric regions and their interfaces;
  no resolution choices.
- `RpaMeshDiscretization`: radial distributions, column search and generation-only
  grids; no Gmsh calls. `GmshMesher` installs these grids with shared boundary
  nodes and uses Gmsh only for the remaining transition. No connectivity walk
  or post-generation overwrite is needed to recover prescribed grid positions.

| rpa_region | Region | hybrid | quadrilateral | triangle |
|---|---|---|---|---|
| 1 | ENGINE | QUAD | QUAD | split structured TRI |
| 2 | BAND | QUAD | QUAD | split structured TRI |
| 3 | PLUME | QUAD | QUAD | split structured TRI |
| 4 | TRANSITION | TRI | recombined QUAD | TRI |
| 5 | FARFIELD | QUAD | QUAD | split structured TRI |

There is one logical plume, not multiple corridors. CAD curve subdivisions
preserve source junctions and do not reset distributions. The TRI mode uses the
same geometric/discretization pipeline and splits structured cells; it does not
retain a second RPA generator.

## Discretization

ENGINE has max(8, ceil(exit_radius/baseCellSize)) radial cells; BAND has
max(4, ceil(band_thickness/baseCellSize)). Physical layers may increase these
global counts, never change them longitudinally. Automatic distributions have
a smooth 4:1 ratio between largest and smallest radial intervals, finer at walls.
No physical y+ is inferred.

Optional existing wall syntax remains: `innerWallLayers true`,
`innerWallFirstCellHeight`, `innerWallLayerCount`, `innerWallGrowthRate`,
`innerWallThickness`, and the corresponding `outer...` entries.
Growth is in [1,1.2]; height/count/growth/thickness must be consistent.
Prescribed normal offsets use the validated miter construction. Only those
near-wall nodes may move in x; core columns remain vertical. Incompatible
prescriptions fail rather than removing layers or enlarging BAND.

PLUME inherits ENGINE exit + physical lip + BAND exit, sharing node identities.
The lip uses one cell, including pure QUAD: it is not a transition boundary and
requires no extra parity split. Outlet rows form one growing geometric law.
Every row is interpolated index-to-index:
`r_j(x)=(1-t)*r_left,j+t*r_right,j`.
It need not be monotone in j-spacing on the inherited left boundary.

ENGINE/BAND columns advance using the arithmetic mean of **all** radial cell
heights, including physical wall layers:
`dx_target = (r_max - r_min) / radial_cell_count`.
At each mandatory CAD junction, throat or endpoint, widths are redistributed
over the interval rather than leaving a tiny final cell. Pure QUAD's even
CAD-curve counts are distributed over the same interval. No curvature cap,
chord-driven subdivision, shape-driven merging or quality-driven spacing
repair remains. Invalid columns are rejected without changing their widths.

This is an experimental spacing rule, not a geometric-accuracy guarantee.
The DXF import tolerance is unchanged, but linear wall faces can exceed it.
The preview reports the maximum vertical wall-chord error separately for
ENGINE/BAND and explicitly flags exceedances. A topologically valid diagnostic
mesh is not thereby certified against the wall-approximation tolerance.
Thin prescribed layers can fit the region heights yet be incompatible with
coarse linear wall chords: cell centres can then lie in the CAD shell. The
post-build preview/domain diagnostics reject that case; RawMeshData's angular
screen alone is not a CAD-domain certificate. No distribution is repaired to
hide it. The extended tests retain this combination as an expected rejection.

PLUME has one progressive x distribution. FARFIELD has independent uniform x
and growing radial spacing; its nominal cells are coarser. TRANSITION absorbs
different boundary counts. Generic interface-size extension provides smooth
unstructured sizes without special throat/exit/plume refinement envelopes.
FARFIELD's uniform x spacing targets `sqrt(dr_min * dr_max)` using its existing
radial heights; it is not globally capped by the last PLUME interval. For pure
QUAD, a bounded, aspect-sensitive relocation near rejected TRANSITION cells
can improve the outlet recombination. All structured and physical-boundary
nodes remain pinned, and the actual angular threshold is checked afterwards.
The unequal prescribed edge lengths at the acute outlet still limit its shape
quality: coarser FARFIELD spacing does not preserve every old interface metric.

Optional `externalSizeField` uses a plain whitespace-separated Cartesian grid,
with **no magic/version token**:

```text
nx nr xmin xmax rmin rmax
h_0_0 h_1_0 ... h_(nx-1)_0
h_0_1 h_1_1 ... h_(nx-1)_1
...
```

`nx` and `nr` are integers at least 2; coordinates and sizes are in metres.
Extents must be finite and increasing. Exactly `nx*nr` finite, strictly positive
sizes follow, with x varying fastest. Extra tokens are rejected. Sampling is
bilinear inside the grid and constantly extended from its edges outside.
For example, `2 2 0 1 0 1` followed by `0.01 0.02 0.03 0.04` is valid;
`CFD_SIZE_GRID_V1` is not a supported header. Paths are case-relative.
The file and callable `MeshSizeField` API constrain TRANSITION,
not ENGINE/BAND's exclusive mean-height rule; fixed radial counts cannot promise arbitrary local
isotropic refinement. `targetCellBudget` is an optional hard ceiling: it never
relaxes geometry, layers or quality and reports an incompatible budget.

## Quality and output

Geometric validity and quality are separate: every corner must have positive
orientation, areas must be positive, and interfaces must be conforming.
`minCellQuality` defaults to **0.05**, in (0,1]. It controls the **minimum signed
corner sine for both TRI and QUAD** (scaled corner Jacobian). A positively
oriented 1/89/90-degree triangle has quality about 0.01745, not 1.
Orthogonal QUAD anisotropy leaves its angles unchanged, whereas a thin TRI has
small angles and can fail the threshold. The aspect-sensitive `cell_quality`
exported by Mesh is a separate shape measure. Neither indicator alone guarantees
CFD accuracy; non-orthogonality, skewness and aspect ratios are reported separately.
The demo requests 0.05. The actual generated cells are checked;
a failed threshold is an error, never an automatic threshold reduction.
With plume factor 3 and downstream factor 20, the outlet wedge angle is about
4.23 degrees (sine 0.0738): threshold 0.10 is incompatible, including TRI/HYBRID
under the minimum-angle contract. Tests use explicitly smaller thresholds for
such configurations or expect a quality rejection. The demo's plume factor 6
opens this wedge; the demo uses HYBRID at 0.01 m with threshold 0.05.
A separate synthetic QUAD regression at 0.005 m retains threshold 0.10.
Before generation, every topology requires an outlet factor at least
`Fmin = Rin/Re + (L/Re)*q/sqrt(1-q*q)`, where `Re` is the inner exit radius,
`Rin` the plume inlet radius, `L` its downstream length and `q = minCellQuality`.
At `q=1`, no finite outlet factor is admissible. The minimum radius must also
remain strictly below the outer box radius. Incompatible inputs report the
requested/minimum factors, downstream factor and quality threshold; they are
never silently adjusted. For the default 80/65 box and `q=0.05`, `Fmin` is
about 5.5250, so factor 6 clears this necessary bound; factor 3
is rejected. At downstream factor 20, factor 3 remains admissible for `q=0.05`,
but not `q=0.10`. Clearing this tip bound does not guarantee quality throughout
TRANSITION: the actual generated-cell checks remain mandatory.
The default factor 16 was selected from 6/7/8/10/12/16 using the reference CAD,
the 80/65 box and `q=0.05`, without changing any mesh distribution. It is the
only tested factor passing both QUAD and HYBRID at 0.020 m and 0.005 m.
Minimum corner sines were respectively 0.17811/0.06120 and 0.08847/0.07026
(QUAD/HYBRID). This measured margin selects a default, not a new user rejection
threshold or a guarantee for other resolutions/profiles. Smaller boxes may
require an explicit outlet factor strictly below their outer radius.
TRI/HYBRID may also be rejected by the actual generated angles;
topological validity alone does not imply acceptable angular quality.

The eight physical groups remain `inlet`, `axis`, `inner_wall`, `outer_wall`,
`nozzle_lip`, `farfield_upstream`, `farfield_radial`, `farfield_downstream`.
Exit interfaces between fluids are internal, while the lip remains physical.
The usual preview exports a VTU with `rpa_region` and existing quality arrays:

```sh
./build-release/cfd_mesh_preview cases/rpaNozzle
```

The layout/diagnostics reject cell centers in the fictitious shell.
RawMeshData, Mesh topology, finite-volume geometry and CFD solvers are unchanged.

## Targeted tests

With CAD dependencies available, `rpa_nozzle` runs reader, parameter, layout,
structured-grid and metric tests without Gmsh generation. `rpa_nozzle_integration`
runs a restricted full-generation subset (CTest label `integration`). The extended
one-factor parameter matrix is opt-in, not a Cartesian product:

```sh
./build-release/cfd_rpa_nozzle_tests --parametric
# Or configure -DBUILD_RPA_PARAMETRIC_TESTS=ON, then:
ctest --test-dir build-release -L parametric --output-on-failure
```

Each configuration uses an isolated temporary case and the normal meshDict
reader/generator. The reference dictionary and DXF are never modified.
The `reference` case reads the tracked demo dictionary and checks its size,
cell type, Hybrid mode and quality threshold against explicit expectations.
Synthetic QUAD cases, including `quad_005`, are independent of that dictionary.

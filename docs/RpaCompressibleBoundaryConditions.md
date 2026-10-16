# RPA boundary data: preparation, not a reacting simulation

## Provenance

The four original user documents were found and read, without moving them,
in `MECH6251_Project_GENET_JACOB_GUILLAUMOT` in the user's technical library:
`Rapport_propulsion_GOx-GCH4.pdf`, `04_Results_RPA.pdf`,
`05_Input_NASA-CEA.txt`, `06_Output_NASA-CEA.txt`.
The entries below are verified against those files, not solely copied from
the task specification. Printed page numbers refer to the propulsion report.

| Context / quantity | Value / units | Static/total; absolute/gauge | Source | Euler BC / availability |
|---|---|---|---|---|
| Fluent CH4 inlet | 0.1032 kg/s; T0=300 K | inward group mass flow; total T | report §2.8.3 p26, table 2.11 | MassFlowInlet for a **single gas**; actual CH4 species needs future model |
| Fluent O2 inlet | 0.3198 kg/s; T0=300 K | inward group mass flow; total T | same | same limitation |
| Both inlet initial pressures | 1.44 MPa | initial gauge, **not imposed P0** | table 2.11 | initialization information, not an extra mass-inlet constraint |
| Fluent operating pressure | 1 Pa | reference absolute pressure | table 2.10 p26 | conversion convention, not gas inlet data |
| Four Fluent outlets | 1 Pa gauge | static, giving **2 Pa absolute** if `pabs=pgauge+poperating` | table 2.12 p26 | PressureOutlet only on supported outgoing branch; ignored when normally supersonic |
| Fluent backflow | 2.7 K; species “none” | backflow T, static/total not established by table | table 2.12 | not mapped automatically; current outlet rejects inflow |
| Near-space environment | 0.0001 Pa | ambient absolute | report §2.6.1, eq. 2.31 | distinct environment; not forced on supersonic faces |
| Optimum-expansion environment | 1000 Pa | ambient absolute | report §3.3.2 | distinct from Fluent's 2 Pa and near-vacuum setting |
| RPA injector/chamber pressure | 1.2000 MPa | static at zero Mach injector | RPA table 2 | not automatically nozzle inlet P0 |
| RPA nozzle inlet | p=1.1700 MPa; T=3343.5658 K; M=0.1508 | static absolute p/T, total Mach | RPA table 2 | approximate equivalent-gas total inlet possible |
| Inlet gas properties | gamma=1.1770; R=414.2 J/(kg K); rho=0.8448 kg/m³ | local RPA properties | RPA table 2; R converted from kJ/(kg K) | freeze gamma/R only as an approximation |
| RPA throat / exit Mach | 1.0000 / 4.5630 | local total Mach | RPA table 2 | normal regime still depends on face orientation |
| Nominal total flow | 0.4230 kg/s | CH4+O2 | RPA performance block after table 8; report table 2.7 | future inlet-group target, not imposed together with P0 |
| O/F | 3.1 | mass ratio | RPA table 1; CEA input/output | future multicomponent model |
| CEA setup | 12 bar; area ratio=100; equilibrium; CH4/O2 | chamber pressure / expansion specification | CEA input and echoed output | not a constant-gamma reference solution |

For a calorically perfect equivalent hot gas, use gamma=1.1770, R=414.2
and **derive** `T0=T*(1+(gamma-1)*M²/2)` and
`P0=p*(T0/T)^(gamma/(gamma-1))`. These reconstructed values are not original
RPA total-condition outputs. Recompute density from `p/(R*T)`; the rounded
RPA density is not an additional independent constraint. RPA uses variable
thermochemistry (even its tabulated isentropic exponent differs from gamma).
Constant properties cannot reproduce its whole nozzle expansion.

The demonstration uses this already-hot equivalent gas in a straight tube,
not cold CH4 and O2 mixed to fabricate combustion products. It transports only
rho, rho*ux, rho*ur, rhoE: no species, reactions, flame or combustion heating.

## Actual geometry groups

The existing RPA generator defines exactly `inlet`, `axis`, `inner_wall`,
`outer_wall`, `nozzle_lip`, `farfield_upstream`, `farfield_radial`,
`farfield_downstream` (see the existing RPA documentation/generator).
There are **not** separate CH4/O2 injector groups, nor a one-to-one set of
the four Fluent pressure outlets. No fictitious groups are added.

An eventual hot-gas nozzle calculation could use the existing `inlet` with
total conditions or a compatible group mass inlet, `axis` with Axis, and the
three solid groups with SlipWall. Ambient farfield boundaries can involve
subsonic inflow, backflow or mixed normal regimes: a general reservoir/farfield
model is still needed before assigning these new outlet types indiscriminately.
The complete RPA geometry is **not** simulated in lot E. Injector geometry,
species transport, reacting thermodynamics, viscous walls and appropriate
ambient inflow data are prerequisites for the report's actual Fluent problem.

#!/usr/bin/env python3
"""CAD-only adapter: DXF -> ordered meridian polyline in metres (no meshing).

stdout protocol: CFD_RPA_PROFILE_V1, point count, x r pairs, then
CURVE_JUNCTIONS, count and point indices. Failure diagnostics use stdout
with nonzero exit status so the C++ caller can retain the cause.
Curve flattening uses the supplied geometric tolerance, independently of mesh size.
"""

import argparse
import math
import sys


def import_profile(filename, tolerance, units_override=None, return_junctions=False):
    try:
        import ezdxf
        from ezdxf import units
    except ImportError as error:
        raise ValueError("DXF import requires ezdxf: install tools/requirements-cad.txt in the configured Python") from error

    doc = ezdxf.readfile(filename)
    unit_names = {"mm": units.MM, "cm": units.CM, "m": units.M, "in": units.IN}
    declared = int(doc.header.get("$INSUNITS", 0))
    if units_override:
        specified = unit_names[units_override]
        if declared and specified != declared:
            raise ValueError("dxfUnits conflicts with DXF $INSUNITS")
        declared = specified
    if not declared:
        raise ValueError("Absent/ambiguous DXF units; specify geometry dxfUnits")
    try:
        scale = units.conversion_factor(declared, units.M)
    except (ValueError, TypeError, KeyError, IndexError) as error:
        raise ValueError("Unsupported/ambiguous DXF units") from error
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid DXF unit conversion")

    annotations = {"DIMENSION", "TEXT", "MTEXT", "ATTRIB", "ATTDEF", "LEADER", "MLEADER",
                   "POINT", "HATCH", "SOLID", "VIEWPORT"}

    def expand(entities, depth=0):
        if depth > 32:
            raise ValueError("Recursive DXF block nesting exceeds 32")
        for entity in entities:
            kind = entity.dxftype()
            if kind in annotations:
                continue
            if kind == "INSERT":
                def skipped(child, reason):
                    if child.dxftype() not in annotations:
                        raise ValueError(f"Cannot transform required DXF entity {child.dxftype()}: {reason}")
                inserts = entity.multi_insert() if entity.mcount > 1 else (entity,)
                for insert in inserts:
                    yield from expand(insert.virtual_entities(skipped_entity_callback=skipped), depth + 1)
            elif kind in {"POLYLINE", "LWPOLYLINE"}:
                if kind == "POLYLINE" and not entity.is_2d_polyline:
                    raise ValueError("Only planar DXF polylines are supported")
                yield from expand(entity.virtual_entities(), depth + 1)
            else:
                yield entity

    entities = list(expand(doc.modelspace()))
    axes = [e for e in entities if e.dxftype() == "LINE" and
            (e.dxf.linetype.upper().startswith("CENTER") or e.dxf.layer.lower() in {"axis", "centerline"})]
    origin = (0.0, 0.0)
    direction = (1.0, 0.0)
    if axes:
        if len(axes) != 1:
            raise ValueError("Ambiguous DXF axis: expected one centreline")
        start, end = axes[0].dxf.start, axes[0].dxf.end
        dx, dy = end.x - start.x, end.y - start.y
        length = math.hypot(dx, dy)
        if not length > 0:
            raise ValueError("Degenerate DXF axis")
        origin = (start.x, start.y)
        direction = (dx / length, dy / length)

    def project(point):
        if abs(point.z * scale) > tolerance:
            raise ValueError("Non-planar DXF contour")
        dx, dy = point.x - origin[0], point.y - origin[1]
        result = ((dx * direction[0] + dy * direction[1]) * scale,
                  (-dx * direction[1] + dy * direction[0]) * scale)
        if not all(math.isfinite(v) for v in result):
            raise ValueError("Non-finite DXF coordinate")
        return result

    chains = []
    lower_chains = []
    for entity in entities:
        kind = entity.dxftype()
        if entity in axes:
            continue
        if kind == "LINE":
            points = [entity.dxf.start, entity.dxf.end]
        elif kind in {"ARC", "ELLIPSE"}:
            # ELLIPSE can arise from non-uniform scaling of a block ARC.
            points = list(entity.flattening(tolerance / scale))
        elif kind == "SPLINE":
            points = list(entity.construction_tool().flattening(tolerance / scale, segments=8))
        else:
            raise ValueError(f"Unsupported required DXF entity: {kind}")
        chain = [project(p) for p in points]
        if len(chain) < 2:
            raise ValueError("Degenerate DXF curve")
        if all(abs(p[1]) <= tolerance for p in chain):
            continue  # A plain axis line is not a wall.
        if all(p[1] < -tolerance for p in chain):
            lower_chains.append(chain)
            continue  # Optional mirrored lower meridian.
        if any(p[1] <= 0 for p in chain):
            raise ValueError("DXF wall radius must be strictly positive")
        chains.append(chain)
    if not chains and lower_chains:
        # A reflected block can reverse the normal of its transformed axis.
        chains = [[(x, -radius) for x, radius in chain] for chain in lower_chains]
    if not chains:
        raise ValueError("DXF contains no positive-radius meridian")

    # Endpoint graph: a unique open chain, not file order or entity handles.
    endpoints = []
    connections = []
    join_tolerance = 2 * tolerance
    def vertex(point):
        for index, existing in enumerate(endpoints):
            if math.dist(point, existing) <= join_tolerance:
                return index
        endpoints.append(point)
        connections.append([])
        return len(endpoints) - 1
    edges = []
    for index, chain in enumerate(chains):
        first, last = vertex(chain[0]), vertex(chain[-1])
        if first == last:
            raise ValueError("Closed/self-intersecting wall contour is not r(x)")
        edges.append((first, last))
        connections[first].append(index)
        connections[last].append(index)
    ends = [i for i, attached in enumerate(connections) if len(attached) == 1]
    if len(ends) != 2 or any(len(attached) not in (1, 2) for attached in connections):
        raise ValueError("Discontinuous, branched or ambiguous DXF contour")
    current = min(ends, key=lambda i: endpoints[i][0])
    used = set()
    ordered = []
    junctions = [0]
    while True:
        available = [i for i in connections[current] if i not in used]
        if not available:
            break
        index = available[0]
        first, last = edges[index]
        chain = chains[index] if first == current else list(reversed(chains[index]))
        ordered.extend(chain if not ordered else chain[1:])
        junctions.append(len(ordered) - 1)
        current = last if first == current else first
        used.add(index)
    if len(used) != len(chains):
        raise ValueError("Disconnected DXF wall contour")
    if any(b[0] <= a[0] for a, b in zip(ordered, ordered[1:])):
        # Strict x-monotonicity also excludes crossings of the flattened wall.
        raise ValueError("Self-intersection/backtracking: wall is not a single-valued r(x)")
    throat = min(range(len(ordered)), key=lambda i: ordered[i][1])
    if throat == 0 or throat == len(ordered) - 1:
        raise ValueError("DXF must contain an interior throat followed by a divergent nozzle")
    return (ordered, junctions) if return_junctions else ordered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file")
    parser.add_argument("tolerance", type=float)
    parser.add_argument("units", nargs="?", default="auto", choices=("auto", "mm", "cm", "m", "in"))
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        raise ValueError("Geometric tolerance must be finite and positive")
    profile, junctions = import_profile(args.file, args.tolerance, None if args.units == "auto" else args.units, True)
    print("CFD_RPA_PROFILE_V1")
    print(len(profile))
    for x, radius in profile:
        print(f"{x:.17g} {radius:.17g}")
    print("CURVE_JUNCTIONS", len(junctions), *junctions)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"RPA DXF import: {error}")
        sys.exit(1)

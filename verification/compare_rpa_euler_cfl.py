#!/usr/bin/env python3
"""Compare short Euler RPA runs, using conservative checkpoints, not rho*u from VTU.

The first directory is a temporal comparison reference, not an exact solution.
Norms are volume-weighted and separately normalized in Engine and exterior.
"""
import csv
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from validate_rpa_euler_output import validated_polygons


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def settings(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def conservative_values(path):
    data = path.read_bytes()
    if data.splitlines()[0] != b"CFD_EULER_CHECKPOINT 1":
        raise ValueError("Unsupported conservative checkpoint version")
    offset = data.rindex(b"FNV1A64 ")
    checksum = 14695981039346656037
    for byte in data[:offset]:
        checksum = ((checksum ^ byte) * 1099511628211) & ((1 << 64) - 1)
    if data[offset:].split() != [b"FNV1A64", str(checksum).encode()]:
        raise ValueError("Checkpoint integrity failure")
    line_end = data.index(b"\n") + 1
    metadata_end = data.index(b"\n", line_end)
    size = int(data[line_end:metadata_end])
    metadata = data[metadata_end + 1:metadata_end + 1 + size].decode()
    lines = metadata.splitlines()
    numerical_options = lines[0].split()
    del numerical_options[2]  # Only CFL may differ in this campaign.
    compatibility = " ".join(numerical_options) + "\n" + "\n".join(lines[1:])
    tokens = data[metadata_end + 1 + size:offset].split()
    count = int(tokens[22])  # Version 1 accepted-step diagnostics precede the field.
    values = list(map(float, tokens[23:23 + 4 * count]))
    if len(values) != 4 * count or not all(map(math.isfinite, values)):
        raise ValueError("Invalid conservative field")
    return {name: values[k::4] for k, name in enumerate(("rho", "rho_ux", "rho_ur", "rho_E"))}, compatibility


def load(root):
    cfg, status = settings(root / "assumptions.txt"), settings(root / "status.txt")
    piece = ET.parse(root / f"rpa_{status['steps']}.vtu").find(".//Piece")
    fields = {a.attrib["Name"]: list(map(float, a.text.split())) for a in piece.find("CellData")}
    xyz = list(map(float, piece.find("Points")[0].text.split()))
    topology = {a.attrib["Name"]: list(map(int, a.text.split())) for a in piece.find("Cells")}
    volumes = []
    for ids in validated_polygons(topology, len(fields["rho"]), len(xyz)//3):
        xy = [(xyz[3 * n], xyz[3 * n + 1]) for n in ids]
        volumes.append(math.pi / 3 * math.fsum(
            (xy[k][1] + xy[(k + 1) % len(ids)][1]) *
            (xy[k][0] * xy[(k + 1) % len(ids)][1] - xy[(k + 1) % len(ids)][0] * xy[k][1]) for k in range(len(ids))))
    suffix = "" if int(status["rejected_steps"]) == 0 else "_r" + status["rejected_steps"]
    conservatives, compatibility = conservative_values(root / f"checkpoint_{status['steps']}{suffix}.chk")
    fields.update(conservatives)
    return {"root": str(root), "configuration": cfg, "status": status, "fields": fields,
            "volumes": volumes, "mesh": [xyz, topology], "sections": table(root / "sections.csv")[-3:],
            "cell_diagnostics": table(root / "cell_diagnostics.csv")[-1],
            "convergence": table(root / "convergence.csv")[-1], "compatibility": compatibility}


def relative(numerator, denominator):
    return numerator / denominator if denominator else (0.0 if numerator == 0 else math.inf)


def compare(reference, run):
    if reference["mesh"] != run["mesh"] or reference["volumes"] != run["volumes"]:
        raise ValueError("Meshes differ")
    if reference["status"]["time_seconds"] != run["status"]["time_seconds"]:
        raise ValueError("Physical comparison times differ")
    if reference["compatibility"] != run["compatibility"]:
        raise ValueError("A numerical/physical option other than CFL differs")
    for key in ("step_cap", "maximum_time", "maximum_wall_seconds", "maximum_rejected", "diag_every", "save_every"):
        if reference["configuration"][key] != run["configuration"][key]:
            raise ValueError("Different campaign limit/schedule: " + key)
    errors = {}
    for region in ("Engine", "exterior"):
        ids = [c for c, label in enumerate(reference["fields"]["rpa_region"])
               if (label == 1) == (region == "Engine")]
        volume = math.fsum(reference["volumes"][c] for c in ids)
        errors[region] = {}
        for field in ("rho", "rho_ux", "rho_ur", "rho_E", "p", "Mach"):
            a, b, v = reference["fields"][field], run["fields"][field], reference["volumes"]
            l1 = math.fsum(v[c] * abs(b[c] - a[c]) for c in ids)
            l2 = math.sqrt(math.fsum(v[c] * (b[c] - a[c]) ** 2 for c in ids))
            scale1 = math.fsum(v[c] * abs(a[c]) for c in ids)
            scale2 = math.sqrt(math.fsum(v[c] * a[c] ** 2 for c in ids))
            maximum = max(abs(b[c] - a[c]) for c in ids)
            errors[region][field] = {"L1_absolute_mean": l1 / volume, "L1_relative": relative(l1, scale1),
                                    "L2_absolute_rms": l2 / math.sqrt(volume), "L2_relative": relative(l2, scale2),
                                    "maximum_absolute": maximum,
                                    "maximum_scaled": relative(maximum, max(abs(a[c]) for c in ids))}
    return errors


def main():
    runs = [load(Path(name)) for name in sys.argv[1:]]
    if not runs:
        raise ValueError("Usage: compare_rpa_euler_cfl.py REFERENCE [CANDIDATE ...]")
    report = []
    print("CFL steps rejects loop_s M_col M_exit Q_col Q_exit Engine_kg conservation max_relative_L2")
    for run in runs:
        history = table(Path(run["root"]) / "steps.csv")
        times = {float(line.split("=", 1)[1]) for line in
                 (Path(run["root"]) / "assumptions.txt").read_text().splitlines()
                 if line.startswith("requested_output_time=")}
        dt = [float(r["dt"]) for r in history]
        natural = [float(r["dt"]) for r in history if float(r["time"]) not in times]
        errors = compare(runs[0], run)
        sections = {r["section"]: r for r in run["sections"]}
        compact = {k: v for k, v in run.items() if k not in ("fields", "volumes", "mesh", "compatibility")}
        compact.update(errors=errors, dt_mean=math.fsum(dt) / len(dt), dt_min_natural=min(natural), dt_max=max(dt),
                       seconds_per_step=float(run["status"]["loop_s"]) / len(dt))
        report.append(compact)
        print(run["configuration"]["CFL"], len(history), run["status"]["rejected_steps"], run["status"]["loop_s"],
              *[sections[n]["mean_cell_Mach"] for n in ("throat", "exit")],
              *[sections[n]["numerical_mass_flow"] for n in ("throat", "exit")],
              run["cell_diagnostics"]["engine_mass"],
              max(float(run["status"][f"balance_relative_{k}"]) for k in range(4)),
              max(e["L2_relative"] for region in errors.values() for e in region.values()))
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

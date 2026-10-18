#!/usr/bin/env python3
"""Independent output/ledger checks, not validation of a rarefied plume model."""
import csv
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET


def settings(path):
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def require_time_increment(previous, current, step):
    # The 17-digit CSV round-trips binary64. Check the actual accumulation;
    # subtracting two growing times amplifies rounding relative to a tiny dt.
    require(step > 0 and current > previous and current == previous + step, "Time increment mismatch")


def require_history_endpoint(history, time, initial_time=0.0):
    require(time == (float(history[-1]["time"]) if history else initial_time),
            "Last snapshot/status time does not match the accepted history endpoint")


def validate(root):
    config = settings(root / "assumptions.txt")
    status = settings(root / "status.txt")
    gamma, gas = float(config["gamma"]), float(config["R"])
    initial_time = float(config.get("segment_initial_time", "0"))
    initial_step = int(config.get("segment_initial_step", "0"))
    series = ET.parse(root / "rpa.pvd").findall(".//DataSet")
    require(series, "Empty PVD")
    times = [float(item.attrib["timestep"]) for item in series]
    require(times[0] == initial_time and all(math.isfinite(t) for t in times), "Invalid initial/time values")
    require(all(b > a for a, b in zip(times, times[1:])), "PVD chronology")
    with (root / "steps.csv").open() as stream:
        history = list(csv.DictReader(stream))
    require(len(history) == int(status["steps"]) - initial_step, "Accepted-step count mismatch")
    for step, row in enumerate(history, 1):
        require(all(math.isfinite(float(value)) for value in row.values()), "Nonfinite history")
        require(int(row["step"]) == initial_step + step and float(row["dt"]) > 0, "Step numbering/dt")
        require(float(row["minimum_density"]) > 0 and float(row["minimum_pressure"]) > 0, "Accepted positivity")
        previous = float(history[step-2]["time"]) if step > 1 else initial_time
        require_time_increment(previous, float(row["time"]), float(row["dt"]))
    require(times[-1] == float(status["time_seconds"]), "Last snapshot is not last accepted state")
    require_history_endpoint(history, times[-1], initial_time)
    accepted_times = {float(row["time"]) for row in history}
    require(all(t == initial_time or t in accepted_times for t in times), "Snapshot is not an accepted-step time")
    if config.get("physical_exports") == "1":
        requested = [float(line.split("=",1)[1]) for line in (root/"assumptions.txt").read_text().splitlines()
                     if line.startswith("requested_output_time=")]
        require(0 <= len(requested) <= 6 and all(t > initial_time for t in requested) and
                all(b > a for a,b in zip(requested,requested[1:])), "Output schedule")
        require(times[1:-1] == requested[:max(0,len(times)-2)], "Physical export prefix mismatch")
        require(times[-1] <= requested[min(len(times)-2,len(requested)-1)] if len(times) > 1 and requested else True, "Skipped physical export")
        if status["status"] == "TARGET_REACHED":
            require(times == [initial_time]+requested and times[-1] == float(config["maximum_time"]), "Incomplete target run")
        require((status["target_reached"] == "1") == (times[-1] >= float(config["maximum_time"])), "Target status")
    section_rows = {}
    cell_rows = {}
    if (root/"cell_diagnostics.csv").exists():
        with (root/"cell_diagnostics.csv").open() as stream:
            cell_rows = {float(row["time"]):row for row in csv.DictReader(stream)}
    if (root/"sections.csv").exists():
        with (root/"sections.csv").open() as stream:
            for row in csv.DictReader(stream):
                key = (float(row["time"]),row["section"])
                require(key not in section_rows, "Duplicate section observation")
                section_rows[key] = row
    reference_mesh = None
    initial_integrals = None
    initial_pressure = None
    regions_reference = None
    for entry in series:
        filename = Path(entry.attrib["file"])
        require(not filename.is_absolute() and ".." not in filename.parts, "Unsafe snapshot path")
        piece = ET.parse(root / filename).find(".//Piece")
        require(piece is not None, "Missing VTU Piece")
        nc, np = int(piece.attrib["NumberOfCells"]), int(piece.attrib["NumberOfPoints"])
        require(nc == int(config["cells"]), "Cell count mismatch")
        fields = {a.attrib["Name"]: list(map(float, a.text.split())) for a in piece.find("CellData")}
        names = ["rho", "p", "T", "ux", "ur", "velocity_magnitude", "Mach", "rhoE", "log10_p_Pa", "log10_rho", "rpa_region"]
        require(all(n in fields and len(fields[n]) == nc and all(map(math.isfinite, fields[n])) for n in names), "Invalid fields")
        require(len(fields["U"]) == 3*nc and all(map(math.isfinite, fields["U"])), "Invalid velocity vector")
        points = list(map(float, piece.find("Points")[0].text.split()))
        cells = {a.attrib["Name"]: list(map(int, a.text.split())) for a in piece.find("Cells")}
        require(len(points) == 3*np and all(map(math.isfinite, points)), "Invalid points")
        require(min(points[1::3]) >= 0, "Negative radius")
        require(cells["types"] == [9]*nc and cells["offsets"] == list(range(4,4*nc+1,4)), "Expected QUAD connectivity")
        require(len(cells["connectivity"]) == 4*nc, "Connectivity size")
        mesh = (points, cells)
        if reference_mesh is None:
            reference_mesh = mesh
        require(mesh == reference_mesh, "Mesh changed between snapshots")
        require(all(r in (1,2,3,4,5) for r in fields["rpa_region"]), "Invalid region label")
        if regions_reference is None:
            regions_reference = fields["rpa_region"]
        require(fields["rpa_region"] == regions_reference, "Regions changed")
        if section_rows:
            # Independent reconstruction from VTU polygons. These are controlled
            # Engine-cell estimates on exact sections, not physical trace states.
            edges = {}
            for c in range(nc):
                ids = cells["connectivity"][4*c:4*c+4]
                for i in range(4):
                    edges.setdefault(tuple(sorted((ids[i],ids[(i+1)%4]))),[]).append(c)
            for name in ("inlet","throat","exit"):
                row = section_rows.get((float(entry.attrib["timestep"]),name))
                if row is None and (float(entry.attrib["timestep"]) == initial_time or
                                    (status["status"] in ("SEGMENT_STOP", "WALL_LIMIT") and
                                     float(entry.attrib["timestep"]) == times[-1])):
                    continue  # No synthetic monitor observation is inserted at a segment boundary.
                require(row is not None, "Missing section observation at exported time")
                x, radius = float(config[name+"_x"]),float(config[name+"_r"])
                rings = []
                mass_terms = []
                for (a,b), neighbors in edges.items():
                    if max(abs(points[3*a]-x),abs(points[3*b]-x)) < 1e-10 and max(points[3*a+1],points[3*b+1]) <= radius+1e-10:
                        inside = [c for c in neighbors if fields["rpa_region"][c] == 1]
                        if inside:
                            ring_area = math.pi*abs(points[3*a+1]**2-points[3*b+1]**2)
                            rings.append((ring_area,inside))
                            if name != "inlet":
                                require(len(neighbors) == 2, "Expected an internal throat/exit face")
                                left,right = sorted(neighbors,key=lambda c:math.fsum(points[3*i] for i in cells["connectivity"][4*c:4*c+4]))
                                rho_l,rho_r = fields["rho"][left],fields["rho"][right]
                                u_l,u_r = fields["ux"][left],fields["ux"][right]
                                wave = max(abs(u_l)+math.sqrt(gamma*fields["p"][left]/rho_l),
                                           abs(u_r)+math.sqrt(gamma*fields["p"][right]/rho_r))
                                mass_terms.append(ring_area*(.5*(rho_l*u_l+rho_r*u_r)-.5*wave*(rho_r-rho_l)))
                area = math.fsum(a for a,_ in rings)
                require(math.isclose(area,math.pi*radius**2,rel_tol=1e-11), "Analytical section coverage")
                require(int(row["faces"]) == len(rings), "Section face count")
                if name != "inlet":
                    require(math.isclose(float(row["numerical_mass_flow"]),math.fsum(mass_terms),rel_tol=1e-11),
                            "Independent internal numerical section flow mismatch")
                for field, column in [("Mach","mean_cell_Mach"),("p","mean_cell_pressure"),("T","mean_cell_temperature")]:
                    expected = math.fsum(a*math.fsum(fields[field][c] for c in ids)/len(ids) for a,ids in rings)/area
                    require(math.isclose(float(row[column]),expected,rel_tol=1e-10), "Independent section estimate mismatch")
        volumes = []
        for cell in range(nc):
            ids = cells["connectivity"][4*cell:4*cell+4]
            require(len(set(ids)) == 4 and all(0 <= node < np for node in ids), "Invalid cell node IDs")
            cross = [points[3*ids[i]]*points[3*ids[(i+1)%4]+1]-points[3*ids[(i+1)%4]]*points[3*ids[i]+1] for i in range(4)]
            require(math.fsum(cross) > 0, "Inverted cell")
            volume = math.pi/3*math.fsum((points[3*ids[i]+1]+points[3*ids[(i+1)%4]+1])*cross[i] for i in range(4))
            require(volume > 0, "Invalid revolution volume")
            volumes.append(volume)
            rho,p,t,u,v,speed,mach,energy,logp,logrho,_ = [fields[n][cell] for n in names]
            require(min(rho,p,t) > 0, "Snapshot positivity")
            expected = [p/(rho*gas),math.hypot(u,v),math.hypot(u,v)/math.sqrt(gamma*p/rho),p/(gamma-1)+.5*rho*(u*u+v*v),math.log10(p),math.log10(rho)]
            for actual, target in zip([t,speed,mach,energy,logp,logrho],expected):
                require(math.isclose(actual,target,rel_tol=1e-10,abs_tol=1e-13 if target == 0 else 0), "EOS/derived field mismatch")
            require(fields["U"][3*cell:3*cell+3] == [u,v,0], "Vector/component mismatch")
        integrals = [math.fsum(volumes[c]*value for c,value in enumerate(component)) for component in
                     [fields["rho"], [r*u for r,u in zip(fields["rho"],fields["ux"])],
                      [r*v for r,v in zip(fields["rho"],fields["ur"])], fields["rhoE"]]]
        cell_row = cell_rows.get(float(entry.attrib["timestep"]))
        if cell_row is not None and "engine_mass" in cell_row:
            engine_mass = math.fsum(volumes[c]*fields["rho"][c] for c in range(nc) if fields["rpa_region"][c] == 1)
            require(math.isclose(float(cell_row["engine_mass"]),engine_mass,rel_tol=1e-12), "Engine mass mismatch")
        if initial_integrals is None:
            initial_integrals, initial_pressure = integrals, fields["p"]
            hot_count = sum(r == 1 for r in fields["rpa_region"])
            require(hot_count == int(config["hot_cells"]), "Engine partition mismatch")
            hot_p, hot_t = float(config["p_inlet_static"]), float(config["T_inlet_static"])
            hot_u = float(config["Mach_inlet"])*math.sqrt(gamma*gas*hot_t)
            for c, region in enumerate(fields["rpa_region"] if initial_time == 0 else []):
                expected = ([hot_p, hot_t, hot_u, 0] if region == 1 else
                            [float(config["p_ambient_nominal"]), float(config["T_ambient_ASSUMED"]),
                             float(config["ux_ambient_ASSUMED"]), float(config["ur_ambient_ASSUMED"])])
                for actual, target in zip([fields[n][c] for n in ("p","T","ux","ur")], expected):
                    require(math.isclose(actual,target,rel_tol=1e-12,abs_tol=0), "Initial region/state mismatch")
        last_integrals = integrals
    if history:
        boundary = ["mass_outflow_rate","boundary_mx","boundary_mr","boundary_energy"]
        for k in range(4):
            rates = [float(r["dt"])*(float(r[boundary[k]])-(float(r["radial_source"]) if k == 2 else 0)) for r in history]
            defect = last_integrals[k]-initial_integrals[k]+math.fsum(rates)
            scale = abs(last_integrals[k])+abs(initial_integrals[k])+math.fsum(map(abs,rates))
            require(abs(defect) <= 512*sys.float_info.epsilon*(len(history)+1)*scale, "Independent integrated conservation failure")
            print(f"balance[{k}]={defect:.9g}; scaled={abs(defect)/scale:.3g}")
    require(status["conserved"] == "1", "Driver conservation gate failed")
    if len(history) >= 200:
        require(max(abs(a-b) for a,b in zip(initial_pressure,fields["p"])) > 1, "No nontrivial solved evolution")
    print(f"PASS: {len(series)} snapshots; {nc} QUAD; {len(history)} accepted steps; status={status['status']}")
    if initial_step:
        print("Restart segment: independent ledger covers this segment only; checkpoint preserves the cumulative ledger.")
    print("Structural/EOS/ledger checks only; not a physical validation or proof of stationarity.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 verification/validate_rpa_euler_output.py OUTPUT")
    try:
        validate(Path(sys.argv[1]))
    except (ValueError,OSError,KeyError,ET.ParseError) as error:
        sys.exit(f"FAIL: {error}")

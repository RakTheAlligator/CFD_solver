"""Deterministic CAD fixtures/tests; this module never generates a mesh."""
import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
import ezdxf

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rpa_cad", ROOT / "tools/import_rpa_dxf.py")
cad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cad)


def make_document(points=None, cubic=False, block=False, units=6):
    doc = ezdxf.new("R2000")
    doc.units = units
    layout = doc.blocks.new("MERIDIAN") if block else doc.modelspace()
    if cubic:
        layout.add_line((0, 0.3), (0.5, 0.3))
        layout.add_line((0.5, 0.3), (1, 0.1))
        layout.add_open_spline([(1, 0.1), (1.2, 0.1), (1.6, 0.35), (2, 0.45)], degree=3)
    else:
        layout.add_lwpolyline(points or [(0, 0.3), (0.5, 0.3), (1, 0.1), (2, 0.4)])
    layout.add_line((0, 0), (2, 0), dxfattribs={"linetype": "CENTER"})
    layout.add_text("ignore annotations", dxfattribs={"height": 0.05})
    if block:
        doc.modelspace().add_blockref("MERIDIAN", (3, 4), dxfattribs={"rotation": 23, "xscale": 2, "yscale": 2})
    return doc


def write_fixtures(directory):
    directory.mkdir(parents=True, exist_ok=True)
    make_document().saveas(directory / "conical.dxf")
    make_document(cubic=True).saveas(directory / "cubic.dxf")


class ImportTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="cfd-rpa-cad-")
        self.file = Path(self.directory.name) / "test.dxf"

    def tearDown(self):
        self.directory.cleanup()

    def profile(self, doc, override=None):
        doc.saveas(self.file)
        return cad.import_profile(self.file, 1e-6, override)

    def test_real_rpa(self):
        points = cad.import_profile(ROOT / "cases/rpaNozzle/constant/geometry/engine.dxf", 1e-6)
        throat = min(points, key=lambda p: p[1])
        exit = points[-1]
        self.assertAlmostEqual(2 * throat[1], 0.02858, delta=3e-6)
        self.assertAlmostEqual(2 * exit[1], 0.28582, delta=5e-6)
        self.assertAlmostEqual((exit[1] / throat[1]) ** 2, 100, delta=0.01)
        self.assertAlmostEqual(exit[0] - throat[0], 0.38458, delta=3e-6)

    def test_conical_polyline_and_cubic_spline(self):
        for doc in (make_document(), make_document(cubic=True)):
            points = self.profile(doc)
            self.assertAlmostEqual(min(p[1] for p in points), 0.1)
            self.assertTrue(all(a[0] < b[0] for a, b in zip(points, points[1:])))

    def test_polyline_and_bulge(self):
        doc = ezdxf.new()
        doc.units = 6
        layout = doc.modelspace()
        layout.add_polyline2d([(0, 0.3), (0.5, 0.3), (1, 0.1), (2, 0.4)])
        self.assertAlmostEqual(min(p[1] for p in self.profile(doc)), 0.1)
        doc = make_document()
        polyline = doc.modelspace().query("LWPOLYLINE")[0]
        polyline.set_points([(0, 0.3, 0), (0.5, 0.3, -0.05), (1, 0.1, 0), (2, 0.4, 0)], format="xyb")
        self.assertGreater(len(self.profile(doc)), 4)

    def test_insert_transform(self):
        points = self.profile(make_document(cubic=True, block=True))
        self.assertAlmostEqual(points[0][0], 0, delta=1e-12)
        self.assertAlmostEqual(points[0][1], 0.6)
        self.assertAlmostEqual(points[-1][0], 4)
        self.assertAlmostEqual(points[-1][1], 0.9)

    def test_units(self):
        with self.assertRaisesRegex(ValueError, "units"):
            self.profile(make_document(units=0))
        self.assertAlmostEqual(self.profile(make_document(units=0), "mm")[-1][0], 0.002)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.profile(make_document(), "mm")

    def test_nonuniform_insert_arc(self):
        doc = ezdxf.new()
        doc.units = 6
        block = doc.blocks.new("ARC_NOZZLE")
        block.add_line((0, 0.3), (1, 0.1))
        block.add_arc((1, 0.3), 0.2, 270, 315)
        block.add_line((1 + 0.2 / math.sqrt(2), 0.3 - 0.2 / math.sqrt(2)), (2, 0.4))
        block.add_line((0, 0), (2, 0), dxfattribs={"linetype": "CENTER"})
        insert = doc.modelspace().add_blockref("ARC_NOZZLE", (3, 4),
            dxfattribs={"rotation": 23, "xscale": 2, "yscale": 1.5})
        points = self.profile(doc)
        self.assertAlmostEqual(min(p[1] for p in points), 0.15)
        self.assertAlmostEqual(points[-1][0], 4)
        self.assertAlmostEqual(points[-1][1], 0.6)
        insert.dxf.xscale = -2
        reflected = self.profile(doc)
        self.assertAlmostEqual(min(p[1] for p in reflected), 0.15)
        self.assertAlmostEqual(reflected[-1][0], 4)
        self.assertAlmostEqual(reflected[-1][1], 0.6)

    def test_bad_contours(self):
        for points in ([(0, 0.3), (1, 0.1), (0.5, 0.4), (2, 0.4)],
                       [(0, 0.3), (1, 0), (2, 0.4)],
                       [(0, 0.3), (2, 0.4)]):
            with self.assertRaises(ValueError):
                self.profile(make_document(points=points))
        doc = make_document()
        doc.modelspace().add_line((3, 0.5), (4, 0.6))
        with self.assertRaisesRegex(ValueError, "Discontinuous|Disconnected"):
            self.profile(doc)
        doc = make_document()
        doc.modelspace().add_circle((3, 0.5), 0.1)
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            self.profile(doc)

    def test_missing_file(self):
        with self.assertRaises(OSError):
            cad.import_profile(self.file, 1e-6)

    def test_flattening_accuracy(self):
        doc = make_document(cubic=True)
        spline = doc.modelspace().query("SPLINE")[0].construction_tool()
        points = self.profile(doc)
        cubic = [p for p in points if p[0] >= 1]
        # Dense independent spline evaluation checks the chord approximation.
        for t in (i / 1000 for i in range(1001)):
            point = spline.point(t)
            distances = []
            for a, b in zip(cubic, cubic[1:]):
                dx, dy = b[0] - a[0], b[1] - a[1]
                projection = max(0, min(1, ((point.x - a[0]) * dx + (point.y - a[1]) * dy) / (dx * dx + dy * dy)))
                distances.append(math.hypot(point.x - a[0] - projection * dx, point.y - a[1] - projection * dy))
            self.assertLessEqual(min(distances), 1.1e-6)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--write-fixtures":
        write_fixtures(Path(sys.argv[2]))
    else:
        unittest.main()

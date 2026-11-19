"""Independent geometry expectations for the mesh-only report."""
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location(
    "rpa_mesh_report", Path(__file__).resolve().parents[1] / "verification/rpa_mesh_topology_report.py")
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def write_rectangles(path):
    root = ET.Element("VTKFile")
    piece = ET.SubElement(root, "Piece")
    points = ET.SubElement(piece, "Points")
    ET.SubElement(points, "DataArray").text = " ".join(str(v) for x in range(6) for r in (0, 1) for v in (x, r, 0))
    cells = ET.SubElement(piece, "Cells")
    for name, values in {
        "connectivity": [node for i in range(5) for node in (2*i, 2*i+2, 2*i+3, 2*i+1)],
        "offsets": [4, 8, 12, 16, 20], "types": [9]*5,
    }.items():
        ET.SubElement(cells, "DataArray", Name=name).text = " ".join(map(str, values))
    data = ET.SubElement(piece, "CellData")
    for name, values in {"axisymmetric_volume": [math.pi]*5, "rpa_region": [1, 2, 3, 4, 5],
                         "cell_quality": [1]*5}.items():
        ET.SubElement(data, "DataArray", Name=name).text = " ".join(map(str, values))
    ET.ElementTree(root).write(path)


class ReportTests(unittest.TestCase):
    def test_unit_cylinders_and_shared_faces(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mesh.vtu"
            write_rectangles(path)
            result, cells = report.audit(path)
        self.assertEqual(result["cells"], 5)
        self.assertEqual(result["faces"], 16)
        for cell in cells:
            self.assertAlmostEqual(cell["area"], 1)
            self.assertAlmostEqual(cell["volume"], math.pi)
            self.assertAlmostEqual(cell["h_axi"], .25)
            self.assertAlmostEqual(cell["aspect"], 1)
            self.assertAlmostEqual(cell["angular_quality"], 1)
            self.assertAlmostEqual(cell["nonorthogonality"], 0)
            self.assertAlmostEqual(cell["gradation"], 1)

    def test_inverted_cells_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mesh.vtu"
            write_rectangles(path)
            tree = ET.parse(path)
            node = tree.find('.//Cells/DataArray[@Name="connectivity"]')
            node.text = " ".join(node.text.split()[3::-1] + node.text.split()[4:])
            tree.write(path)
            with self.assertRaises(AssertionError):
                report.audit(path)

    def test_visualization_subset_preserves_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'mesh.vtu'
            write_rectangles(path)
            _,cells=report.audit(path)
            report.write_subsets(path,cells,Path(directory))
            piece=ET.parse(Path(directory)/'engine.vtu').find('.//Piece')
            self.assertEqual(piece.attrib['NumberOfCells'],'1')
            self.assertEqual(piece.attrib['NumberOfPoints'],'4')
            self.assertEqual(float(piece.find('CellData/DataArray[@Name="axisymmetric_volume"]').text),math.pi)


if __name__ == "__main__":
    unittest.main()

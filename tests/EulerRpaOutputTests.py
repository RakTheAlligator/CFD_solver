"""Regression checks for the independent accepted-state output validator."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    "rpa_validator", Path(__file__).resolve().parents[1] / "verification/validate_rpa_euler_output.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


class AcceptedHistoryTests(unittest.TestCase):
    def test_mixed_cell_offsets_and_common_edge(self):
        cells = dict(types=[9, 5, 5], offsets=[4, 7, 10],
                     connectivity=[0, 1, 4, 3, 1, 2, 4, 2, 5, 4])
        polygons = validator.validated_polygons(cells, 3, 6)
        self.assertEqual(list(map(len, polygons)), [4, 3, 3])
        edges = {}
        for cell, ids in enumerate(polygons):
            for i, node in enumerate(ids):
                edges.setdefault(tuple(sorted((node, ids[(i+1) % len(ids)]))), []).append(cell)
        self.assertEqual(edges[(1, 4)], [0, 1])
        for bad in (dict(cells, offsets=[4, 8, 10]), dict(cells, types=[9, 9, 5]),
                    dict(cells, connectivity=cells['connectivity']+[0]),
                    dict(cells, connectivity=[9]+cells['connectivity'][1:])):
            with self.assertRaises(ValueError):
                validator.validated_polygons(bad, 3, 6)

    def test_endpoint_must_not_be_an_older_accepted_time(self):
        history = [{"time": "0.1"}, {"time": "0.2"}]
        validator.require_history_endpoint(history, 0.2)
        with self.assertRaisesRegex(ValueError, "history endpoint"):
            validator.require_history_endpoint(history, 0.1)
        validator.require_history_endpoint([], 0.0)

    def test_time_accumulation_and_corruption(self):
        previous, step = 0.000244, 1.6e-9
        validator.require_time_increment(previous, previous + step, step)
        with self.assertRaises(ValueError):
            validator.require_time_increment(previous, previous + step * 2, step)


if __name__ == "__main__":
    unittest.main()

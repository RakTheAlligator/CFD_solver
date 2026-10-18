"""Regression checks for the independent accepted-state output validator."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    "rpa_validator", Path(__file__).resolve().parents[1] / "verification/validate_rpa_euler_output.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


class AcceptedHistoryTests(unittest.TestCase):
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

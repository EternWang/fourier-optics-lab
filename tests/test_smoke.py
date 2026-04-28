from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class FourierOpticsSmokeTest(unittest.TestCase):
    def test_rebuilds_processed_results(self) -> None:
        subprocess.run(
            [sys.executable, "analysis/analyze.py"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        results_path = ROOT / "data" / "processed" / "results.json"
        comparison_path = ROOT / "analysis" / "output" / "grating_method_comparison.png"

        self.assertTrue(results_path.exists())
        self.assertGreater(comparison_path.stat().st_size, 0)

        results = json.loads(results_path.read_text(encoding="utf-8"))
        for key in ["screen", "camera", "slit", "abbe"]:
            self.assertIn(key, results)
        self.assertGreater(results["screen"]["r2"], 0.99)
        self.assertGreater(results["slit"]["d_um"], results["camera"]["d_um"])


if __name__ == "__main__":
    unittest.main()

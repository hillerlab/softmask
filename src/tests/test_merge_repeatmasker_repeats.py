"""Regression tests for merge_repeatmasker_repeats.sh."""

from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


SRC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = SRC_ROOT / "bin" / "merge_repeatmasker_repeats.sh"


class MergeRepeatMaskerRepeatsTests(unittest.TestCase):
    def test_ignores_repeatmasker_diagnostic_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            diagnostic = root / "all_n.out"
            diagnostic.write_text(
                "RepeatMasker quit because the file all_n.fasta only contains "
                "ambiguous bases, if any.\n",
                encoding="utf-8",
            )
            fofn = root / "inputs.fofn"
            fofn.write_text(f"{diagnostic}\n", encoding="utf-8")

            result = subprocess.run(
                [str(SCRIPT), str(fofn)],
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.reader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][0], "SW_score")


if __name__ == "__main__":
    unittest.main()

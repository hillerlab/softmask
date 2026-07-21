"""Tests for bin/group_mask_chunks.py."""

from __future__ import annotations

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "group_mask_chunks.py"


class GroupMaskChunksCliTests(unittest.TestCase):
    """Exercise the public command-line behavior."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_masked(self, name: str, content: str) -> Path:
        """Write one input file and return its path."""
        path = self.input_dir / name
        path.write_text(content, encoding="utf-8")
        return path

    def run_script(
        self,
        *,
        input_path: Path | None = None,
        output_path: Path | None = None,
        gzip_output: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        """Run the CLI with common arguments."""
        command = [
            sys.executable,
            str(SCRIPT),
            "--input",
            str(input_path or self.input_dir),
            "--output",
            str(output_path or (self.root / "result.fasta")),
        ]
        if gzip_output:
            command.append("--gzip")
        return subprocess.run(command, text=True, capture_output=True, check=False)

    def assert_failed(
        self, result: subprocess.CompletedProcess[str], text: str
    ) -> None:
        """Assert that a CLI run failed with a useful error."""
        self.assertNotEqual(result.returncode, 0, result)
        self.assertIn(text, result.stderr)

    def test_merges_unordered_chunks_with_complex_filename_prefix(self) -> None:
        self.write_masked(
            "FA-12_bat2@chr1_3-4_something_else.masked",
            ">chr1:3-4\nGT\n",
        )
        self.write_masked(
            "FA-12_bat2@chr1_1-2_something_else.masked",
            ">chr1:1-2\nAC\n",
        )

        output = self.root / "merged.fasta"
        result = self.run_script(output_path=output)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output.read_text(encoding="utf-8"), ">chr1\nAC\nGT\n")

    def test_writes_gzip_output_and_appends_suffix(self) -> None:
        self.write_masked("sample.masked", ">sample\nACGT\n")
        output = self.root / "sample.fasta"

        result = self.run_script(output_path=output, gzip_output=True)

        gz_output = self.root / "sample.fasta.gz"
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(output.exists())
        with gzip.open(gz_output, "rt", encoding="utf-8") as handle:
            self.assertEqual(handle.read(), ">sample\nACGT\n")

    def test_allows_coordinate_free_singleton(self) -> None:
        self.write_masked("mock1r_test_assembly.masked", ">mock1r\nAC\nGT\n")
        output = self.root / "singleton.fasta"

        result = self.run_script(output_path=output)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output.read_text(encoding="utf-8"), ">mock1r\nAC\nGT\n")

    def test_rejects_missing_input_and_empty_directory(self) -> None:
        missing = self.root / "missing"
        result = self.run_script(input_path=missing)
        self.assert_failed(result, "input directory not found")

        result = self.run_script()
        self.assert_failed(result, "no .masked files found")

    def test_rejects_multiple_groups(self) -> None:
        self.write_masked("first_1-2.masked", ">first:1-2\nAC\n")
        self.write_masked("second_3-4.masked", ">second:3-4\nGT\n")

        result = self.run_script()

        self.assert_failed(result, "multiple FASTA groups")

    def test_rejects_coordinate_free_file_mixed_with_chunks(self) -> None:
        self.write_masked("sample.masked", ">sample\nAC\n")
        self.write_masked("sample_1-2.masked", ">sample:1-2\nAC\n")

        result = self.run_script()

        self.assert_failed(result, "only allowed as a singleton")

    def test_rejects_malformed_and_multi_record_fasta(self) -> None:
        self.write_masked("bad.masked", "ACGT\n")
        result = self.run_script()
        self.assert_failed(result, "does not begin with a header")

        (self.input_dir / "bad.masked").write_text(
            ">first\nAC\n>second\nGT\n",
            encoding="utf-8",
        )
        result = self.run_script()
        self.assert_failed(result, "more than one record")

    def test_rejects_ambiguous_filename_ranges(self) -> None:
        self.write_masked("sample_1-2_extra_3-4.masked", ">sample:1-2\nAC\n")

        result = self.run_script()

        self.assert_failed(result, "ambiguous filename coordinates")

    def test_rejects_header_filename_coordinate_mismatch(self) -> None:
        self.write_masked("sample_1-2.masked", ">sample:2-3\nAC\n")

        result = self.run_script()

        self.assert_failed(result, "filename/header coordinate mismatch")

    def test_rejects_non_contiguous_chunks(self) -> None:
        cases = [
            (
                [
                    ("sample_1-2.masked", ">sample:1-2\nAC\n"),
                    ("sample_4-5.masked", ">sample:4-5\nGT\n"),
                ],
                "gap in chunk coordinates",
            ),
            (
                [
                    ("sample_1-3.masked", ">sample:1-3\nACG\n"),
                    ("sample_3-4.masked", ">sample:3-4\nGT\n"),
                ],
                "overlapping chunk coordinates",
            ),
            (
                [
                    ("a_1-2.masked", ">sample:1-2\nAC\n"),
                    ("b_1-2.masked", ">sample:1-2\nAC\n"),
                ],
                "duplicate chunk coordinates",
            ),
        ]
        for files, expected in cases:
            with self.subTest(expected=expected):
                for path in self.input_dir.iterdir():
                    path.unlink()
                for name, content in files:
                    self.write_masked(name, content)
                result = self.run_script()
                self.assert_failed(result, expected)

    def test_rejects_reversed_ranges_and_incorrect_lengths(self) -> None:
        self.write_masked("sample_2-1.masked", ">sample:2-1\nAC\n")
        result = self.run_script()
        self.assert_failed(result, "reversed chunk coordinates")

        (self.input_dir / "sample_2-1.masked").unlink()
        self.write_masked("sample_1-3.masked", ">sample:1-3\nAC\n")
        result = self.run_script()
        self.assert_failed(result, "sequence length does not match coordinates")

    def test_rejects_group_that_does_not_start_at_one(self) -> None:
        self.write_masked("sample_2-3.masked", ">sample:2-3\nAC\n")

        result = self.run_script()

        self.assert_failed(result, "must start at 1")

    def test_rejects_existing_output_without_changing_it(self) -> None:
        self.write_masked("sample.masked", ">sample\nACGT\n")
        output = self.root / "result.fasta"
        output.write_text("existing\n", encoding="utf-8")

        result = self.run_script(output_path=output)

        self.assert_failed(result, "output path already exists")
        self.assertEqual(output.read_text(encoding="utf-8"), "existing\n")


if __name__ == "__main__":
    unittest.main()

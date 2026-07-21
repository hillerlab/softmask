"""Docker-backed end-to-end test for ambiguity-only FASTA chunks."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
FIXTURES = SRC_ROOT / "assets" / "test_data"
TEST_CONFIG = Path(__file__).with_name("pipeline_ambiguity.config")
RUN_INTEGRATION = os.environ.get("SOFTMASK_RUN_INTEGRATION") == "1"
DEFAULT_IMAGE = "softmask:test"


def read_fasta(path: Path) -> dict[str, str]:
    """Read a small FASTA fixture or result into an ID-to-sequence mapping."""
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise AssertionError(f"duplicate FASTA identifier: {current}")
                records[current] = []
            elif current is None:
                raise AssertionError(f"sequence before FASTA header in {path}")
            else:
                records[current].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


@unittest.skipUnless(
    RUN_INTEGRATION,
    "set SOFTMASK_RUN_INTEGRATION=1 to run the Docker-backed pipeline test",
)
class PipelineAmbiguityIntegrationTests(unittest.TestCase):
    maxDiff = None

    def run_checked(
        self,
        command: list[str],
        *,
        env: dict[str, str] | None = None,
        timeout: int = 1800,
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        self.assertEqual(
            result.returncode,
            0,
            f"command failed: {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def require_executable(self, name: str) -> None:
        self.assertIsNotNone(
            shutil.which(name),
            f"{name} is required when SOFTMASK_RUN_INTEGRATION=1",
        )

    def ensure_image(self, image: str) -> None:
        inspected = subprocess.run(
            ["docker", "image", "inspect", image],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        if inspected.returncode == 0:
            return
        self.run_checked(
            [
                "docker",
                "build",
                "--build-arg",
                "MODULE_VERSION=test",
                "-f",
                str(SRC_ROOT / "assets" / "image" / "Dockerfile"),
                "-t",
                image,
                str(SRC_ROOT),
            ]
        )

    def test_pipeline_preserves_ambiguity_and_coordinates(self) -> None:
        self.require_executable("nextflow")
        self.require_executable("docker")
        self.run_checked(["docker", "info"], timeout=60)

        image = os.environ.get("SOFTMASK_TEST_IMAGE", DEFAULT_IMAGE)
        self.ensure_image(image)

        with tempfile.TemporaryDirectory(prefix="softmask-integration-") as tempdir:
            root = Path(tempdir)
            results = root / "results"
            environment = os.environ.copy()
            environment.update(
                {
                    "NXF_ANSI_LOG": "false",
                    "NXF_CONTAINER_IMAGE": image,
                    "NXF_HOME": str(root / "nxf-home"),
                }
            )
            self.run_checked(
                [
                    "nextflow",
                    "-c",
                    str(TEST_CONFIG),
                    "run",
                    str(SRC_ROOT / "main.nf"),
                    "-profile",
                    "docker",
                    "-work-dir",
                    str(root / "work"),
                    "--genome",
                    str(FIXTURES / "ambiguity_genome.fa"),
                    "--assembly_prefix",
                    "ambiguity_test",
                    "--repeat_library",
                    str(FIXTURES / "ambiguity_repeats.fa"),
                    "--chunks",
                    "120",
                    "--output_format",
                    "fasta",
                    "--outdir",
                    str(results),
                    "--use_container",
                    "true",
                ],
                env=environment,
            )

            masked_path = results / "05_final" / "ambiguity_test.masked.fasta"
            repeats_path = results / "05_final" / "ambiguity_test.repeats.tsv"
            self.assertTrue(masked_path.is_file(), masked_path)
            self.assertTrue(repeats_path.is_file(), repeats_path)

            source = read_fasta(FIXTURES / "ambiguity_genome.fa")
            masked = read_fasta(masked_path)
            self.assertEqual(masked.keys(), source.keys())

            for name, input_sequence in source.items():
                with self.subTest(record=name):
                    output_sequence = masked[name]
                    self.assertEqual(len(output_sequence), len(input_sequence))
                    self.assertEqual(output_sequence.upper(), input_sequence.upper())
                    self.assertEqual(
                        [i for i, base in enumerate(output_sequence) if base.upper() == "N"],
                        [i for i, base in enumerate(input_sequence) if base.upper() == "N"],
                    )

            self.assertEqual(masked["all_n"], source["all_n"])
            self.assertEqual(masked["ambiguity"], source["ambiguity"])
            self.assertEqual(masked["repeat_free"], source["repeat_free"])
            self.assertTrue(
                any(base.islower() for base in masked["normal_with_gap"][:120])
            )
            self.assertTrue(
                any(
                    base.islower()
                    for base in masked["mixed"]
                    if base.upper() != "N"
                )
            )

            with repeats_path.open(encoding="utf-8", newline="") as handle:
                repeat_rows = list(csv.DictReader(handle, delimiter="\t"))

            queries = [row["query"] for row in repeat_rows]
            self.assertTrue(any(query.startswith("normal_with_gap:") for query in queries))
            self.assertIn("mixed", queries)
            forbidden_queries = {
                "normal_with_gap:121-240",
                "all_n",
                "ambiguity",
                "repeat_free",
            }
            self.assertTrue(forbidden_queries.isdisjoint(queries), queries)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Reconstruct one masked FASTA record from coordinate-labelled chunk files.

The input directory must contain one or more plain-text ``*.masked`` FASTA
files. Chunk filenames contain one underscore-delimited ``<start>-<end>``
coordinate token, and chunk FASTA headers use ``<id>:<start>-<end>``. Files are
validated, sorted by their 1-based inclusive coordinates, and their sequences
are written as one FASTA record.

A single coordinate-free ``*.masked`` file is also accepted and is emitted as
one record using its complete FASTA header.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/hillerlab/softmask"
__version__ = "0.0.1"

_LOG: logging.Logger = logging.getLogger("group_mask_chunks")

# A coordinate token must be a complete underscore-delimited filename field.
_FILENAME_COORD_RE: re.Pattern[str] = re.compile(
    r"(?:^|_)(?P<start>\d+)-(?P<end>\d+)(?=_|$)"
)
_HEADER_COORD_RE: re.Pattern[str] = re.compile(
    r"^(?P<base>.+):(?P<start>\d+)-(?P<end>\d+)$"
)


@dataclass(frozen=True)
class FastaInfo:
    """Metadata collected while validating one FASTA file."""

    header: str
    sequence_length: int


@dataclass(frozen=True)
class Chunk:
    """One validated coordinate-labelled FASTA chunk."""

    path: Path
    base: str
    start: int
    end: int
    sequence_length: int


def _filename_coordinates(path: Path) -> tuple[int, int] | None:
    """Extract the one underscore-delimited coordinate token from ``path``.

    Returns ``None`` when no coordinate token exists and raises ``ValueError``
    when more than one candidate token exists.
    """
    stem = path.name.removesuffix(".masked")
    matches = list(_FILENAME_COORD_RE.finditer(stem))
    if len(matches) > 1:
        raise ValueError(
            f"ambiguous filename coordinates in {path.name!r}; "
            "expected exactly one underscore-delimited start-end token"
        )
    if not matches:
        return None
    match = matches[0]
    return int(match.group("start")), int(match.group("end"))


def _validate_header(header: str, path: Path) -> None:
    """Validate an identifier-only FASTA header."""
    if not header:
        raise ValueError(f"empty FASTA header in {path}")
    if any(char.isspace() for char in header):
        raise ValueError(
            f"FASTA header in {path} must be an identifier without whitespace: "
            f"{header!r}"
        )


def _iter_sequence_lines(path: Path) -> Iterator[str]:
    """Yield validated, non-empty sequence lines from one-record FASTA ``path``."""
    with path.open("rt", encoding="utf-8") as handle:
        raw_header = handle.readline()
        if not raw_header:
            raise ValueError(f"empty FASTA file: {path}")

        header_line = raw_header.rstrip("\r\n")
        if not header_line.startswith(">"):
            raise ValueError(f"FASTA file does not begin with a header: {path}")
        _validate_header(header_line[1:], path)

        for line_number, raw in enumerate(handle, start=2):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                raise ValueError(
                    f"FASTA file contains more than one record at {path}:{line_number}"
                )
            if not line:
                continue
            if any(char.isspace() for char in line):
                raise ValueError(
                    f"whitespace in FASTA sequence at {path}:{line_number}"
                )
            yield line


def inspect_fasta(path: Path) -> FastaInfo:
    """Validate one plain-text, one-record FASTA and return its metadata."""
    with path.open("rt", encoding="utf-8") as handle:
        raw_header = handle.readline()
        if not raw_header:
            raise ValueError(f"empty FASTA file: {path}")
        header_line = raw_header.rstrip("\r\n")
        if not header_line.startswith(">"):
            raise ValueError(f"FASTA file does not begin with a header: {path}")
        header = header_line[1:]
        _validate_header(header, path)

    sequence_length = sum(len(line) for line in _iter_sequence_lines(path))
    if sequence_length == 0:
        raise ValueError(f"FASTA record has no sequence: {path}")
    return FastaInfo(header=header, sequence_length=sequence_length)


def _parse_chunk(path: Path, coordinates: tuple[int, int]) -> Chunk:
    """Validate one coordinate-labelled chunk and return sortable metadata."""
    filename_start, filename_end = coordinates
    if filename_start < 1:
        raise ValueError(
            f"chunk coordinates must be 1-based and start at 1 or greater: "
            f"{path.name!r}"
        )
    if filename_end < filename_start:
        raise ValueError(f"reversed chunk coordinates in {path.name!r}")

    fasta = inspect_fasta(path)
    match = _HEADER_COORD_RE.fullmatch(fasta.header)
    if match is None:
        raise ValueError(
            f"chunk FASTA header must end with :start-end in {path}: {fasta.header!r}"
        )

    header_start = int(match.group("start"))
    header_end = int(match.group("end"))
    if (header_start, header_end) != (filename_start, filename_end):
        raise ValueError(
            f"filename/header coordinate mismatch in {path}: "
            f"{filename_start}-{filename_end} != {header_start}-{header_end}"
        )

    expected_length = filename_end - filename_start + 1
    if fasta.sequence_length != expected_length:
        raise ValueError(
            f"sequence length does not match coordinates in {path}: "
            f"expected {expected_length}, found {fasta.sequence_length}"
        )

    return Chunk(
        path=path,
        base=match.group("base"),
        start=filename_start,
        end=filename_end,
        sequence_length=fasta.sequence_length,
    )


def _masked_files(input_dir: Path) -> list[Path]:
    """Return immediate-child ``*.masked`` files in deterministic order."""
    files = sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.name.endswith(".masked")
    )
    if not files:
        raise ValueError(f"no .masked files found in input directory: {input_dir}")
    return files


def prepare_input(input_dir: Path) -> tuple[str, list[Path], int]:
    """Validate and order one input group.

    Returns:
        The output FASTA header, ordered input paths, and total sequence length.
    """
    paths = _masked_files(input_dir)
    coordinates = [(path, _filename_coordinates(path)) for path in paths]

    if len(paths) == 1 and coordinates[0][1] is None:
        fasta = inspect_fasta(paths[0])
        return fasta.header, paths, fasta.sequence_length

    missing = [path.name for path, coord in coordinates if coord is None]
    if missing:
        raise ValueError(
            "coordinate-free .masked files are only allowed as a singleton: "
            + ", ".join(repr(name) for name in missing)
        )

    chunks = [
        _parse_chunk(path, coord) for path, coord in coordinates if coord is not None
    ]
    chunks.sort(key=lambda chunk: (chunk.start, chunk.end, chunk.path.name))

    bases = {chunk.base for chunk in chunks}
    if len(bases) != 1:
        raise ValueError(
            "input directory contains multiple FASTA groups: "
            + ", ".join(sorted(repr(base) for base in bases))
        )

    first = chunks[0]
    if first.start != 1:
        raise ValueError(
            f"first chunk for {first.base!r} must start at 1, found {first.start}"
        )

    previous = first
    for chunk in chunks[1:]:
        expected_start = previous.end + 1
        if chunk.start == previous.start and chunk.end == previous.end:
            raise ValueError(
                f"duplicate chunk coordinates for {chunk.base!r}: "
                f"{chunk.start}-{chunk.end}"
            )
        if chunk.start <= previous.end:
            raise ValueError(
                f"overlapping chunk coordinates for {chunk.base!r}: "
                f"{previous.start}-{previous.end} and {chunk.start}-{chunk.end}"
            )
        if chunk.start != expected_start:
            raise ValueError(
                f"gap in chunk coordinates for {chunk.base!r}: "
                f"expected {expected_start}, found {chunk.start}"
            )
        previous = chunk

    return (
        first.base,
        [chunk.path for chunk in chunks],
        sum(chunk.sequence_length for chunk in chunks),
    )


@contextmanager
def _open_output(path: Path, gzip_output: bool) -> Iterator[IO[str]]:
    """Open ``path`` exclusively for plain or gzip-compressed text output."""
    handle: IO[str] = (
        gzip.open(path, "xt", encoding="utf-8", newline="\n")
        if gzip_output
        else path.open("xt", encoding="utf-8", newline="\n")
    )
    try:
        yield handle
    finally:
        handle.close()


def _cleanup(path: Path) -> None:
    """Remove a partially written output file, ignoring cleanup failures."""
    try:
        path.unlink(missing_ok=True)
    except OSError:
        _LOG.warning("could not remove partial output %s", path)


def write_grouped_fasta(
    output_path: Path,
    header: str,
    ordered_paths: Sequence[Path],
    *,
    gzip_output: bool,
) -> None:
    """Write one reconstructed FASTA record from validated ``ordered_paths``."""
    created_output = False
    try:
        with _open_output(output_path, gzip_output) as output:
            created_output = True
            output.write(f">{header}\n")
            for path in ordered_paths:
                for line in _iter_sequence_lines(path):
                    output.write(line)
                    output.write("\n")
    except BaseException:
        if created_output:
            _cleanup(output_path)
        raise


def group_mask_chunks(
    input_dir: Path,
    output_path: Path,
    *,
    gzip_output: bool,
) -> tuple[str, int, int]:
    """Validate, order, and write one group of masked FASTA chunks.

    Returns:
        The FASTA header, number of input files, and total sequence length.
    """
    if output_path.exists():
        raise FileExistsError(f"output path already exists: {output_path}")
    if not output_path.parent.is_dir():
        raise NotADirectoryError(
            f"output parent directory does not exist: {output_path.parent}"
        )

    header, ordered_paths, sequence_length = prepare_input(input_dir)
    write_grouped_fasta(
        output_path,
        header,
        ordered_paths,
        gzip_output=gzip_output,
    )
    return header, len(ordered_paths), sequence_length


def _output_path(path: Path, gzip_output: bool) -> Path:
    """Append ``.gz`` to ``path`` when gzip output was requested."""
    if gzip_output and path.suffix != ".gz":
        return path.with_name(path.name + ".gz")
    return path


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="group_mask_chunks.py",
        description=(
            "Sort coordinate-labelled .masked FASTA chunks and reconstruct "
            "one FASTA record."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="directory containing one group of plain-text .masked FASTA files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="exact output FASTA path",
    )
    parser.add_argument(
        "-G",
        "--gzip",
        action="store_true",
        help="gzip the output and append .gz to --output when needed",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose logging")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not args.input.is_dir():
        _LOG.error("input directory not found: %s", args.input)
        return 2

    output_path = _output_path(args.output, args.gzip)
    try:
        header, n_files, sequence_length = group_mask_chunks(
            args.input,
            output_path,
            gzip_output=args.gzip,
        )
    except (OSError, UnicodeError, ValueError) as exc:
        _LOG.error("%s", exc)
        return 1

    _LOG.info(
        "grouped %d file(s) into %s (%s; %d bases)",
        n_files,
        output_path,
        header,
        sequence_length,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Rename FASTA sequence headers to clean, unique identifiers.

The script auto-detects a single *header class* from the first record and
verifies that every header conforms to it. Pass ``--enforce-rename`` to skip
detection and assign sequential ``<prefix>_<n>`` names instead. The input is
streamed exactly once: it is read (gzip is detected by magic bytes, not by
extension), the renamed FASTA is written, and a TSV mapping of
``old_header<TAB>new_id`` is emitted. Memory use is independent of genome size.

By default, dots in extracted identifiers are replaced with underscores
(e.g. ``EQ202275.1`` -> ``EQ202275_1``); pass ``--keep-dots`` to preserve them.

Header classes (numbers match the original ``lgFormatStep.pl``):

* class 2 - single token, or the first token before whitespace.
            ``>chr1 /length=...``           -> ``chr1``
* class 3 - single token with exactly one dot. Kept for fidelity, but in
            practice unreachable: class 2 already matches it.
* class 4 - old NCBI ``>gi|...|gb|ACC.1| ...`` -> ``ACC_1``
* class 5 - ``token|junk``                   -> ``token``
* class 6 - header longer than ``--max-header-length`` with no whitespace, or
            forced via ``--enforce-rename``   -> ``<prefix>_<00000001>``

Note: header validation is strict.
If any header does not conform to the detected class, the run errors out (the
original sometimes proceeded with a partial mapping). Partial outputs are
removed on any failure.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import IO

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/hillerlab/softmask"
__version__ = "0.0.1"

_LOG: logging.Logger = logging.getLogger("rename_fasta_headers")

# IUPAC nucleotide alphabet (bases + ambiguity codes + N); anything else fails.
_NON_IUPAC_RE: re.Pattern[str] = re.compile(r"[^ACGTRYSWKMBDHVN]", re.IGNORECASE)

# Header-class patterns. Explicit character classes intentionally mirror the
# original Perl regexes rather than using ``\w`` (which is Unicode-aware).
_CLASS2_WORD_REST_RE: re.Pattern[str] = re.compile(r"^>([A-Za-z0-9_.\-]+)[ \t].*$")
_CLASS2_WORD_RE: re.Pattern[str] = re.compile(r"^>([A-Za-z0-9_.\-]+)$")
_CLASS3_RE: re.Pattern[str] = re.compile(r"^>([A-Za-z0-9_\-]+\.[A-Za-z0-9_\-]+)$")
_CLASS4_RE: re.Pattern[str] = re.compile(r"^>gi\|.+\|gb\|([A-Za-z0-9._\-]+)\|.*$")
_CLASS5_RE: re.Pattern[str] = re.compile(r"^>([A-Za-z0-9_\-]+)\|[A-Za-z0-9_\-:. ]+$")
_WHITESPACE_RE: re.Pattern[str] = re.compile(r"[ \t]")

_GZIP_MAGIC: bytes = b"\x1f\x8b"


class HeaderClass(IntEnum):
    """FASTA header categories recognised by the renamer.

    Values match the class numbers used in the original ``lgFormatStep.pl``.
    """

    SIMPLE = 2  # single token, or first token before whitespace
    SINGLE_DOT = 3  # one-dot token (kept for fidelity; effectively dead)
    OLD_NCBI = 4  # >gi|...|gb|ACC| ...
    PIPE = 5  # token|junk
    RUNNING_NUMBER = 6  # <prefix>_<n> (too-long headers, or --enforce-rename)


def detect_class(
    header: str,
    max_header_length: int,
    enforce_rename: bool,
) -> HeaderClass:
    """Infer the header class from the first record.

    Args:
        header: The first header line, including the leading ``>``.
        max_header_length: Headers longer than this (and containing no
            whitespace) are treated as class 6.
        enforce_rename: If ``True``, always return :attr:`HeaderClass.RUNNING_NUMBER`.

    Returns:
        The detected :class:`HeaderClass`.

    Raises:
        ValueError: If no class matches and ``enforce_rename`` is ``False``.
    """
    if enforce_rename:
        return HeaderClass.RUNNING_NUMBER
    if len(header) > max_header_length and not _WHITESPACE_RE.search(header):
        return HeaderClass.RUNNING_NUMBER
    if _CLASS2_WORD_RE.match(header) or _CLASS2_WORD_REST_RE.match(header):
        return HeaderClass.SIMPLE
    if _CLASS3_RE.match(header):
        return HeaderClass.SINGLE_DOT
    if _CLASS4_RE.match(header):
        return HeaderClass.OLD_NCBI
    if _CLASS5_RE.match(header):
        return HeaderClass.PIPE
    raise ValueError(
        f"could not determine header class for {header!r}; "
        "fix the headers or pass --enforce-rename"
    )


@dataclass
class HeaderRenamer:
    """Stateful renamer for a single, fixed header class.

    The class is decided once (from the first header) and then applied to every
    record. Uniqueness of the resulting identifiers is enforced as records are
    seen; for class 6 the running counter guarantees uniqueness.

    Attributes:
        header_class: The class to apply to every header.
        prefix: Assembly prefix used to build class-6 names (e.g. ``HLrouAeg4``).
            Required only for :attr:`HeaderClass.RUNNING_NUMBER`.
        replace_dots: If ``True`` (default), replace ``.`` with ``_`` in extracted
            identifiers (e.g. ``EQ202275.1`` -> ``EQ202275_1``). Class-6 running
            numbers are never affected.
    """

    header_class: HeaderClass
    prefix: str | None = None
    replace_dots: bool = True
    _counter: int = field(default=0, init=False, repr=False)
    _seen: dict[str, str] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.header_class is HeaderClass.RUNNING_NUMBER and not self.prefix:
            raise ValueError(
                "class 6 renaming requires an assembly prefix (--prefix / -p)"
            )

    def rename(self, header: str) -> str:
        """Return the new identifier for one raw header line.

        Args:
            header: A header line including the leading ``>``.

        Returns:
            The new identifier, without a leading ``>``.

        Raises:
            ValueError: If the header does not conform to the detected class,
                or if it would yield an identifier already used by another header.
        """
        new = self._extract(header)
        clash = self._seen.get(new)
        if clash is not None:
            raise ValueError(
                f"duplicate identifier {new!r} produced by {header!r} and "
                f"{clash!r}; fix the headers or pass --enforce-rename"
            )
        self._seen[new] = header
        return new

    def _extract(self, header: str) -> str:
        """Compute the new identifier for ``header`` per :attr:`header_class`."""
        if self.header_class is HeaderClass.RUNNING_NUMBER:
            self._counter += 1
            assert self.prefix is not None  # guaranteed by __post_init__
            return f"{self.prefix}_{self._counter:08d}"

        match self.header_class:
            case HeaderClass.SIMPLE:
                m = _CLASS2_WORD_REST_RE.match(header) or _CLASS2_WORD_RE.match(header)
            case HeaderClass.SINGLE_DOT:
                m = _CLASS3_RE.match(header)
            case HeaderClass.OLD_NCBI:
                m = _CLASS4_RE.match(header)
            case HeaderClass.PIPE:
                m = _CLASS5_RE.match(header)
            case _:  # pragma: no cover - all classes handled above
                m = None

        if m is None:
            raise ValueError(
                f"header {header!r} does not conform to detected class "
                f"{int(self.header_class)}; fix the headers or pass --enforce-rename"
            )
        name = m.group(1)
        # By default dots become underscores (EQ202275.1 -> EQ202275_1).
        return name.replace(".", "_") if self.replace_dots else name


@contextmanager
def open_fasta_in(path: Path) -> Iterator[IO[str]]:
    """Open a possibly-gzipped FASTA for text reading.

    Compression is detected from the first two bytes (gzip magic), so the
    extension does not have to be accurate.

    Args:
        path: Path to the input FASTA (plain or gzipped).

    Yields:
        A text-mode file handle.
    """
    with path.open("rb") as probe:
        is_gzip = probe.read(2) == _GZIP_MAGIC
    handle: IO[str] = gzip.open(path, "rt") if is_gzip else path.open("rt")
    try:
        yield handle
    finally:
        handle.close()


@contextmanager
def open_fasta_out(
    path: Path,
    gzip_output: bool,
    compresslevel: int,
) -> Iterator[IO[str]]:
    """Open the output FASTA for text writing, optionally gzip-compressed.

    Args:
        path: Output path.
        gzip_output: If ``True``, write gzip-compressed data.
        compresslevel: gzip level (1-9) used when ``gzip_output`` is ``True``.

    Yields:
        A text-mode file handle.
    """
    handle: IO[str] = (
        gzip.open(path, "wt", compresslevel=compresslevel)
        if gzip_output
        else path.open("wt")
    )
    try:
        yield handle
    finally:
        handle.close()


def _cleanup(*paths: Path) -> None:
    """Remove partially written output files, ignoring missing ones."""
    for p in paths:
        try:
            p.unlink(missing_ok=True)
        except OSError:
            _LOG.warning("could not remove partial output %s", p)


def rename_fasta(
    input_path: Path,
    output_path: Path,
    dict_path: Path,
    *,
    prefix: str | None,
    enforce_rename: bool,
    replace_dots: bool,
    max_header_length: int,
    validate_iupac: bool,
    gzip_output: bool,
    compresslevel: int,
) -> int:
    """Stream ``input_path`` and write the renamed FASTA plus an old->new TSV.

    The header class is determined from the first record and applied to all
    others. The input is read once; sequence lines are passed through verbatim
    (optionally validated against the IUPAC alphabet).

    Args:
        input_path: Input FASTA (plain or gzipped).
        output_path: Destination for the renamed FASTA.
        dict_path: Destination for the ``old_header<TAB>new_id`` TSV.
        prefix: Assembly prefix for class-6 names; required only if class 6 applies.
        enforce_rename: Force class 6 (sequential ``<prefix>_<n>``) for all records.
        replace_dots: If ``True``, replace ``.`` with ``_`` in extracted identifiers.
        max_header_length: Length threshold (incl. ``>``) for auto class-6 detection.
        validate_iupac: If ``True``, error out on non-IUPAC characters in sequence.
        gzip_output: If ``True``, gzip-compress the output FASTA.
        compresslevel: gzip level used when ``gzip_output`` is ``True``.

    Returns:
        The number of renamed records.

    Raises:
        ValueError: On an empty/headerless file, a class mismatch, a duplicate
            identifier, or (when enabled) a non-IUPAC sequence character.
    """
    renamer: HeaderRenamer | None = None
    n_records = 0
    try:
        with (
            open_fasta_in(input_path) as fin,
            open_fasta_out(output_path, gzip_output, compresslevel) as fout,
            dict_path.open("wt") as fdict,
        ):
            for raw in fin:
                line = raw.rstrip("\r\n")
                if line.startswith(">"):
                    if renamer is None:
                        hclass = detect_class(line, max_header_length, enforce_rename)
                        _LOG.info("detected header class %d", int(hclass))
                        renamer = HeaderRenamer(
                            hclass, prefix, replace_dots=replace_dots
                        )
                    new = renamer.rename(line)
                    fout.write(f">{new}\n")
                    fdict.write(f"{line[1:]}\t{new}\n")  # strip leading '>'
                    n_records += 1
                else:
                    if validate_iupac and line and _NON_IUPAC_RE.search(line):
                        raise ValueError(f"non-IUPAC characters in sequence: {line!r}")
                    fout.write(line)
                    fout.write("\n")
            if renamer is None:
                raise ValueError("no FASTA headers found in input")
    except BaseException:
        _cleanup(output_path, dict_path)
        raise
    return n_records


def _strip_gz(name: str) -> str:
    """Return ``name`` without a trailing ``.gz`` (if present)."""
    return name[:-3] if name.endswith(".gz") else name


def _default_output(input_path: Path, gzip_output: bool) -> Path:
    """Derive a default output path next to the input (``<stem>.fasta[.gz]``)."""
    stem = Path(_strip_gz(input_path.name)).stem
    name = f"{stem}.fasta{'.gz' if gzip_output else ''}"
    return input_path.with_name(name)


def _default_dict(output_path: Path) -> Path:
    """Derive a default mapping path from the output (``<stem>.rename.tsv``)."""
    stem = Path(_strip_gz(output_path.name)).stem
    return output_path.with_name(f"{stem}.rename.tsv")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    p = argparse.ArgumentParser(
        prog="rename_fasta_headers.py",
        description="Rename FASTA headers to clean, unique identifiers (FASTA-only).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", type=Path, help="input FASTA (.fa/.fasta, optionally .gz)")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="output FASTA path (default: <input-stem>.fasta[.gz])",
    )
    p.add_argument(
        "-d",
        "--dict",
        dest="dict_path",
        type=Path,
        default=None,
        help="output TSV path for the old<TAB>new mapping "
        "(default: <output-stem>.rename.tsv)",
    )
    p.add_argument(
        "-p",
        "--prefix",
        default=None,
        help="assembly prefix for class-6 / --enforce-rename names (e.g. HLrouAeg4); "
        "required only when class 6 is used",
    )
    p.add_argument(
        "-G",
        "--gzip",
        action="store_true",
        help="gzip the output FASTA (.fasta.gz); default writes plain .fasta",
    )
    p.add_argument(
        "--compresslevel",
        type=int,
        default=6,
        help="gzip compression level (1-9) when --gzip is set",
    )
    p.add_argument(
        "-e",
        "--enforce-rename",
        action="store_true",
        help="force sequential <prefix>_<n> names (class 6) for every record",
    )
    p.add_argument(
        "--keep-dots",
        action="store_true",
        help="keep dots in extracted IDs (default: replace '.' with '_', "
        "e.g. EQ202275.1 -> EQ202275_1)",
    )
    p.add_argument(
        "-m",
        "--max-header-length",
        type=int,
        default=30,
        help="headers longer than this (and without whitespace) trigger class 6",
    )
    p.add_argument(
        "--validate-iupac",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="error out on non-IUPAC characters in sequence lines",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="verbose logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not args.input.is_file():
        _LOG.error("input file not found: %s", args.input)
        return 2

    output: Path = args.output or _default_output(args.input, args.gzip)
    if args.gzip and output.suffix != ".gz":
        output = output.with_name(output.name + ".gz")
    dict_path: Path = args.dict_path or _default_dict(output)

    try:
        n = rename_fasta(
            args.input,
            output,
            dict_path,
            prefix=args.prefix,
            enforce_rename=args.enforce_rename,
            replace_dots=not args.keep_dots,
            max_header_length=args.max_header_length,
            validate_iupac=args.validate_iupac,
            gzip_output=args.gzip,
            compresslevel=args.compresslevel,
        )
    except ValueError as exc:
        _LOG.error("%s", exc)
        return 1

    _LOG.info("renamed %d records -> %s", n, output)
    _LOG.info("mapping -> %s", dict_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

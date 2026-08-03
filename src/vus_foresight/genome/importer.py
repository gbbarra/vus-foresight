"""Populate a gene config's coordinates from MANE Select resources.

Run once per gene, by a curator, with the reference files in hand. The result is
written back into the gene YAML as an explicit, reviewable block -- coordinates
belong in version control where a change to them shows up in a diff, not in a
cache that silently refreshes.

The CDS boundaries are *derived and checked*, not asked for: given the declared
CDS length, there is normally exactly one open reading frame in the transcript
that starts with ATG, ends with a terminator, and contains no internal stop. If
there is not exactly one, the import fails rather than picking.

The write-back edits the derived keys in place instead of reserialising the
file. A gene config carries curation notes -- why an exon label is skipped, why
the exon order is not sorted -- and ``yaml.safe_dump`` would silently delete
every one of them on each import. The splice is verified by reparsing, so a
config whose layout defeats it fails the import rather than being corrupted.
"""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

import yaml

from .reference import GeneConfig, ReferenceUnavailable, read_fasta
from .sequence import translate_codon
from .transcript import Exon

__all__ = [
    "ImportedReference",
    "read_exon_table",
    "locate_cds",
    "import_reference",
    "splice_derived",
]


class ImportedReference(dict):
    """The block written back into the gene YAML."""


def read_exon_table(path: str | Path) -> tuple[Exon, ...]:
    """Read a TSV of ``label<TAB>start<TAB>end`` in **transcript** order.

    Transcript order, not genomic order: for a minus-strand gene the coordinates
    descend, and silently sorting them would produce a transcript that is
    internally consistent and biologically wrong.
    """
    path = Path(path)
    if not path.exists():
        raise ReferenceUnavailable(f"exon table {path} is missing")
    exons: list[Exon] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0].lower() in {"label", "exon"}:
                continue
            label, start, end = row[0], int(row[1]), int(row[2])
            exons.append(Exon(label=label, start=start, end=end))
    if not exons:
        raise ReferenceUnavailable(f"exon table {path} contains no exons")
    return tuple(exons)


def locate_cds(sequence: str, cds_length: int) -> tuple[int, int]:
    """Find the unique ORF of exactly ``cds_length`` bases. 1-based, inclusive.

    Raises when zero or more than one candidate exists -- both mean the declared
    length disagrees with the sequence, and guessing would be exactly the kind of
    plausible-looking error this project is built to avoid.
    """
    candidates: list[int] = []
    for start in range(1, len(sequence) - cds_length + 2):
        if sequence[start - 1 : start + 2] != "ATG":
            continue
        window = sequence[start - 1 : start - 1 + cds_length]
        if translate_codon(window[-3:]) != "*":
            continue
        internal = any(
            translate_codon(window[i : i + 3]) == "*" for i in range(0, cds_length - 3, 3)
        )
        if internal:
            continue
        candidates.append(start)
    if len(candidates) != 1:
        raise ReferenceUnavailable(
            f"expected exactly one open reading frame of {cds_length} nt, found "
            f"{len(candidates)} at {candidates[:5]}. Check the declared cds_length "
            "and that the FASTA is the spliced transcript, not the genomic locus."
        )
    start = candidates[0]
    return start, start + cds_length - 1


DEFAULT_INDENT = "  "

#: ``(section, key, block key)`` for every scalar the importer derives.
DERIVED_SCALARS: tuple[tuple[str, str, str], ...] = (
    ("transcript", "cds_start_tx", "cds_start_tx"),
    ("transcript", "cds_end_tx", "cds_end_tx"),
    ("sequence", "sha256", "sequence_sha256"),
    ("sequence", "expected_length", "sequence_length"),
)


def _render_scalar(value: int | str) -> str:
    """Render one scalar the way YAML will read it back unchanged.

    Exon labels are strings that often look like integers (``'1'``), so the
    quoting is decided by round-tripping rather than by type inspection.
    """
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise TypeError(f"unsupported derived scalar {value!r}")
    if isinstance(value, int):
        return str(value)
    if value != "" and yaml.safe_load(value) == value and "#" not in value:
        return value
    return "'" + value.replace("'", "''") + "'"


def _indent_of(line: str) -> int:
    return len(line) - len(line.lstrip())


def _section(lines: list[str], name: str) -> tuple[int, int] | None:
    """Bounds ``(header, end)`` of a top-level mapping key, or ``None``."""
    header = next(
        (i for i, line in enumerate(lines) if re.match(rf"^{re.escape(name)}:\s*$", line)),
        None,
    )
    if header is None:
        return None
    end = len(lines)
    for j in range(header + 1, len(lines)):
        if lines[j].strip() and not lines[j][0].isspace():
            end = j
            break
    return header, end


def _body_indent(lines: list[str], header: int, end: int) -> str:
    for i in range(header + 1, end):
        stripped = lines[i].strip()
        if stripped and not stripped.startswith("#"):
            return lines[i][: _indent_of(lines[i])]
    return DEFAULT_INDENT


def _append_point(lines: list[str], header: int, end: int) -> int:
    """The line a new key should go on: after the body, before trailing noise."""
    k = end
    while k > header + 1:
        stripped = lines[k - 1].strip()
        if stripped and not stripped.startswith("#"):
            break
        k -= 1
    return k


def _trailing_comment(rest: str) -> str:
    """The ``  # ...`` tail of a value line, so a curator's note survives.

    Skipped for quoted values, where a ``#`` cannot be told from content.
    """
    if rest.strip().startswith(("'", '"')):
        return ""
    index = rest.find(" #")
    return rest[index:] if index != -1 else ""


def _set_scalar(lines: list[str], section: str, key: str, value: int | str) -> list[str]:
    bounds = _section(lines, section)
    if bounds is None:
        return lines + [f"{section}:", f"{DEFAULT_INDENT}{key}: {_render_scalar(value)}"]
    header, end = bounds
    indent = _body_indent(lines, header, end)
    pattern = re.compile(rf"^{indent}{re.escape(key)}:(?P<rest>.*)$")
    for i in range(header + 1, end):
        match = pattern.match(lines[i])
        if match:
            lines = list(lines)
            lines[i] = (
                f"{indent}{key}: {_render_scalar(value)}"
                f"{_trailing_comment(match.group('rest'))}"
            )
            return lines
    lines = list(lines)
    lines.insert(_append_point(lines, header, end), f"{indent}{key}: {_render_scalar(value)}")
    return lines


def _block_end(lines: list[str], start: int, end: int, base: int) -> int:
    """End of a block sequence whose items sit at ``base`` indent, as PyYAML emits."""
    k = start
    while k < end:
        line = lines[k]
        if not line.strip():
            k += 1
            continue
        indent = _indent_of(line)
        if indent > base or (indent == base and line.lstrip().startswith("- ")):
            k += 1
            continue
        break
    return k


def _set_exons(lines: list[str], exons: list[dict[str, int | str]]) -> list[str]:
    bounds = _section(lines, "transcript")
    if bounds is None:
        raise ReferenceUnavailable("the gene config has no 'transcript:' section")
    header, end = bounds
    indent = _body_indent(lines, header, end)

    def rendered(item_indent: str) -> list[str]:
        out = [f"{indent}exons:"]
        for exon in exons:
            out.append(f"{item_indent}- label: {_render_scalar(exon['label'])}")
            out.append(f"{item_indent}  start: {_render_scalar(exon['start'])}")
            out.append(f"{item_indent}  end: {_render_scalar(exon['end'])}")
        return out

    for i in range(header + 1, end):
        if lines[i].startswith(f"{indent}exons:"):
            stop = _block_end(lines, i + 1, end, len(indent))
            # Keep whatever indent the file already uses for the items: PyYAML
            # writes them flush with the key, a hand-written config usually does
            # not, and neither should be rewritten as a side effect of an import.
            item_indent = next(
                (
                    lines[j][: _indent_of(lines[j])]
                    for j in range(i + 1, stop)
                    if lines[j].lstrip().startswith("- ")
                ),
                indent,
            )
            return lines[:i] + rendered(item_indent) + lines[stop:]
    at = _append_point(lines, header, end)
    return lines[:at] + rendered(indent) + lines[at:]


def _preflight(lines: list[str], text: str) -> None:
    """Refuse layouts where an edit would append a duplicate key instead of replacing.

    YAML lets the last of two identical keys win, so a splice that misses the
    existing one still parses to the right values while leaving a contradictory
    file behind. Reparsing cannot catch that; this can.
    """
    parsed = yaml.safe_load(text) or {}
    for section in {name for name, _, _ in DERIVED_SCALARS}:
        if section in parsed and _section(lines, section) is None:
            raise ReferenceUnavailable(
                f"the derived block could not be spliced into the config: '{section}' "
                "is not a plain block mapping on its own line. The file's layout is "
                "outside what the splice handles; fix it by hand rather than letting "
                "the import reserialise it."
            )
    transcript = parsed.get("transcript")
    if isinstance(transcript, dict) and "exons" in transcript:
        header, end = _section(lines, "transcript")
        indent = _body_indent(lines, header, end)
        if not any(
            lines[i].startswith(f"{indent}exons:") for i in range(header + 1, end)
        ):
            raise ReferenceUnavailable(
                "the derived block could not be spliced into the config: 'exons' is "
                "declared but not as a key on its own line under 'transcript'. Fix "
                "the layout by hand rather than letting the import reserialise it."
            )


def splice_derived(text: str, block: ImportedReference) -> str:
    """Write the derived keys into ``text``, leaving every other byte alone.

    Comments inside the exon list do not survive -- the list is replaced whole,
    because it is entirely derived. Everything else in the file, including the
    comments that explain the curated fields, is untouched.
    """
    lines = text.splitlines()
    _preflight(lines, text)
    for section, key, source in DERIVED_SCALARS:
        lines = _set_scalar(lines, section, key, block[source])
    lines = _set_exons(lines, block["exons"])
    spliced = "\n".join(lines) + "\n"

    parsed = yaml.safe_load(spliced)
    observed = {
        "cds_start_tx": parsed.get("transcript", {}).get("cds_start_tx"),
        "cds_end_tx": parsed.get("transcript", {}).get("cds_end_tx"),
        "exons": parsed.get("transcript", {}).get("exons"),
        "sequence_sha256": parsed.get("sequence", {}).get("sha256"),
        "sequence_length": parsed.get("sequence", {}).get("expected_length"),
    }
    if observed != dict(block):
        differing = sorted(k for k in observed if observed[k] != block[k])
        raise ReferenceUnavailable(
            f"the derived block could not be spliced into the config: {differing} did "
            "not read back as written. The file's layout is outside what the splice "
            "handles; fix it by hand rather than letting the import reserialise it."
        )
    return spliced


def import_reference(
    config: GeneConfig,
    *,
    sequence_path: str | Path,
    exon_table_path: str | Path,
    config_path: str | Path | None = None,
) -> ImportedReference:
    """Derive the transcript block and, optionally, write it back to the YAML."""
    sequence = read_fasta(Path(sequence_path))
    exons = read_exon_table(exon_table_path)

    tc = config.transcript
    if len(exons) != tc.exon_count:
        raise ReferenceUnavailable(
            f"{config.gene}: exon table has {len(exons)} exons, config declares "
            f"{tc.exon_count}"
        )
    if tc.exon_labels and tuple(e.label for e in exons) != tc.exon_labels:
        raise ReferenceUnavailable(
            f"{config.gene}: exon labels in the table do not match the config. "
            f"Table: {[e.label for e in exons]}"
        )
    spliced_length = sum(e.length for e in exons)
    if spliced_length != len(sequence):
        raise ReferenceUnavailable(
            f"{config.gene}: exons sum to {spliced_length} nt but the FASTA holds "
            f"{len(sequence)} nt"
        )

    cds_start, cds_end = locate_cds(sequence, tc.cds_length)

    block = ImportedReference(
        {
            "cds_start_tx": cds_start,
            "cds_end_tx": cds_end,
            "exons": [
                {"label": e.label, "start": e.start, "end": e.end} for e in exons
            ],
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "sequence_length": len(sequence),
        }
    )

    if config_path is not None:
        path = Path(config_path)
        original = path.read_text(encoding="utf-8")
        path.write_text(splice_derived(original, block), encoding="utf-8")
    return block

"""Enforce the per-module coverage floor that the project standard declares.

CLAUDE.md section 4.3 sets two floors: 80% for the project and 95% for the
modules whose errors "produce a number that looks clinically interpretable and
is wrong". coverage.py enforces the first and has no notion of the second, so a
critical module can rot to 60% while the total sails past 80% on the back of
everything else. This closes that.

The module list is configuration, not code: it lives in pyproject under
``[tool.vus-foresight]`` so that adding a critical module is a config change a
reviewer sees, rather than an edit buried in a script.

    python tools/check_module_coverage.py            # reads coverage.json
    python tools/check_module_coverage.py --json PATH
"""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_policy(pyproject: Path) -> tuple[list[str], int]:
    with pyproject.open("rb") as handle:
        config = tomllib.load(handle)
    section = config.get("tool", {}).get("vus-foresight", {})
    modules = section.get("critical_modules")
    floor = section.get("critical_coverage_floor")
    if not modules or floor is None:
        raise SystemExit(
            "pyproject declares no [tool.vus-foresight] critical_modules or "
            "critical_coverage_floor; the per-module floor cannot be enforced "
            "against a policy that does not exist."
        )
    return list(modules), int(floor)


def branch_coverage(entry: dict) -> float:
    """Branch coverage for one file, as coverage.py's own totals define it.

    Recomputed here rather than taken from ``percent_covered`` because that key
    reports line coverage when branch measurement is off, and a floor that
    silently switches metric is worse than no floor.
    """
    summary = entry["summary"]
    statements = summary["num_statements"]
    branches = summary.get("num_branches", 0)
    covered = summary["covered_lines"] + summary.get("covered_branches", 0)
    total = statements + branches
    return 100.0 * covered / total if total else 100.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=REPO_ROOT / "coverage.json")
    parser.add_argument("--pyproject", type=Path, default=REPO_ROOT / "pyproject.toml")
    args = parser.parse_args()

    if not args.json.exists():
        raise SystemExit(
            f"{args.json} not found. Produce it with:\n"
            "  pytest -m 'not validation_study' --cov --cov-report=json"
        )

    modules, floor = load_policy(args.pyproject)
    report = json.loads(args.json.read_text(encoding="utf-8"))
    files = report["files"]

    print(f"per-module floor: {floor}% branch coverage ({len(modules)} critical modules)\n")
    failures: list[tuple[str, float]] = []
    missing: list[str] = []

    for module in modules:
        key = next((k for k in files if k.replace("\\", "/").endswith(module)), None)
        if key is None:
            # A module named in the policy and absent from the report is a
            # silent hole: it would otherwise pass by never being checked.
            missing.append(module)
            print(f"  {module:36s}   ABSENT FROM REPORT")
            continue
        pct = branch_coverage(files[key])
        ok = pct >= floor
        if not ok:
            failures.append((module, pct))
        print(f"  {module:36s} {pct:6.1f}%   {'ok' if ok else 'BELOW FLOOR'}")

    if missing:
        print(f"\n{len(missing)} module(s) named in the policy are not in the report:")
        for module in missing:
            print(f"  {module}")
    if failures:
        print(f"\n{len(failures)} critical module(s) below the {floor}% floor:")
        for module, pct in sorted(failures, key=lambda item: item[1]):
            print(f"  {module:36s} {pct:6.1f}%  needs {floor - pct:.1f} more points")
    if missing or failures:
        return 1

    print(f"\nall {len(modules)} critical modules at or above {floor}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())

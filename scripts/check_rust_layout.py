#!/usr/bin/env python3
"""Check human-readable repository Rust sources, without building or expanding macros.

The limits measure reader load in authored source; compiler-expanded output is
irrelevant. Useful macros do not need expanding to satisfy this audit.
"""

from collections import Counter
from pathlib import Path
import subprocess
import sys


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    paths = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", "*.rs"],
        cwd=root,
    ).split(b"\0")
    sources = sorted({Path(path.decode("utf-8")) for path in paths if path})
    counts: Counter[Path] = Counter()
    failures = 0
    warnings = 0
    existing = 0
    for relative in sources:
        source = root / relative
        # Deleted tracked files are not part of the current implementation.
        if not source.exists():
            continue
        if source.is_symlink() or not source.is_file():
            print(f"ERROR {relative}: expected a regular Rust source file")
            failures += 1
            continue
        existing += 1
        counts[relative.parent] += 1
        data = source.read_bytes()
        lines = data.count(b"\n") + bool(data and not data.endswith(b"\n"))
        if lines > 5000:
            print(f"ERROR {relative}: {lines:,} lines (maximum 5,000)")
            failures += 1
        elif lines >= 3000:
            print(f"WARNING {relative}: {lines:,} lines (ownership review at 3,000)")
            warnings += 1
    for directory, count in sorted(counts.items()):
        if count > 20:
            print(f"ERROR {directory}/: {count} direct Rust files (maximum 20)")
            failures += 1
    print(
        f"Checked {existing} Rust files in {len(counts)} directories: "
        f"{failures} violations, {warnings} size warnings."
    )
    return int(failures != 0)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Run the desktop workflow suites against the isolated local QA browser/preview.

Requires the prerequisites documented in README.md. Each suite clears test-origin
storage. Performance measurements remain a separate, machine-specific check.
"""
from pathlib import Path
import os
import subprocess
import sys
import time


SUITES = (
    "desktop_mvp_browser", "browser", "format_browser", "biology_browser",
    "primer_source_browser", "cloning_gel_browser", "workspace_browser",
    "reads_oligos_browser", "annotation_library_browser", "physical_browser",
    "clipboard_legacy_browser", "search_tracks_browser", "alignment_browser",
    "history_graph_browser", "protein_browser", "methods_browser",
    "run_history_browser", "codon_browser", "dotplot_browser", "contig_browser",
    "methylation_browser", "checked_restriction_browser", "crispr_browser",
    "genome_browser",
)


def main():
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.pop("COLLAGENE_QA_MOBILE", None)
    for suite in SUITES:
        started = time.monotonic()
        print(f"RUN {suite}_qa.py", flush=True)
        result = subprocess.run(
            [sys.executable, "-u", str(root / "scripts" / f"{suite}_qa.py")],
            cwd=root, env=env,
        )
        print(f"RESULT {suite}: {result.returncode} ({time.monotonic()-started:.1f}s)", flush=True)
        if result.returncode:
            return result.returncode
    print(f"PASS all {len(SUITES)} desktop workflow suites", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

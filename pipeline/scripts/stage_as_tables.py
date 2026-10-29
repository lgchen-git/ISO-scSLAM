#!/usr/bin/env python3
"""Run the AS core calculation in scratch and publish only result tables."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()
    if args.scratch.exists():
        shutil.rmtree(args.scratch)
    if args.output.exists():
        shutil.rmtree(args.output)
    args.scratch.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    core_work = args.scratch.parent / (args.scratch.name + "_work")
    if core_work.exists():
        shutil.rmtree(core_work)
    command = [sys.executable, str(Path(__file__).with_name("as_core.py")), "--config", str(args.config), "--samples", str(args.samples), "--results-dir", str(args.scratch), "--work-dir", str(core_work)]
    with args.log.open("w", encoding="utf-8") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    for name in ("tables", "provenance"):
        source = args.scratch / name
        if source.is_dir():
            shutil.copytree(source, args.output / name)
    (args.output / "stage_complete.txt").write_text("AS tables complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()

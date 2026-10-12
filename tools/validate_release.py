#!/usr/bin/env python3
"""Check repository structure and syntax without executing analysis scripts."""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    failures = []
    counts = {"Python": 0, "R": 0, "Bash": 0}
    for path in sorted((ROOT / "pipelines").rglob("*")):
        if not path.is_file():
            continue
        if path.name in {".DS_Store", ".Rhistory"}:
            failures.append(f"Workstation file: {path.relative_to(ROOT)}")
        if path.suffix == ".py":
            try:
                ast.parse(path.read_text(), filename=str(path))
                counts["Python"] += 1
            except SyntaxError as exc:
                failures.append(str(exc))
        elif path.suffix == ".sh":
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            counts["Bash"] += 1
            if result.returncode:
                failures.append(result.stderr)
        elif path.suffix == ".R":
            counts["R"] += 1
    for path in [ROOT / "run.py", *sorted((ROOT / "tools").glob("*.py")), *sorted((ROOT / "tests").glob("*.py"))]:
        ast.parse(path.read_text(), filename=str(path))
    r_files = sorted((ROOT / "pipelines").rglob("*.R")) + sorted((ROOT / "tools").glob("*.R"))
    if shutil.which("Rscript"):
        result = subprocess.run(["Rscript", "-e", 'for (f in commandArgs(TRUE)) parse(file=f)', *map(str, r_files)], capture_output=True, text=True)
        if result.returncode:
            failures.append(result.stderr)
    elif r_files:
        failures.append("Rscript is required to validate R source files.")
    stages = json.loads((ROOT / "config/workflows.json").read_text())
    for stage, commands in stages.items():
        for command in commands:
            script = ROOT / command[0]
            if not script.is_file() or ROOT not in script.resolve().parents:
                failures.append(f"Invalid script in {stage}: {command[0]}")
    for name, count in counts.items():
        print(f"{name} pipeline files: {count}")
    print(f"Workflow stages: {len(stages)}")
    for message in failures:
        print(message, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

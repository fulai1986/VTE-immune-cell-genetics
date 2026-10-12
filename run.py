#!/usr/bin/env python3
"""List or run manuscript analysis stages with an explicit workspace."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

REPO = Path(__file__).resolve().parent


def load_stages():
    return json.loads((REPO / "config/workflows.json").read_text())


def main():
    stages = load_stages()
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    run = sub.add_parser("run")
    run.add_argument("stage", choices=stages)
    run.add_argument("--project", type=Path, default=os.environ.get("PROJECT_DIR"))
    run.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.command == "list":
        for name, commands in stages.items():
            print(f"{name}: {len(commands)} command(s)")
        return
    if args.project is None:
        parser.error("Provide --project or set PROJECT_DIR.")
    workspace = args.project.expanduser().resolve()
    if not workspace.is_dir():
        parser.error(f"Workspace does not exist: {workspace}")
    env = os.environ.copy()
    env["PROJECT_DIR"] = str(workspace)
    env["PYTHON"] = sys.executable
    env.setdefault("MPLBACKEND", "Agg")
    for item in stages[args.stage]:
        script = REPO / item[0]
        if not script.is_file():
            raise FileNotFoundError(script)
        interpreter = {".py": sys.executable, ".R": "Rscript", ".sh": "bash"}[script.suffix]
        command = [interpreter, str(script), *item[1:]]
        print(shlex.join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=workspace, env=env, check=True)


if __name__ == "__main__":
    main()

"""Command line entry points. Exit 0: complete; 1: requirement failure; 2: error."""

from __future__ import annotations

import argparse
import json
import sys
import uuid
from pathlib import Path

from .benchmark import benchmark, compare_benchmarks
from .campaigns import campaign, compare
from .config import document, load
from .execution import (
    evaluate_run,
    execute,
    matlab_process,
    project_root,
    reproduce,
    verify_artifacts,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="controlbench")
    commands = parser.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor")
    doctor.add_argument("--output", type=Path, default=Path("results"))
    for name in ("run", "sweep", "benchmark"):
        p = commands.add_parser(name)
        p.add_argument("config", type=Path)
        p.add_argument("--output", type=Path, default=Path("results"))
        if name == "sweep":
            p.add_argument("--resume", type=Path)
        if name == "benchmark":
            p.add_argument("--repetitions", type=int, default=5)
    p = commands.add_parser("compare-benchmarks")
    p.add_argument("baseline", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("--relative", type=float, default=0.15)
    p = commands.add_parser("analyze")
    p.add_argument("directory", type=Path)
    p = commands.add_parser("reproduce")
    p.add_argument("directory", type=Path)
    p.add_argument("--output", type=Path, default=Path("results"))
    p.add_argument("--allow-mismatch", action="store_true")
    p = commands.add_parser("compare")
    p.add_argument("baseline", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("--rules", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "doctor":
            out = args.output.resolve() / ("doctor-" + str(uuid.uuid4()))
            out.mkdir(parents=True)
            matlab_process(project_root(), out, "doctor", 120)
            print((out / "stdout.log").read_text())
            return 0
        if args.command == "run":
            out, code = execute(load(args.config), args.output)
        elif args.command == "benchmark":
            out, code = benchmark(load(args.config), args.output, args.repetitions)
        elif args.command == "compare-benchmarks":
            v = compare_benchmarks(args.baseline, args.candidate, args.relative)
            print(json.dumps(v, indent=2))
            return 2 if v["overall"] == "ERROR" else 1 if v["overall"] == "FAIL" else 0
        elif args.command == "sweep":
            out, code = campaign(document(args.config), args.output, args.resume)
        elif args.command == "reproduce":
            out, code = reproduce(args.directory, args.output, args.allow_mismatch)
        elif args.command == "analyze":
            verify_artifacts(args.directory)
            # Reanalysis is a new artifact set, preserving the original run unchanged.
            import shutil

            out = args.directory.parent / ("analysis-" + str(uuid.uuid4()))
            shutil.copytree(args.directory, out)
            v = evaluate_run(out)
            from .data import checksum, read_json, write_json

            m = read_json(out / "manifest.json")
            m["parent"] = m["run_id"]
            m["run_id"] = out.name
            from . import __version__
            from .data import checksum as file_checksum

            m["analysis_provenance"] = {
                "controlbench_version": __version__,
                "sources": {p.name: file_checksum(p) for p in Path(__file__).parent.glob("*.py")},
            }
            m["validation"] = v["overall"]
            m["checksums"] = {
                p.relative_to(out).as_posix(): checksum(p)
                for p in out.rglob("*")
                if p.is_file() and p.name != "manifest.json"
            }
            write_json(out / "manifest.json", m)
            code = 2 if v["overall"] == "ERROR" else 1 if v["overall"] == "FAIL" else 0
        else:
            v = compare(args.baseline, args.candidate, document(args.rules))
            print(json.dumps(v, indent=2))
            return 2 if v["overall"] == "ERROR" else 1 if v["overall"] == "FAIL" else 0
        print(f"Artifacts: {out}")
        if (out / "report.html").exists():
            print(f"Report: {out / 'report.html'}")
        elif (out / "benchmark.json").exists():
            print(f"Benchmark summary: {out / 'benchmark.json'}")
        return code
    except (Exception, KeyboardInterrupt) as exc:
        print(f"ControlBench error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

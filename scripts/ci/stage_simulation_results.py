"""Verify one ControlBench run and copy its immutable result bundle for CI publishing."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from controlbench.data import checksum, read_json
from controlbench.execution import verify_artifacts


def completed_runs(results_root: Path) -> list[Path]:
    return sorted(
        (path for path in results_root.iterdir() if (path / "manifest.json").is_file()),
        key=lambda path: path.stat().st_mtime_ns,
        reverse=True,
    )


def stage_run(results_root: Path, staging_root: Path) -> Path:
    candidates = completed_runs(results_root)
    if not candidates:
        raise ValueError(f"No ControlBench run found in {results_root}")

    source = candidates[0]
    verified = verify_artifacts(source)
    validation = read_json(source / "validation.json")
    if verified["status"] != "complete" or validation["overall"] != "PASS":
        raise ValueError(
            f"Refusing to publish run {source.name}: "
            f"status={verified['status']}, validation={validation['overall']}"
        )

    destination = staging_root / source.name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)

    metadata = {
        "run_id": source.name,
        "validation": validation["overall"],
        "manifest_sha256": checksum(destination / "manifest.json"),
    }
    (staging_root / "ci-metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return destination


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--staging", type=Path, required=True)
    args = parser.parse_args()
    args.staging.mkdir(parents=True, exist_ok=True)
    staged = stage_run(args.results.resolve(), args.staging.resolve())
    print(f"Validated simulation bundle: {staged}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

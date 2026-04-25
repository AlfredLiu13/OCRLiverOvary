#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

from enrichment_analysis.config import (
    INPUT_BED_DIR,
    GREAT_OUTPUT_DIR,
    GREAT_R_SCRIPT,
    BED_FILES,
)


def check_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def run_great_for_bed(label: str, bed_file: Path, genome: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_file = output_dir / f"{bed_file.stem}.metadata.txt"

    if metadata_file.exists():
        print(f"[SKIP] GREAT already completed for {label}")
        return

    cmd = [
        "Rscript",
        str(GREAT_R_SCRIPT),
        str(bed_file),
        genome,
        str(output_dir),
    ]

    print(f"[RUN] {label}")
    print(f"      BED: {bed_file}")
    print(f"      Genome: {genome}")
    print(f"      Output: {output_dir}")

    subprocess.run(cmd, check=True)


def main() -> None:
    check_file(GREAT_R_SCRIPT, "GREAT R script")

    if not INPUT_BED_DIR.exists():
        raise FileNotFoundError(f"Input BED directory not found: {INPUT_BED_DIR}")

    GREAT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for label, info in BED_FILES.items():
        bed_file = INPUT_BED_DIR / info["file"]
        genome = info["genome"]
        output_dir = GREAT_OUTPUT_DIR / label

        check_file(bed_file, f"Input BED for {label}")
        run_great_for_bed(label, bed_file, genome, output_dir)

    print("\nGREAT analysis complete.")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] GREAT failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)
    except Exception as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
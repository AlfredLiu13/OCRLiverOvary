#!/usr/bin/env python3

from pathlib import Path
import argparse
import subprocess
import sys

DEFAULT_BEDS = [
    "human_promoters.bed",
    "human_enhancers.bed",
    "human_mapped_from_mouse_promoters.bed",
    "human_mapped_from_mouse_enhancers.bed",
    "human_shared_from_mouse.bed",
]

def check_exists(path_str, label):
    path = Path(path_str)
    if not path.exists():
        print(f"[ERROR] {label} not found: {path}", file=sys.stderr)
        sys.exit(1)
    return path

def main():
    parser = argparse.ArgumentParser(
        description="Batch-run HOMER findMotifsGenome.pl on selected BED files."
    )
    parser.add_argument(
        "--bed-dir",
        default="results/classification_results/raw_results",
        help="Directory containing BED files."
    )
    parser.add_argument(
        "--outdir",
        default="motif_analysis/findmotifs_results",
        help="Base output directory for HOMER motif results."
    )
    parser.add_argument(
        "--genome",
        required=True,
        help="Genome name or path to genome FASTA (e.g. /path/to/hg38.fa)."
    )
    parser.add_argument(
        "--homer-bin",
        default="findMotifsGenome.pl",
        help="Path to HOMER findMotifsGenome.pl executable."
    )
    parser.add_argument(
        "--beds",
        nargs="+",
        default=DEFAULT_BEDS,
        help="BED filenames to process."
    )
    parser.add_argument(
        "--size",
        default="given",
        help='Motif search size, e.g. "given", 200, 500.'
    )
    parser.add_argument(
        "--len",
        default="8,10,12",
        help='Motif lengths for HOMER, e.g. "8,10,12".'
    )
    parser.add_argument(
        "--threads",
        default="4",
        help="Number of CPU threads for HOMER (-p)."
    )
    parser.add_argument(
        "--mask",
        action="store_true",
        help="Use repeat masking (-mask)."
    )
    parser.add_argument(
        "--bg",
        default=None,
        help="Optional background BED file."
    )

    args = parser.parse_args()

    bed_dir = check_exists(args.bed_dir, "BED directory")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome = args.genome
    homer_bin = args.homer_bin

    for bed_name in args.beds:
        bed_path = bed_dir / bed_name
        if not bed_path.exists():
            print(f"[WARN] BED file missing, skipping: {bed_path}", file=sys.stderr)
            continue

        prefix = bed_path.stem
        sample_outdir = outdir / prefix
        sample_outdir.mkdir(parents=True, exist_ok=True)

        cmd = [
            homer_bin,
            str(bed_path),
            genome,
            str(sample_outdir),
            "-size", str(args.size),
            "-len", str(args.len),
            "-p", str(args.threads),
        ]

        if args.mask:
            cmd.append("-mask")

        if args.bg:
            cmd.extend(["-bg", str(args.bg)])

        print(f"\n[INFO] Running findMotifsGenome.pl for: {bed_path}", file=sys.stderr)
        print(f"[INFO] Output directory: {sample_outdir}", file=sys.stderr)
        print(f"[INFO] Command: {' '.join(cmd)}", file=sys.stderr)

        try:
            subprocess.run(cmd, check=True)
            print(f"[INFO] Finished: {prefix}", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] HOMER failed for {prefix} with exit code {e.returncode}", file=sys.stderr)

    print("\n[INFO] Batch motif analysis complete.", file=sys.stderr)

if __name__ == "__main__":
    main()
    
#!/usr/bin/env python3
#This wrapper runs HOMER annotatePeaks.pl on the classified human OCR BED files and writes one annotation table per BED file to enrichment_analysis/annotated/. By default it assumes a standard HOMER installation with the hg38 genome available, but it also supports custom FASTA + GTF/GFF3 annotation when needed.


from pathlib import Path
import argparse
import subprocess
import sys

DEFAULT_BEDS = [
    "human_promoters.bed",
    "human_enhancers.bed",
    "human_shared_from_mouse.bed",
    "human_mapped_from_mouse_promoters.bed",
    "human_mapped_from_mouse_enhancers.bed",
]

def run_annotate_peaks(homer_bin, bed_path, genome, out_path, gtf=None, gff3=None):
    cmd = [homer_bin, str(bed_path), genome]

    if gtf and gff3:
        raise ValueError("Use only one of --gtf or --gff3.")

    if gtf:
        cmd.extend(["-gtf", gtf])
    elif gff3:
        cmd.extend(["-gff3", gff3])

    with open(out_path, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)

def main():
    parser = argparse.ArgumentParser(
        description="Batch wrapper for HOMER annotatePeaks.pl on classified OCR BED files."
    )
    parser.add_argument(
        "--bed-dir",
        default="classification/results/raw_results",
        help="Directory containing BED inputs."
    )
    parser.add_argument(
        "--outdir",
        default="enrichment_analysis/annotated",
        help="Directory for annotation outputs."
    )
    parser.add_argument(
        "--genome",
        default="hg38",
        help="HOMER genome name or custom genome FASTA path (default: hg38)."
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="Optional GTF file for custom annotation."
    )
    parser.add_argument(
        "--gff3",
        default=None,
        help="Optional GFF3 file for custom annotation."
    )
    parser.add_argument(
        "--beds",
        nargs="+",
        default=DEFAULT_BEDS,
        help="BED filenames to annotate."
    )
    parser.add_argument(
        "--homer-bin",
        default="annotatePeaks.pl",
        help="Path to HOMER annotatePeaks.pl executable."
    )

    args = parser.parse_args()

    bed_dir = Path(args.bed_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for bed_name in args.beds:
        bed_path = bed_dir / bed_name
        if not bed_path.exists():
            print(f"[WARN] Missing BED file, skipping: {bed_path}", file=sys.stderr)
            continue

        out_path = outdir / f"{bed_path.stem}_annotated.txt"
        print(f"[INFO] Annotating {bed_path} -> {out_path}", file=sys.stderr)

        run_annotate_peaks(
            homer_bin=args.homer_bin,
            bed_path=bed_path,
            genome=args.genome,
            out_path=out_path,
            gtf=args.gtf,
            gff3=args.gff3,
        )

    print("[INFO] annotatePeaks batch complete.", file=sys.stderr)

if __name__ == "__main__":
    main()
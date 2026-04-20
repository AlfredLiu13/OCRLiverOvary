#!/usr/bin/env python3

# main entry point for the LiverOCRAnalysis pipeline

from pathlib import Path
import argparse
import subprocess
import sys


def run_command(cmd, step_name):
    print(f"\n[INFO] starting step: {step_name}")
    print(f"[INFO] command: {' '.join(str(x) for x in cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"[INFO] finished step: {step_name}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] step failed: {step_name}")
        print(f"[ERROR] exit code: {e.returncode}")
        sys.exit(e.returncode)


def check_file_exists(path_str, label):
    path = Path(path_str)
    if not path.exists():
        print(f"[ERROR] {label} not found: {path}")
        sys.exit(1)
    return path


def ensure_dir(path_str):
    path = Path(path_str)
    path.mkdir(parents=True, exist_ok=True)
    return path


# alignment
def run_alignment(args):
    script = check_file_exists(args.script, "HALPER alignment script")
    human_peaks = check_file_exists(args.human_peaks, "Human peaks file")
    mouse_peaks = check_file_exists(args.mouse_peaks, "Mouse peaks file")
    hal_file = check_file_exists(args.hal_file, "HAL alignment file")
    output_dir = ensure_dir(args.outdir)

    human_to_mouse_dir = output_dir / "human_to_mouse"
    mouse_to_human_dir = output_dir / "mouse_to_human"
    human_to_mouse_dir.mkdir(parents=True, exist_ok=True)
    mouse_to_human_dir.mkdir(parents=True, exist_ok=True)

    human_cmd = [
        "sbatch",
        str(script),
        "-b", str(human_peaks),
        "-o", str(human_to_mouse_dir),
        "-s", "Human",
        "-t", "Mouse",
        "-c", str(hal_file),
        "-min_len", str(args.min_len),
        "-protect_dist", str(args.protect_dist),
        "-max_frac", str(args.max_frac),
    ]

    mouse_cmd = [
        "sbatch",
        str(script),
        "-b", str(mouse_peaks),
        "-o", str(mouse_to_human_dir),
        "-s", "Mouse",
        "-t", "Human",
        "-c", str(hal_file),
        "-min_len", str(args.min_len),
        "-protect_dist", str(args.protect_dist),
        "-max_frac", str(args.max_frac),
    ]

    if args.hal_liftover_path:
        human_cmd.extend(["--halPath", str(args.hal_liftover_path)])
        mouse_cmd.extend(["--halPath", str(args.hal_liftover_path)])

    if args.keep_chr_prefix:
        human_cmd.extend(["--keepChrPrefix", str(args.keep_chr_prefix)])
        mouse_cmd.extend(["--keepChrPrefix", str(args.keep_chr_prefix)])

    if args.preserve:
        preserve_arg = ",".join(args.preserve)
        human_cmd.extend(["-preserve", preserve_arg])
        mouse_cmd.extend(["-preserve", preserve_arg])

    run_command(human_cmd, "alignment_human_to_mouse")
    run_command(mouse_cmd, "alignment_mouse_to_human")


# motif
def run_motif(args):
    script = check_file_exists(args.script, "Motif analysis script")
    run_command(["bash", str(script)], "motif_analysis")


# classification
def run_classification(args):
    config = check_file_exists(args.config, "Classification config")
    script = check_file_exists(
        args.script,
        "Classification script",
    )

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        args.log_level,
    ]
    run_command(cmd, "classification")


# annotate
def run_annotate(args):
    script = check_file_exists(args.script, "annotatePeaks wrapper")
    ensure_dir(args.outdir)

    cmd = [
        sys.executable,
        str(script),
        "--bed-dir",
        str(args.bed_dir),
        "--outdir",
        str(args.outdir),
        "--genome",
        str(args.genome),
        "--homer-bin",
        str(args.homer_bin),
    ]

    if args.gtf:
        cmd.extend(["--gtf", str(args.gtf)])
    if args.gff3:
        cmd.extend(["--gff3", str(args.gff3)])
    if args.beds:
        cmd.extend(["--beds"] + args.beds)

    run_command(cmd, "annotate_peaks")


# great
def run_great(args):
    batch_script = check_file_exists(args.script, "GREAT batch script")
    check_file_exists(args.r_script, "GREAT R script")
    ensure_dir(args.outdir)

    cmd = [
        sys.executable,
        str(batch_script),
        "--bed-dir",
        str(args.bed_dir),
        "--outdir",
        str(args.outdir),
        "--species",
        str(args.species),
        "--rscript-bin",
        str(args.rscript_bin),
        "--great-r-script",
        str(args.r_script),
    ]

    if args.beds:
        cmd.extend(["--beds"] + args.beds)

    run_command(cmd, "great_analysis")


# full pipeline
def run_full(args):
    alignment_args = argparse.Namespace(
        script=args.alignment_script,
        human_peaks=args.human_peaks,
        mouse_peaks=args.mouse_peaks,
        hal_file=args.hal_file,
        outdir=args.alignment_outdir,
        hal_liftover_path=args.hal_liftover_path,
        keep_chr_prefix=args.keep_chr_prefix,
        min_len=args.min_len,
        protect_dist=args.protect_dist,
        max_frac=args.max_frac,
        preserve=args.preserve,
    )

    motif_args = argparse.Namespace(
        script=args.motif_script,
    )

    classification_args = argparse.Namespace(
        script=args.classification_script,
        config=args.config,
        log_level=args.log_level,
    )

    annotate_args = argparse.Namespace(
        script=args.annotate_script,
        bed_dir=args.bed_dir,
        outdir=args.annotate_outdir,
        genome=args.genome,
        homer_bin=args.homer_bin,
        gtf=args.gtf,
        gff3=args.gff3,
        beds=args.annotate_beds,
    )

    great_args = argparse.Namespace(
        script=args.great_script,
        r_script=args.great_r_script,
        bed_dir=args.bed_dir,
        outdir=args.great_outdir,
        species=args.great_species,
        beds=args.great_beds,
        rscript_bin=args.rscript_bin,
    )

    run_alignment(alignment_args)
    run_motif(motif_args)
    run_classification(classification_args)
    run_annotate(annotate_args)
    run_great(great_args)

    print("\n[INFO] full pipeline finished")


def main():
    parser = argparse.ArgumentParser(
        description="CLI for the LiverOCRAnalysis pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # alignment parser
    alignment_parser = subparsers.add_parser("alignment", help="run alignment")
    alignment_parser.add_argument(
        "--script",
        default="alignment/halper_map_peak_orthologs.sh",
        help="path to the HALPER slurm script",
    )
    alignment_parser.add_argument(
        "--human-peaks",
        required=True,
        help="path to human narrowPeak or BED file",
    )
    alignment_parser.add_argument(
        "--mouse-peaks",
        required=True,
        help="path to mouse narrowPeak or BED file",
    )
    alignment_parser.add_argument(
        "--hal-file",
        required=True,
        help="path to HAL alignment file",
    )
    alignment_parser.add_argument(
        "--outdir",
        default="results/alignment",
        help="output directory for alignment results",
    )
    alignment_parser.add_argument(
        "--hal-liftover-path",
        default=None,
        help="optional path to halLiftover binary",
    )
    alignment_parser.add_argument(
        "--keep-chr-prefix",
        default=None,
        help="optional chromosome prefix filter",
    )
    alignment_parser.add_argument(
        "--min-len",
        type=int,
        default=50,
        help="minimum ortholog length",
    )
    alignment_parser.add_argument(
        "--protect-dist",
        type=int,
        default=10,
        help="summit protection distance",
    )
    alignment_parser.add_argument(
        "--max-frac",
        type=float,
        default=1.5,
        help="maximum fraction of original peak length",
    )
    alignment_parser.add_argument(
        "--preserve",
        nargs="+",
        default=None,
        help="optional narrowPeak columns to preserve",
    )

    # motif parser
    motif_parser = subparsers.add_parser("motif", help="run motif analysis")
    motif_parser.add_argument(
        "--script",
        default="motif_analysis/run_homer.sh",
        help="path to motif analysis shell script",
    )

    # classification parser
    classification_parser = subparsers.add_parser("classification", help="run classification")
    classification_parser.add_argument(
        "--script",
        default="classification/classification.py",
        help="path to classification script",
    )
    classification_parser.add_argument(
        "--config",
        required=True,
        help="path to classification YAML config",
    )
    classification_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # annotate parser
    annotate_parser = subparsers.add_parser("annotate", help="run annotatePeaks")
    annotate_parser.add_argument(
        "--script",
        default="enrichment_analysis/run_annotatepeaks.py",
        help="path to annotatePeaks wrapper script",
    )
    annotate_parser.add_argument(
        "--bed-dir",
        default="results",
        help="directory containing BED files",
    )
    annotate_parser.add_argument(
        "--outdir",
        default="results/annotated",
        help="output directory for annotatePeaks results",
    )
    annotate_parser.add_argument(
        "--genome",
        default="hg38",
        help="genome name or fasta path",
    )
    annotate_parser.add_argument(
        "--homer-bin",
        default="annotatePeaks.pl",
        help="path to HOMER annotatePeaks.pl",
    )
    annotate_parser.add_argument(
        "--gtf",
        default=None,
        help="optional GTF file",
    )
    annotate_parser.add_argument(
        "--gff3",
        default=None,
        help="optional GFF3 file",
    )
    annotate_parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional BED filenames",
    )

    # great parser
    great_parser = subparsers.add_parser("great", help="run GREAT analysis")
    great_parser.add_argument(
        "--script",
        default="enrichment_analysis/run_great_batch.py",
        help="path to GREAT batch script",
    )
    great_parser.add_argument(
        "--r-script",
        default="enrichment_analysis/run_great_online.R",
        help="path to GREAT R script",
    )
    great_parser.add_argument(
        "--bed-dir",
        default="results",
        help="directory containing BED files",
    )
    great_parser.add_argument(
        "--outdir",
        default="results/great",
        help="output directory for GREAT results",
    )
    great_parser.add_argument(
        "--species",
        default="hg38",
        help="species or genome assembly for GREAT",
    )
    great_parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional BED filenames",
    )
    great_parser.add_argument(
        "--rscript-bin",
        default="Rscript",
        help="path to Rscript executable",
    )

    # full parser
    full_parser = subparsers.add_parser("full", help="run full pipeline")
    full_parser.add_argument(
        "--alignment-script",
        default="alignment/halper_map_peak_orthologs.sh",
        help="path to the HALPER slurm script",
    )
    full_parser.add_argument(
        "--human-peaks",
        required=True,
        help="path to human narrowPeak or BED file",
    )
    full_parser.add_argument(
        "--mouse-peaks",
        required=True,
        help="path to mouse narrowPeak or BED file",
    )
    full_parser.add_argument(
        "--hal-file",
        required=True,
        help="path to HAL alignment file",
    )
    full_parser.add_argument(
        "--alignment-outdir",
        default="results/alignment",
        help="output directory for alignment results",
    )
    full_parser.add_argument(
        "--hal-liftover-path",
        default=None,
        help="optional path to halLiftover binary",
    )
    full_parser.add_argument(
        "--keep-chr-prefix",
        default=None,
        help="optional chromosome prefix filter",
    )
    full_parser.add_argument(
        "--min-len",
        type=int,
        default=50,
        help="minimum ortholog length",
    )
    full_parser.add_argument(
        "--protect-dist",
        type=int,
        default=10,
        help="summit protection distance",
    )
    full_parser.add_argument(
        "--max-frac",
        type=float,
        default=1.5,
        help="maximum fraction of original peak length",
    )
    full_parser.add_argument(
        "--preserve",
        nargs="+",
        default=None,
        help="optional narrowPeak columns to preserve",
    )
    full_parser.add_argument(
        "--motif-script",
        default="motif_analysis/run_homer.sh",
        help="path to motif analysis shell script",
    )
    full_parser.add_argument(
        "--classification-script",
        default="classification/classification.py",
        help="path to classification script",
    )
    full_parser.add_argument(
        "--config",
        required=True,
        help="path to classification YAML config",
    )
    full_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )
    full_parser.add_argument(
        "--annotate-script",
        default="enrichment_analysis/run_annotatepeaks.py",
        help="path to annotatePeaks wrapper script",
    )
    full_parser.add_argument(
        "--bed-dir",
        default="results",
        help="directory containing BED files",
    )
    full_parser.add_argument(
        "--annotate-outdir",
        default="results/annotated",
        help="output directory for annotatePeaks results",
    )
    full_parser.add_argument(
        "--genome",
        default="hg38",
        help="genome name or fasta path",
    )
    full_parser.add_argument(
        "--homer-bin",
        default="annotatePeaks.pl",
        help="path to HOMER annotatePeaks.pl",
    )
    full_parser.add_argument(
        "--gtf",
        default=None,
        help="optional GTF file",
    )
    full_parser.add_argument(
        "--gff3",
        default=None,
        help="optional GFF3 file",
    )
    full_parser.add_argument(
        "--annotate-beds",
        nargs="+",
        default=None,
        help="optional BED filenames for annotatePeaks",
    )
    full_parser.add_argument(
        "--great-script",
        default="enrichment_analysis/run_great_batch.py",
        help="path to GREAT batch script",
    )
    full_parser.add_argument(
        "--great-r-script",
        default="enrichment_analysis/run_great_online.R",
        help="path to GREAT R script",
    )
    full_parser.add_argument(
        "--great-outdir",
        default="results/great",
        help="output directory for GREAT results",
    )
    full_parser.add_argument(
        "--great-species",
        default="hg38",
        help="species or genome assembly for GREAT",
    )
    full_parser.add_argument(
        "--great-beds",
        nargs="+",
        default=None,
        help="optional BED filenames for GREAT",
    )
    full_parser.add_argument(
        "--rscript-bin",
        default="Rscript",
        help="path to Rscript executable",
    )

    args = parser.parse_args()

    if args.command == "alignment":
        run_alignment(args)
    elif args.command == "motif":
        run_motif(args)
    elif args.command == "classification":
        run_classification(args)
    elif args.command == "annotate":
        run_annotate(args)
    elif args.command == "great":
        run_great(args)
    elif args.command == "full":
        run_full(args)
    else:
        print("[ERROR] unknown command")
        sys.exit(1)


if __name__ == "__main__":
    main()
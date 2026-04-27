#!/usr/bin/env python3

"""
Simple CLI for the LiverOCRAnalysis pipeline.

This file is meant to be the main entry point for running the project modules.
It does not re-implement each analysis step. Instead, it calls the existing
module scripts in a consistent order.

Example usage:
    python main.py alignment
    python main.py preprocess --config classification/sample_config.yaml
    python main.py classification --config classification/sample_config.yaml
    python main.py motif-prepare --config motif_config.yaml
    python main.py motif-run --genome /path/to/hg38.fa
    python main.py enrichment
    python main.py full --config classification/sample_config.yaml
"""

from pathlib import Path
import argparse
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def run_command(cmd, step_name, cwd=ROOT):
    """Run one shell/Python command and stop if it fails."""
    print(f"\n[INFO] Starting: {step_name}")
    print(f"[INFO] Command: {' '.join(str(x) for x in cmd)}")

    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {step_name} failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    print(f"[INFO] Finished: {step_name}")


def require_file(path, label):
    """Check that a required file exists."""
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path

    if not path.exists():
        print(f"[ERROR] {label} not found: {path}")
        sys.exit(1)

    return path


def optional_script(*possible_paths):
    """Return the first script path that exists."""
    for path in possible_paths:
        candidate = ROOT / path
        if candidate.exists():
            return candidate

    print("[ERROR] Could not find any of these scripts:")
    for path in possible_paths:
        print(f"  - {ROOT / path}")
    sys.exit(1)


def processed_config_path(config_path):
    """Return the expected config.processed.yaml path."""
    config_path = Path(config_path)
    return config_path.with_stem(config_path.stem + ".processed")


# -------------------------
# pipeline steps
# -------------------------

def run_alignment(args):
    """Submit the HALPER alignment SLURM script."""
    script = require_file(args.script, "Alignment SLURM script")

    if args.local:
        run_command(["bash", str(script)], "alignment")
    else:
        run_command(["sbatch", str(script)], "alignment")


def run_preprocess(args):
    """Run BEDTools preprocessing and create config.processed.yaml."""
    config = require_file(args.config, "Config file")

    script = optional_script(
        "classification/bedtools_preprocessing.py",
        "bedtools_preprocessing.py",
    )

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        args.log_level,
    ]

    run_command(cmd, "bedtools preprocessing")

    processed = processed_config_path(config)
    print(f"\n[INFO] Processed config expected at: {processed}")


def run_classification(args):
    """Run preprocessing, then promoter/enhancer classification."""
    config = require_file(args.config, "Config file")

    if not args.skip_preprocess:
        preprocess_args = argparse.Namespace(
            config=config,
            log_level=args.log_level,
        )
        run_preprocess(preprocess_args)
        config = processed_config_path(config)

    config = require_file(config, "Processed classification config")
    script = require_file(args.script, "Classification script")

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        args.log_level,
    ]

    run_command(cmd, "classification")


def run_classification_summary(args):
    """Run the classification plotting/summary script."""
    script = require_file(args.script, "Classification results analysis script")
    run_command([sys.executable, str(script)], "classification summary")


def run_motif_prepare(args):
    """Prepare BED/FASTA inputs for motif analysis."""
    config = require_file(args.config, "Motif config file")
    script = require_file(args.script, "Motif input preparation script")

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        args.log_level,
    ]

    run_command(cmd, "motif input preparation")


def run_motif(args):
    """Run HOMER findMotifsGenome.pl on BED files."""
    script = require_file(args.script, "Motif batch script")

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
        "--size",
        str(args.size),
        "--len",
        str(args.lengths),
        "--threads",
        str(args.threads),
    ]

    if args.mask:
        cmd.append("--mask")

    if args.bg:
        cmd.extend(["--bg", str(args.bg)])

    if args.beds:
        cmd.extend(["--beds"] + args.beds)

    run_command(cmd, "motif analysis")


def run_enrichment(args):
    """Run GREAT, summarize results, and generate plots."""
    run_script = require_file(args.run_script, "GREAT runner script")
    summary_script = require_file(args.summary_script, "GREAT summary script")
    plot_script = require_file(args.plot_script, "GREAT plotting script")

    run_command([sys.executable, str(run_script)], "GREAT enrichment")
    run_command([sys.executable, str(summary_script)], "GREAT summary")
    run_command(["Rscript", str(plot_script)], "GREAT plots")


def run_full(args):
    """Run the main downstream pipeline steps."""
    if args.include_alignment:
        alignment_args = argparse.Namespace(
            script=args.alignment_script,
            local=args.local_alignment,
        )
        run_alignment(alignment_args)

    classification_args = argparse.Namespace(
        config=args.config,
        script=args.classification_script,
        log_level=args.log_level,
        skip_preprocess=False,
    )
    run_classification(classification_args)

    if args.run_motif:
        motif_args = argparse.Namespace(
            script=args.motif_script,
            bed_dir=args.motif_bed_dir,
            outdir=args.motif_outdir,
            genome=args.genome,
            homer_bin=args.homer_bin,
            size=args.size,
            lengths=args.lengths,
            threads=args.threads,
            mask=args.mask,
            bg=args.bg,
            beds=args.beds,
        )
        run_motif(motif_args)

    enrichment_args = argparse.Namespace(
        run_script=args.enrichment_run_script,
        summary_script=args.enrichment_summary_script,
        plot_script=args.enrichment_plot_script,
    )
    run_enrichment(enrichment_args)

    print("\n[INFO] Full pipeline finished")


# -------------------------
# CLI parser
# -------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Simple CLI for the LiverOCRAnalysis pipeline"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # alignment
    alignment_parser = subparsers.add_parser(
        "alignment",
        help="submit the HALPER alignment job",
    )
    alignment_parser.add_argument(
        "--script",
        default="alignment/run_halper.sh",
        help="path to alignment SLURM script",
    )
    alignment_parser.add_argument(
        "--local",
        action="store_true",
        help="run with bash instead of sbatch",
    )

    # preprocessing
    preprocess_parser = subparsers.add_parser(
        "preprocess",
        help="clean BED files and create config.processed.yaml",
    )
    preprocess_parser.add_argument(
        "--config",
        required=True,
        help="path to YAML config file",
    )
    preprocess_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # classification
    classification_parser = subparsers.add_parser(
        "classification",
        help="run preprocessing and OCR classification",
    )
    classification_parser.add_argument(
        "--config",
        required=True,
        help="path to raw config YAML file, or processed config if using --skip-preprocess",
    )
    classification_parser.add_argument(
        "--script",
        default="classification/classification.py",
        help="path to classification script",
    )
    classification_parser.add_argument(
        "--skip-preprocess",
        action="store_true",
        help="skip preprocessing and use the config exactly as provided",
    )
    classification_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # classification summary
    class_summary_parser = subparsers.add_parser(
        "classification-summary",
        help="run classification summary plots/statistics",
    )
    class_summary_parser.add_argument(
        "--script",
        default="classification/results_analysis.py",
        help="path to classification results analysis script",
    )

    # motif input preparation
    motif_prepare_parser = subparsers.add_parser(
        "motif-prepare",
        help="prepare BED and FASTA files for motif analysis",
    )
    motif_prepare_parser.add_argument(
        "--config",
        required=True,
        help="path to motif YAML config file",
    )
    motif_prepare_parser.add_argument(
        "--script",
        default="motif_analysis/prepare_motif_inputs.py",
        help="path to motif input preparation script",
    )
    motif_prepare_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # motif analysis
    motif_parser = subparsers.add_parser(
        "motif",
        help="run HOMER motif analysis",
    )
    motif_parser.add_argument(
        "--script",
        default="motif_analysis/run_findmotifs_batch.py",
        help="path to HOMER batch script",
    )
    motif_parser.add_argument(
        "--bed-dir",
        default="results/classification_results/raw_results",
        help="directory containing BED files",
    )
    motif_parser.add_argument(
        "--outdir",
        default="results/findmotifs_results",
        help="output directory for HOMER results",
    )
    motif_parser.add_argument(
        "--genome",
        required=True,
        help="genome name or FASTA path for HOMER",
    )
    motif_parser.add_argument(
        "--homer-bin",
        default="findMotifsGenome.pl",
        help="path to HOMER findMotifsGenome.pl",
    )
    motif_parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional list of BED filenames to process",
    )
    motif_parser.add_argument(
        "--size",
        default="given",
        help='HOMER motif size, for example "given" or 200',
    )
    motif_parser.add_argument(
        "--lengths",
        default="8,10,12",
        help='HOMER motif lengths, for example "8,10,12"',
    )
    motif_parser.add_argument(
        "--threads",
        default="4",
        help="number of HOMER threads",
    )
    motif_parser.add_argument(
        "--mask",
        action="store_true",
        help="use repeat masking in HOMER",
    )
    motif_parser.add_argument(
        "--bg",
        default=None,
        help="optional background BED file",
    )

    # enrichment
    enrichment_parser = subparsers.add_parser(
        "enrichment",
        help="run GREAT enrichment, summary, and plots",
    )
    enrichment_parser.add_argument(
        "--run-script",
        default="enrichment_analysis/run_great.py",
        help="path to GREAT runner script",
    )
    enrichment_parser.add_argument(
        "--summary-script",
        default="enrichment_analysis/summarize_great.py",
        help="path to GREAT summary script",
    )
    enrichment_parser.add_argument(
        "--plot-script",
        default="enrichment_analysis/plot_great.R",
        help="path to GREAT plotting script",
    )

    # full pipeline
    full_parser = subparsers.add_parser(
        "full",
        help="run preprocessing, classification, optional motif analysis, and enrichment",
    )
    full_parser.add_argument(
        "--config",
        required=True,
        help="path to raw YAML config file",
    )
    full_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    full_parser.add_argument(
        "--include-alignment",
        action="store_true",
        help="also submit the alignment step before downstream analysis",
    )
    full_parser.add_argument(
        "--alignment-script",
        default="alignment/run_halper.sh",
        help="path to alignment SLURM script",
    )
    full_parser.add_argument(
        "--local-alignment",
        action="store_true",
        help="run alignment script with bash instead of sbatch",
    )

    full_parser.add_argument(
        "--classification-script",
        default="classification/classification.py",
        help="path to classification script",
    )

    full_parser.add_argument(
        "--run-motif",
        action="store_true",
        help="also run motif analysis after classification",
    )
    full_parser.add_argument(
        "--motif-script",
        default="motif_analysis/run_findmotifs_batch.py",
        help="path to HOMER motif batch script",
    )
    full_parser.add_argument(
        "--motif-bed-dir",
        default="results/classification_results/raw_results",
        help="directory containing motif BED files",
    )
    full_parser.add_argument(
        "--motif-outdir",
        default="results/findmotifs_results",
        help="output directory for motif results",
    )
    full_parser.add_argument(
        "--genome",
        default=None,
        help="genome name or FASTA path for motif analysis; required if --run-motif is used",
    )
    full_parser.add_argument(
        "--homer-bin",
        default="findMotifsGenome.pl",
        help="path to HOMER findMotifsGenome.pl",
    )
    full_parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional BED filenames for motif analysis",
    )
    full_parser.add_argument(
        "--size",
        default="given",
        help="HOMER motif size",
    )
    full_parser.add_argument(
        "--lengths",
        default="8,10,12",
        help="HOMER motif lengths",
    )
    full_parser.add_argument(
        "--threads",
        default="4",
        help="number of HOMER threads",
    )
    full_parser.add_argument(
        "--mask",
        action="store_true",
        help="use HOMER repeat masking",
    )
    full_parser.add_argument(
        "--bg",
        default=None,
        help="optional HOMER background BED file",
    )

    full_parser.add_argument(
        "--enrichment-run-script",
        default="enrichment_analysis/run_great.py",
        help="path to GREAT runner script",
    )
    full_parser.add_argument(
        "--enrichment-summary-script",
        default="enrichment_analysis/summarize_great.py",
        help="path to GREAT summary script",
    )
    full_parser.add_argument(
        "--enrichment-plot-script",
        default="enrichment_analysis/plot_great.R",
        help="path to GREAT plotting script",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "alignment":
        run_alignment(args)
    elif args.command == "preprocess":
        run_preprocess(args)
    elif args.command == "classification":
        run_classification(args)
    elif args.command == "classification-summary":
        run_classification_summary(args)
    elif args.command == "motif-prepare":
        run_motif_prepare(args)
    elif args.command == "motif":
        run_motif(args)
    elif args.command == "enrichment":
        run_enrichment(args)
    elif args.command == "full":
        if args.run_motif and args.genome is None:
            print("[ERROR] --genome is required when using --run-motif")
            sys.exit(1)
        run_full(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

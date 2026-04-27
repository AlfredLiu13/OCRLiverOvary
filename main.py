# main.py - command-line interface for running the LiverOCRAnalysis pipeline
import argparse
import sys

from pipeline import (
    run_alignment,
    run_preprocess,
    run_classification,
    run_classification_summary,
    run_motif_prepare,
    run_motif,
    run_enrichment,
)


def build_parser():
    # main parser
    parser = argparse.ArgumentParser(
        description="Run the LiverOCRAnalysis pipeline"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # alignment step
    alignment_parser = subparsers.add_parser(
        "alignment",
        help="run HALPER alignment",
    )
    alignment_parser.add_argument(
        "--local",
        action="store_true",
        help="run with bash instead of sbatch",
    )

    # bed preprocessing step
    preprocess_parser = subparsers.add_parser(
        "preprocess",
        help="clean BED files before classification",
    )
    preprocess_parser.add_argument(
        "--config",
        required=True,
        help="path to config YAML file",
    )
    preprocess_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # classification step
    classification_parser = subparsers.add_parser(
        "classification",
        help="run promoter/enhancer classification",
    )
    classification_parser.add_argument(
        "--config",
        required=True,
        help="path to config YAML file",
    )
    classification_parser.add_argument(
        "--skip-preprocess",
        action="store_true",
        help="use an already processed config file",
    )
    classification_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # classification summary step
    subparsers.add_parser(
        "classification-summary",
        help="run classification summary plots",
    )

    # motif input preparation step
    motif_prepare_parser = subparsers.add_parser(
        "motif-prepare",
        help="prepare BED and FASTA files for motif analysis",
    )
    motif_prepare_parser.add_argument(
        "--config",
        required=True,
        help="path to motif config YAML file",
    )
    motif_prepare_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level",
    )

    # motif analysis step
    motif_parser = subparsers.add_parser(
        "motif",
        help="run HOMER motif analysis",
    )
    motif_parser.add_argument(
        "--genome",
        required=True,
        help="genome name or genome FASTA path",
    )
    motif_parser.add_argument(
        "--bed-dir",
        default="results/classification_results/raw_results",
        help="directory with BED files",
    )
    motif_parser.add_argument(
        "--outdir",
        default="results/findmotifs_results",
        help="output directory",
    )
    motif_parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional BED filenames to run",
    )
    motif_parser.add_argument(
        "--size",
        default="given",
        help="HOMER search size",
    )
    motif_parser.add_argument(
        "--lengths",
        default="8,10,12",
        help="HOMER motif lengths",
    )
    motif_parser.add_argument(
        "--threads",
        default="4",
        help="number of threads",
    )
    motif_parser.add_argument(
        "--homer-bin",
        default="findMotifsGenome.pl",
        help="path to findMotifsGenome.pl",
    )
    motif_parser.add_argument(
        "--mask",
        action="store_true",
        help="use repeat masking",
    )
    motif_parser.add_argument(
        "--bg",
        default=None,
        help="optional background BED file",
    )

    # enrichment step
    subparsers.add_parser(
        "enrichment",
        help="run GREAT enrichment, summary, and plots",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "alignment":
        run_alignment(local=args.local)

    elif args.command == "preprocess":
        run_preprocess(
            config=args.config,
            log_level=args.log_level,
        )

    elif args.command == "classification":
        run_classification(
            config=args.config,
            log_level=args.log_level,
            skip_preprocess=args.skip_preprocess,
        )

    elif args.command == "classification-summary":
        run_classification_summary()

    elif args.command == "motif-prepare":
        run_motif_prepare(
            config=args.config,
            log_level=args.log_level,
        )

    elif args.command == "motif":
        run_motif(
            genome=args.genome,
            bed_dir=args.bed_dir,
            outdir=args.outdir,
            beds=args.beds,
            homer_bin=args.homer_bin,
            size=args.size,
            lengths=args.lengths,
            threads=args.threads,
            mask=args.mask,
            bg=args.bg,
        )

    elif args.command == "enrichment":
        run_enrichment()

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

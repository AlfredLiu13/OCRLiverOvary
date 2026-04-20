#!/usr/bin/env python3

# main entry point for the LiverOCRAnalysis pipeline

from pathlib import Path
import argparse
import subprocess
import sys


def run_command(cmd, step_name):
    """
    run a shell command and print info about the step being run
    """
    print(f"\n[INFO] Starting step: {step_name}")
    print(f"[INFO] Command: {' '.join(str(x) for x in cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"[INFO] Finished step: {step_name}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Step failed: {step_name}")
        print(f"[ERROR] Exit code: {e.returncode}")
        sys.exit(e.returncode)


def check_file_exists(path_str, label):
    """
    check if a file exists and exit with an error if it doesn't, otherwise return the Path object
    """
    path = Path(path_str)
    if not path.exists():
        print(f"[ERROR] {label} not found: {path}")
        sys.exit(1)
    return path

def run_quality_control():
    """
    Placeholder for QC.
    Right now QC was already completed manually based on your progress report.
    """
    print("\n[INFO] Quality control step is not implemented in the CLI yet.")
    print("[INFO] This step was already completed manually for the current project.")


def run_alignment(script_path):
    """
    Run the alignment module.
    """
    script = check_file_exists(script_path, "Alignment SLURM script")
    run_command(["sbatch", str(script)], "alignment")


def run_classification(config_path, log_level="INFO"):
    """
    Run the classification module using its YAML config
    """
    config = check_file_exists(config_path, "Classification config")
    script = check_file_exists("classification/classification.py", "Classification script")

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        log_level,
    ]
    run_command(cmd, "classification")


def run_motif_analysis(script_path):
    """
    Run motif analysis
    """
    script = check_file_exists(script_path, "Motif analysis script")
    run_command(["bash", str(script)], "motif_analysis")


def run_annotate_peaks(
    bed_dir,
    outdir,
    genome="hg38",
    homer_bin="annotatePeaks.pl",
    gtf=None,
    gff3=None,
    beds=None,
):
    """
    Run HOMER annotatePeaks wrapper.
    """
    script = check_file_exists(
        "enrichment_analysis/run_annotatepeaks.py",
        "annotatePeaks wrapper",
    )

    cmd = [
        sys.executable,
        str(script),
        "--bed-dir",
        str(bed_dir),
        "--outdir",
        str(outdir),
        "--genome",
        str(genome),
        "--homer-bin",
        str(homer_bin),
    ]

    if gtf:
        cmd.extend(["--gtf", str(gtf)])
    if gff3:
        cmd.extend(["--gff3", str(gff3)])
    if beds:
        cmd.extend(["--beds"] + beds)

    run_command(cmd, "enrichment_annotation")


def run_great(script_path, bed_file, species, output_prefix):
    """
    Run GREAT / rGREAT wrapper on one BED file
    """
    script = check_file_exists(script_path, "GREAT wrapper script")
    bed = check_file_exists(bed_file, "BED file for GREAT")

    cmd = [
        "bash",
        str(script),
        str(bed),
        str(species),
        str(output_prefix),
    ]
    run_command(cmd, "great_analysis")




def main():
    parser = argparse.ArgumentParser(
        description="CLI for the LiverOCRAnalysis pipeline"
    )

    parser.add_argument(
        "--step",
        required=True,
        choices=[
            "qc",
            "alignment",
            "classification",
            "motif",
            "annotate",
            "great",
            "benchmark",
            "full",
        ],
        help="Pipeline step to run",
    )

    # common optional args
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level for Python-based steps",
    )

    # alignment args
    parser.add_argument(
        "--alignment-script",
        default="alignment/run_halper.sh",
        help="Path to the SLURM script for HALPER alignment",
    )

    # classification args
    parser.add_argument(
        "--config",
        default="classification/sample_config.yaml",
        help="Path to classification YAML config",
    )

    # motif args
    parser.add_argument(
        "--motif-script",
        default="motif_analysis/run_homer.sh",
        help="Path to motif analysis shell script",
    )

    # annotatePeaks args
    parser.add_argument(
        "--bed-dir",
        default="results",
        help="Directory containing BED files for enrichment annotation",
    )
    parser.add_argument(
        "--outdir",
        default="enrichment_analysis/annotated",
        help="Output directory for enrichment annotation results",
    )
    parser.add_argument(
        "--genome",
        default="hg38",
        help="Genome name or FASTA path for annotatePeaks",
    )
    parser.add_argument(
        "--homer-bin",
        default="annotatePeaks.pl",
        help="Path to HOMER annotatePeaks.pl",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="Optional GTF file for annotatePeaks",
    )
    parser.add_argument(
        "--gff3",
        default=None,
        help="Optional GFF3 file for annotatePeaks",
    )
    parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="Optional BED filenames to annotate",
    )

    # GREAT args
    parser.add_argument(
        "--great-script",
        default="enrichment_analysis/run_great.sh",
        help="Path to GREAT shell wrapper",
    )
    parser.add_argument(
        "--great-bed",
        default=None,
        help="Input BED file for GREAT analysis",
    )
    parser.add_argument(
        "--great-species",
        default=None,
        help='Species for GREAT, e.g. "hg38" or "mm10"',
    )
    parser.add_argument(
        "--great-output-prefix",
        default=None,
        help="Output prefix for GREAT results",
    )

    args = parser.parse_args()

    if args.step == "qc":
        run_quality_control()

    elif args.step == "alignment":
        run_alignment(args.alignment_script)

    elif args.step == "classification":
        run_classification(args.config, args.log_level)

    elif args.step == "motif":
        run_motif_analysis(args.motif_script)

    elif args.step == "annotate":
        run_annotate_peaks(
            bed_dir=args.bed_dir,
            outdir=args.outdir,
            genome=args.genome,
            homer_bin=args.homer_bin,
            gtf=args.gtf,
            gff3=args.gff3,
            beds=args.beds,
        )

    elif args.step == "great":
        if not args.great_bed or not args.great_species or not args.great_output_prefix:
            print("[ERROR] GREAT step requires:")
            print("        --great-bed")
            print("        --great-species")
            print("        --great-output-prefix")
            sys.exit(1)

        run_great(
            script_path=args.great_script,
            bed_file=args.great_bed,
            species=args.great_species,
            output_prefix=args.great_output_prefix,
        )

    elif args.step == "benchmark":
        run_benchmarking()

    elif args.step == "full":
        run_quality_control()
        run_alignment(args.alignment_script)
        run_classification(args.config, args.log_level)
        run_motif_analysis(args.motif_script)
        run_annotate_peaks(
            bed_dir=args.bed_dir,
            outdir=args.outdir,
            genome=args.genome,
            homer_bin=args.homer_bin,
            gtf=args.gtf,
            gff3=args.gff3,
            beds=args.beds,
        )
        print("\n[INFO] Full pipeline wrapper finished.")
        print("[INFO] GREAT is not run automatically in --step full because it needs a specific BED input each time.")

    else:
        print("[ERROR] Unknown step selected.")
        sys.exit(1)


if __name__ == "__main__":
    main()
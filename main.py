#!/usr/bin/env python3

# main entry point for the LiverOCRAnalysis pipeline

from pathlib import Path
import argparse
import subprocess
import sys


def run_command(cmd, step_name):
    """
    run a command and stop if it fails
    """
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
    """
    check that a file exists
    """
    path = Path(path_str)
    if not path.exists():
        print(f"[ERROR] {label} not found: {path}")
        sys.exit(1)
    return path


def ensure_dir(path_str):
    """
    make sure an output directory exists
    """
    path = Path(path_str)
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_alignment(
    halper_script,
    human_peaks,
    mouse_peaks,
    hal_file,
    output_dir,
    hal_liftover_path=None,
    keep_chr_prefix=None,
    min_len=50,
    protect_dist=10,
    max_frac=1.5,
    preserve=None,
):
    """
    submit two HALPER alignment jobs:
    1) human -> mouse
    2) mouse -> human
    """
    script = check_file_exists(halper_script, "HALPER alignment script")
    human_peaks = check_file_exists(human_peaks, "Human peaks file")
    mouse_peaks = check_file_exists(mouse_peaks, "Mouse peaks file")
    hal_file = check_file_exists(hal_file, "HAL alignment file")
    output_dir = ensure_dir(output_dir)

    human_to_mouse_dir = output_dir / "human_to_mouse"
    mouse_to_human_dir = output_dir / "mouse_to_human"
    human_to_mouse_dir.mkdir(parents=True, exist_ok=True)
    mouse_to_human_dir.mkdir(parents=True, exist_ok=True)

    # first mapping: human to mouse
    human_cmd = [
        "sbatch",
        str(script),
        "-b", str(human_peaks),
        "-o", str(human_to_mouse_dir),
        "-s", "Human",
        "-t", "Mouse",
        "-c", str(hal_file),
        "-min_len", str(min_len),
        "-protect_dist", str(protect_dist),
        "-max_frac", str(max_frac),
    ]

    # second mapping: mouse to human
    mouse_cmd = [
        "sbatch",
        str(script),
        "-b", str(mouse_peaks),
        "-o", str(mouse_to_human_dir),
        "-s", "Mouse",
        "-t", "Human",
        "-c", str(hal_file),
        "-min_len", str(min_len),
        "-protect_dist", str(protect_dist),
        "-max_frac", str(max_frac),
    ]

    # add optional arguments if the user passes them
    if hal_liftover_path:
        human_cmd.extend(["--halPath", str(hal_liftover_path)])
        mouse_cmd.extend(["--halPath", str(hal_liftover_path)])

    if keep_chr_prefix:
        human_cmd.extend(["--keepChrPrefix", str(keep_chr_prefix)])
        mouse_cmd.extend(["--keepChrPrefix", str(keep_chr_prefix)])

    if preserve:
        preserve_arg = ",".join(preserve)
        human_cmd.extend(["-preserve", preserve_arg])
        mouse_cmd.extend(["-preserve", preserve_arg])

    run_command(human_cmd, "alignment_human_to_mouse")
    run_command(mouse_cmd, "alignment_mouse_to_human")


def run_motif_analysis(script_path):
    """
    run motif analysis script
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
    run HOMER annotation wrapper
    """
    script = check_file_exists(
        "enrichment_analysis/run_annotatepeaks.py",
        "annotatePeaks wrapper",
    )

    ensure_dir(outdir)

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


def run_classification(config_path, log_level="INFO"):
    """
    run classification step
    """
    config = check_file_exists(config_path, "Classification config")
    script = check_file_exists(
        "classification/classification.py",
        "Classification script",
    )

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        log_level,
    ]
    run_command(cmd, "classification")


def main():
    parser = argparse.ArgumentParser(
        description="CLI for the LiverOCRAnalysis pipeline"
    )

    parser.add_argument(
        "--step",
        required=True,
        choices=[
            "alignment",
            "motif",
            "annotate",
            "classification",
            "full",
        ],
        help="pipeline step to run",
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="logging level for python-based steps",
    )

    # alignment arguments
    parser.add_argument(
        "--alignment-script",
        default="alignment/halper_map_peak_orthologs.sh",
        help="path to the HALPER slurm script",
    )
    parser.add_argument(
        "--human-peaks",
        help="path to human narrowPeak or BED file",
    )
    parser.add_argument(
        "--mouse-peaks",
        help="path to mouse narrowPeak or BED file",
    )
    parser.add_argument(
        "--hal-file",
        help="path to HAL alignment file",
    )
    parser.add_argument(
        "--alignment-outdir",
        default="results/alignment",
        help="output directory for alignment results",
    )
    parser.add_argument(
        "--hal-liftover-path",
        default=None,
        help="optional path to halLiftover binary",
    )
    parser.add_argument(
        "--keep-chr-prefix",
        default=None,
        help="optional chromosome prefix filter",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=50,
        help="minimum ortholog length for HALPER",
    )
    parser.add_argument(
        "--protect-dist",
        type=int,
        default=10,
        help="summit protection distance for HALPER",
    )
    parser.add_argument(
        "--max-frac",
        type=float,
        default=1.5,
        help="maximum fraction of original peak length for HALPER",
    )
    parser.add_argument(
        "--preserve",
        nargs="+",
        default=None,
        help="optional narrowPeak columns to preserve, e.g. signal pValue qValue",
    )

    # motif arguments
    parser.add_argument(
        "--motif-script",
        default="motif_analysis/run_homer.sh",
        help="path to motif analysis shell script",
    )

    # enrichment arguments
    parser.add_argument(
        "--bed-dir",
        default="results",
        help="directory containing BED files for enrichment annotation",
    )
    parser.add_argument(
        "--outdir",
        default="enrichment_analysis/annotated",
        help="output directory for enrichment annotation results",
    )
    parser.add_argument(
        "--genome",
        default="hg38",
        help="genome name or fasta path for annotatePeaks",
    )
    parser.add_argument(
        "--homer-bin",
        default="annotatePeaks.pl",
        help="path to HOMER annotatePeaks.pl",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="optional GTF file for annotation",
    )
    parser.add_argument(
        "--gff3",
        default=None,
        help="optional GFF3 file for annotation",
    )
    parser.add_argument(
        "--beds",
        nargs="+",
        default=None,
        help="optional BED filenames to annotate",
    )

    # classification arguments
    parser.add_argument(
        "--config",
        default="classification/sample_config.yaml",
        help="path to classification YAML config",
    )

    args = parser.parse_args()

    if args.step == "alignment":
        if not args.human_peaks or not args.mouse_peaks or not args.hal_file:
            print("[ERROR] alignment step requires:")
            print("        --human-peaks")
            print("        --mouse-peaks")
            print("        --hal-file")
            sys.exit(1)

        run_alignment(
            halper_script=args.alignment_script,
            human_peaks=args.human_peaks,
            mouse_peaks=args.mouse_peaks,
            hal_file=args.hal_file,
            output_dir=args.alignment_outdir,
            hal_liftover_path=args.hal_liftover_path,
            keep_chr_prefix=args.keep_chr_prefix,
            min_len=args.min_len,
            protect_dist=args.protect_dist,
            max_frac=args.max_frac,
            preserve=args.preserve,
        )

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

    elif args.step == "classification":
        run_classification(args.config, args.log_level)

    elif args.step == "full":
        if not args.human_peaks or not args.mouse_peaks or not args.hal_file:
            print("[ERROR] full pipeline requires alignment inputs:")
            print("        --human-peaks")
            print("        --mouse-peaks")
            print("        --hal-file")
            sys.exit(1)

        # run steps in the order you confirmed
        run_alignment(
            halper_script=args.alignment_script,
            human_peaks=args.human_peaks,
            mouse_peaks=args.mouse_peaks,
            hal_file=args.hal_file,
            output_dir=args.alignment_outdir,
            hal_liftover_path=args.hal_liftover_path,
            keep_chr_prefix=args.keep_chr_prefix,
            min_len=args.min_len,
            protect_dist=args.protect_dist,
            max_frac=args.max_frac,
            preserve=args.preserve,
        )
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
        run_classification(args.config, args.log_level)

        print("\n[INFO] full pipeline finished")

    else:
        print("[ERROR] unknown step selected")
        sys.exit(1)


if __name__ == "__main__":
    main()
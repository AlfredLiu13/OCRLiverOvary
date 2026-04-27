# pipeline.py - main pipeline functions to run each step of the analysis

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

# paths to scripts in the repo
ALIGNMENT_SCRIPT = ROOT / "alignment" / "run_halper.sh"
BEDTOOLS_PREPROCESS_SCRIPT = ROOT / "classification" / "bedtools_preprocessing.py"
BEDTOOLS_PREPROCESS_SCRIPT_ROOT = ROOT / "bedtools_preprocessing.py"
CLASSIFICATION_SCRIPT = ROOT / "classification" / "classification.py"
CLASSIFICATION_SUMMARY_SCRIPT = ROOT / "classification" / "results_analysis.py"
MOTIF_PREPARE_SCRIPT = ROOT / "motif_analysis" / "prepare_motif_inputs.py"
MOTIF_SCRIPT = ROOT / "motif_analysis" / "run_findmotifs_batch.py"
GREAT_RUN_SCRIPT = ROOT / "enrichment_analysis" / "run_great.py"
GREAT_SUMMARY_SCRIPT = ROOT / "enrichment_analysis" / "summarize_great.py"
GREAT_PLOT_SCRIPT = ROOT / "enrichment_analysis" / "plot_great.R"


def run_command(cmd, step_name):
    # run a command and stop if it fails
    print(f"\n[INFO] {step_name}")

    try:
        subprocess.run(cmd, cwd=ROOT, check=True)
    except subprocess.CalledProcessError:
        print(f"[ERROR] failed at: {step_name}")
        sys.exit(1)


def require_file(path, label):
    # make sure file exists before running anything
    path = Path(path)

    if not path.is_absolute():
        path = ROOT / path

    if not path.exists():
        print(f"[ERROR] {label} not found: {path}")
        sys.exit(1)

    return path


def get_preprocess_script():
    # preprocessing script might be in two places
    if BEDTOOLS_PREPROCESS_SCRIPT.exists():
        return BEDTOOLS_PREPROCESS_SCRIPT

    if BEDTOOLS_PREPROCESS_SCRIPT_ROOT.exists():
        return BEDTOOLS_PREPROCESS_SCRIPT_ROOT

    print("[ERROR] bedtools_preprocessing.py not found")
    sys.exit(1)


def get_processed_config_path(config_path):
    # turn config.yaml -> config.processed.yaml
    config_path = Path(config_path)
    return config_path.with_stem(config_path.stem + ".processed")


# -------------------------
# pipeline steps
# -------------------------

def run_alignment(local=False):
    # run halper alignment (slurm or local)
    script = require_file(ALIGNMENT_SCRIPT, "alignment script")

    if local:
        run_command(["bash", str(script)], "alignment")
    else:
        run_command(["sbatch", str(script)], "alignment")


def run_preprocess(config, log_level="INFO"):
    # clean bed files and create processed config
    config = require_file(config, "config file")
    script = get_preprocess_script()

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        log_level,
    ]

    run_command(cmd, "preprocess")

    return get_processed_config_path(config)


def run_classification(config, log_level="INFO", skip_preprocess=False):
    # run preprocessing (unless skipped) then classification
    config = require_file(config, "config file")

    if not skip_preprocess:
        config = run_preprocess(config, log_level)

    config = require_file(config, "processed config")
    script = require_file(CLASSIFICATION_SCRIPT, "classification script")

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        log_level,
    ]

    run_command(cmd, "classification")


def run_classification_summary():
    # run plots / stats for classification
    script = require_file(CLASSIFICATION_SUMMARY_SCRIPT, "summary script")
    run_command([sys.executable, str(script)], "classification summary")


def run_motif_prepare(config, log_level="INFO"):
    # prepare fasta + bed inputs for motif step
    config = require_file(config, "motif config")
    script = require_file(MOTIF_PREPARE_SCRIPT, "motif prepare script")

    cmd = [
        sys.executable,
        str(script),
        "--config",
        str(config),
        "--log-level",
        log_level,
    ]

    run_command(cmd, "motif prepare")


def run_motif(
    genome,
    bed_dir="results/classification_results/raw_results",
    outdir="results/findmotifs_results",
    beds=None,
    homer_bin="findMotifsGenome.pl",
    size="given",
    lengths="8,10,12",
    threads="4",
    mask=False,
    bg=None,
):
    # run homer motif discovery
    script = require_file(MOTIF_SCRIPT, "motif script")

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
        "--size",
        str(size),
        "--len",
        str(lengths),
        "--threads",
        str(threads),
    ]

    if beds:
        cmd.extend(["--beds"] + beds)

    if mask:
        cmd.append("--mask")

    if bg:
        cmd.extend(["--bg", str(bg)])

    run_command(cmd, "motif")


def run_enrichment():
    # run great + summarize + plots
    run_script = require_file(GREAT_RUN_SCRIPT, "great run")
    summary_script = require_file(GREAT_SUMMARY_SCRIPT, "great summary")
    plot_script = require_file(GREAT_PLOT_SCRIPT, "great plot")

    run_command([sys.executable, str(run_script)], "great run")
    run_command([sys.executable, str(summary_script)], "great summary")
    run_command(["Rscript", str(plot_script)], "great plot")

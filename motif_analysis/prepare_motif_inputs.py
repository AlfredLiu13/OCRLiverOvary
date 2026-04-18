### takes bed files from regulatory_comparison/ or results/
"""
Prepare motif analysis inputs from OCR BED files.

Purpose:
    Take already cleaned/sorted BED files (or raw BED/.gz files),
    optionally standardize them to BED3, and extract FASTA sequences
    for downstream motif analysis tools such as HOMER or MEME.

Tasks:
    - Load YAML config
    - Resolve raw or cleaned BED inputs
    - Unzip .gz files if needed (preserving originals)
    - Extract BED3 if requested
    - Sort BED files
    - Extract FASTA sequences using bedtools getfasta
    - Write processed config pointing to motif-ready BED/FASTA files

Usage:
    python prepare_motif_inputs.py --config config.yaml

Expected config keys:
    motif_input_regions:
      open_human: path/to/file.bed
      open_mouse: path/to/file.bed
      shared_open: path/to/file.bed

    genome_fastas:
      human: path/to/human.fa
      mouse: path/to/mouse.fa

Optional config keys:
    motif_output_dir: results/motif_analysis
    motif_extract_bed3: true

Output:
    - config.motif.processed.yaml
    - motif-ready BED files
    - FASTA files for motif analysis
"""

import argparse
import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


# logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def read_yaml_config(config_path: Path) -> Dict[str, Any]:
    """Read YAML config file."""
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded config: {config_path}")
        return config
    except yaml.YAMLError as e:
        logger.error(f"Failed to parse YAML config: {e}")
        raise


def resolve_file_path(file_path: Path, key_name: str) -> Path:
    """Resolve file path, allowing .gz fallback."""
    file_path = Path(file_path)

    if file_path.exists():
        return file_path

    if file_path.suffix != ".gz":
        gz_path = file_path.with_suffix(file_path.suffix + ".gz")
        if gz_path.exists():
            logger.debug(f"{key_name}: using gz version {gz_path}")
            return gz_path

    raise FileNotFoundError(f"{key_name} not found: {file_path} (or .gz variant)")


def gunzip_keep(src: Path, dst: Path) -> None:
    """Unzip .gz file without deleting original."""
    if not src.exists():
        raise FileNotFoundError(f"Source file not found: {src}")

    try:
        with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logger.debug(f"Unzipped {src} -> {dst}")
    except Exception as e:
        logger.error(f"Failed to unzip {src}: {e}")
        raise


def ensure_unzipped(path: Path) -> Path:
    """Return unzipped path if input is .gz."""
    if path.suffix != ".gz":
        return path

    unzipped = path.with_suffix("")
    if unzipped.exists():
        logger.debug(f"Already unzipped: {unzipped}")
        return unzipped

    gunzip_keep(path, unzipped)
    return unzipped


def extract_bed3(input_path: Path, output_path: Path) -> int:
    """Extract first 3 BED columns."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    try:
        with open(output_path, "w") as out:
            subprocess.run(
                ["cut", "-f1-3", str(input_path)],
                stdout=out,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )

        with open(output_path) as f:
            line_count = sum(1 for _ in f)

        logger.debug(f"Extracted BED3: {input_path} -> {output_path} ({line_count} lines)")
        return line_count

    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract BED3 from {input_path}: {e.stderr}")
        raise


def sort_bed(input_path: Path, output_path: Path) -> int:
    """Sort BED file for bedtools."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    try:
        with open(output_path, "w") as out:
            subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(input_path)],
                stdout=out,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )

        with open(output_path) as f:
            line_count = sum(1 for _ in f)

        logger.debug(f"Sorted BED: {input_path} -> {output_path} ({line_count} lines)")
        return line_count

    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to sort BED file {input_path}: {e.stderr}")
        raise


def check_bedtools_installed() -> None:
    """Check if bedtools is available."""
    try:
        subprocess.run(
            ["bedtools", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
    except Exception:
        raise RuntimeError(
            "bedtools is not installed or not available in PATH. "
            "Please load/install bedtools before running this script."
        )


def extract_fasta_from_bed(bed_path: Path, genome_fasta: Path, output_fasta: Path) -> None:
    """Extract sequences from BED regions using bedtools getfasta."""
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fasta}")

    try:
        subprocess.run(
            [
                "bedtools", "getfasta",
                "-fi", str(genome_fasta),
                "-bed", str(bed_path),
                "-fo", str(output_fasta)
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        logger.info(f"Extracted FASTA: {output_fasta}")
    except subprocess.CalledProcessError as e:
        logger.error(f"bedtools getfasta failed for {bed_path}: {e.stderr}")
        raise


def prepare_single_region_file(
    region_name: str,
    region_path: Path,
    genome_fasta: Path,
    output_dir: Path,
    extract_to_bed3: bool = True,
) -> Dict[str, str]:
    """
    Prepare one BED file for motif analysis and extract FASTA.
    """
    logger.info(f"Preparing motif input for: {region_name}")

    region_path = resolve_file_path(region_path, region_name)
    region_path = ensure_unzipped(region_path)

    region_output_dir = output_dir / region_name
    region_output_dir.mkdir(parents=True, exist_ok=True)

    working_bed = region_output_dir / region_path.name

    if extract_to_bed3:
        logger.info(f"  Extracting BED3 for {region_name}")
        extract_bed3(region_path, working_bed)
    else:
        shutil.copy2(region_path, working_bed)
        logger.info(f"  Copied BED without BED3 extraction: {working_bed}")

    sorted_bed = region_output_dir / f"{working_bed.stem}.sorted{working_bed.suffix}"
    sort_bed(working_bed, sorted_bed)

    final_bed = region_output_dir / f"{region_name}.bed"
    shutil.move(str(sorted_bed), str(final_bed))
    logger.info(f"  Final BED ready: {final_bed}")

    fasta_path = region_output_dir / f"{region_name}.fa"
    extract_fasta_from_bed(final_bed, genome_fasta, fasta_path)

    return {
        "bed": str(final_bed),
        "fasta": str(fasta_path),
        "genome_fasta": str(genome_fasta),
    }


def prepare_motif_inputs(config_path: Path) -> Path:
    """
    Main function to prepare motif BED and FASTA inputs.
    """
    config_path = Path(config_path)
    logger.info(f"Starting motif input preparation: {config_path}")

    raw_config = read_yaml_config(config_path)

    motif_regions = raw_config.get("motif_input_regions")
    genome_fastas = raw_config.get("genome_fastas")

    if not motif_regions:
        raise ValueError(
            "Config is missing 'motif_input_regions'. "
            "Expected something like: {'open_human': 'path/to/file.bed'}"
        )

    if not genome_fastas:
        raise ValueError(
            "Config is missing 'genome_fastas'. "
            "Expected something like: {'human': 'path/to/human.fa'}"
        )

    output_dir = Path(raw_config.get("motif_output_dir", "results/motif_analysis"))
    output_dir.mkdir(parents=True, exist_ok=True)

    extract_to_bed3 = raw_config.get("motif_extract_bed3", True)

    check_bedtools_installed()

    processed_config = raw_config.copy()
    processed_config["motif_prepared_inputs"] = {}

    for region_name, region_info in motif_regions.items():
        if isinstance(region_info, dict):
            bed_path = region_info.get("bed")
            genome_key = region_info.get("genome")
        else:
            raise ValueError(
                f"Each entry in motif_input_regions must be a dictionary. "
                f"Problem found in '{region_name}'."
            )

        if not bed_path:
            raise ValueError(f"Missing 'bed' path for motif_input_regions['{region_name}']")
        if not genome_key:
            raise ValueError(f"Missing 'genome' key for motif_input_regions['{region_name}']")
        if genome_key not in genome_fastas:
            raise ValueError(
                f"Genome key '{genome_key}' for region '{region_name}' "
                f"not found in genome_fastas"
            )

        genome_fasta = resolve_file_path(Path(genome_fastas[genome_key]), f"genome_fastas.{genome_key}")

        prepared = prepare_single_region_file(
            region_name=region_name,
            region_path=Path(bed_path),
            genome_fasta=genome_fasta,
            output_dir=output_dir,
            extract_to_bed3=extract_to_bed3,
        )

        processed_config["motif_prepared_inputs"][region_name] = prepared

    processed_path = config_path.with_name(f"{config_path.stem}.motif.processed.yaml")
    try:
        with open(processed_path, "w") as f:
            yaml.dump(processed_config, f, default_flow_style=False, sort_keys=False)
        logger.info(f"Motif processed config written: {processed_path}")
    except Exception as e:
        logger.error(f"Failed to write processed motif config: {e}")
        raise

    logger.info("Motif input preparation complete")
    return processed_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare BED and FASTA inputs for motif analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python prepare_motif_inputs.py --config config.yaml
  python prepare_motif_inputs.py --config config.yaml --log-level DEBUG
        """,
    )

    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to YAML configuration file",
    )

    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level",
    )

    args = parser.parse_args()
    logger.setLevel(getattr(logging, args.log_level))

    try:
        processed_config_path = prepare_motif_inputs(args.config)
        logger.info("✓ Motif input preparation successful")
        logger.info(f"✓ Use this config for downstream motif analysis: {processed_config_path}")
        print(f"\nMotif processed config: {processed_config_path}")
    except Exception as e:
        logger.error(f"✗ Motif input preparation failed: {e}", exc_info=True)
        raise
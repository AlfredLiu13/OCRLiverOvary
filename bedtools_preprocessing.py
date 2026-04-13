"""
Preprocessing module for OCR classification pipeline.
 
Purpose:
    Prepare raw ATAC-seq and HALPER mapping files for bedtools operations.
    
Tasks:
    - Unzip .gz files (preserving originals)
    - Extract first 3 columns (BED3 format)
    - Generate processed config pointing to cleaned files
 
Usage:
    python bedtools_preprocessing.py --config config.yaml
    
Output:
    - config.processed.yaml (updated configuration)
    - Cleaned BED3 files in bedtool_preprocess_output_dir/
    - Original files untouched
"""
 
import logging
import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
import gzip
import shutil
from typing import Dict, Any, Optional
 

#Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


# Configuration Class
@dataclass
class Config:
    """Configuration for preprocessing pipeline."""
    #Species + tissue being compared
    species_1: str      
    species_2: str
    tissue: str

    #Input ATAC-seq file paths
    species_1_peak_file: Path
    species_2_peak_file: Path

    #HALPER mapping file paths
    species_1_to_species_2: Optional[Path]
    species_2_to_species_1: Optional[Path]

    #TSS file paths
    species_1_tss_file: Path
    species_2_tss_file: Path

    #File ouput locations
    output_dir: Path    #Directory to store outputs
    temp_dir: Path      #Directory to store temporsry files and intermediate results --> May be deleted later

    #Check configuration automatically after intialization
    def __post_init__(self):
        """Validate configuration after initialization."""
        #Validate ATAC-seq files
        for path, name in [
            (self.species_1_peak_file, f"{self.species_1} peaks"),
            (self.species_2_peak_file, f"{self.species_2} peaks"),
        ]:
            if not self._file_exists(path):
                raise FileNotFoundError(f"{name} not found: {path}")
        
        #Validate mapping files (if provided)
        for path, name in [
            (self.species_1_to_species_2, f"Mapping {self.species_1}→{self.species_2}"),
            (self.species_2_to_species_1, f"Mapping {self.species_2}→{self.species_1}"),
        ]:
            if path is None:
                continue
            if not self._file_exists(path):
                raise FileNotFoundError(f"{name} not found: {path}")

        #Validate TSS files
        for path, name in [
            (self.species_1_tss_file, f"{self.species_1} TSS"),
            (self.species_2_tss_file, f"{self.species_2} TSS"),
        ]:
            if not path.exists():
                raise FileNotFoundError(f"{name} not found: {path}")
        
        #Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {self.output_dir}")

    @staticmethod
    def _file_exists(path: Path) -> bool:
        """Check if file exists (handles .gz variants)."""
        if path.exists():
            return True
        # Check for .gz variant
        if path.suffix != ".gz":
            gz_path = path.with_suffix(path.suffix + ".gz")
            if gz_path.exists():
                return True
        return False


#HELPER FUNCTIONS

def gunzip_keep(src: Path, dst: Path) -> None:
    """
    Unzip a .gz file without deleting the original.
    
    Parameters:
    -----------
    src : Path -> Path to .gz file
    dst : Path -> Path to output unzipped file
    
    Raises:
    -------
    FileNotFoundError -> If source file doesn't exist
    IOError -> If unzip operation fails
    """
    if not src.exists():
        raise FileNotFoundError(f"Source file not found: {src}")
    
    try:
        with gzip.open(src, "rb") as f_in:
            with open(dst, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.debug(f"Unzipped: {src} → {dst}")
    except Exception as e:
        logger.error(f"Failed to unzip {src}: {e}")
        raise

def ensure_unzipped(path: Path) -> Path:
    """
    Ensure file is unzipped, returning path to unzipped version.
    
    Parameters:
    -----------
    path : Path -> Path to file (may be .gz)
    
    Returns:
    --------
    Path -> Path to unzipped file
    """

    if path.suffix != ".gz":
        return path
    
    unzipped = path.with_suffix("")
    if unzipped.exists():
        logger.debug(f"File already unzipped: {unzipped}")
        return unzipped
    
    gunzip_keep(path, unzipped)
    return unzipped


def resolve_file_path(config_path: Path, key: str) -> Optional[Path]:
    """
    Resolve file path, checking for .gz variant if needed.
    
    Parameters:
    -----------
    config_path : Path -> Path from config file
    key : str -> Config key name (for logging)
    
    Returns:
    --------
    Optional[Path] -> Resolved path, or None if not found
    
    Raises:
    -------
    FileNotFoundError -> If file and .gz variant don't exist
    """
    config_path = Path(config_path)
    
    if config_path.exists():
        return config_path
    
    # Try .gz variant
    if config_path.suffix != ".gz":
        gz_path = config_path.with_suffix(config_path.suffix + ".gz")
        if gz_path.exists():
            logger.debug(f"{key}: Using .gz variant: {gz_path}")
            return gz_path
    
    raise FileNotFoundError(f"{key} not found: {config_path} (or .gz variant)")

#extract_bed3() --> Extracts the first 3 columns of a BED file
#Keeps only chrom, start, end
def extract_bed3(input_path: Path, output_path: Path):
    with open(output_path, "w") as out:
        subprocess.run(
            ["cut", "-f1-3", str(input_path)],
            stdout=out,
            check=True
        )


#load_config() --> Reads in a yaml file and intializes a config object
def load_config(config_path: Path, output_dir_key: str) -> Config:
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return Config(
        species_1=cfg["species_1"],
        species_2=cfg["species_2"],
        tissue=cfg["tissue"],

        output_dir=Path(cfg[output_dir_key]),
        temp_dir=Path(cfg["temp_dir"]),

        species_1_peak_file=Path(cfg["species_1_peak_file_cleaned"]),
        species_2_peak_file=Path(cfg["species_2_peak_file_cleaned"]),

        species_1_to_species_2=Path(cfg["species_1_to_species_2_cleaned"])
            if "species_1_to_species_2_cleaned" in cfg else None,
        species_2_to_species_1=Path(cfg["species_2_to_species_1_cleaned"])
            if "species_2_to_species_1_cleaned" in cfg else None,

        species_1_tss_file=Path(cfg["species_1_TSS_file"]),
        species_2_tss_file=Path(cfg["species_2_TSS_file"]),
    )


# Preprocessing
HALPER_KEYS = [
    "species_1_to_species_2",
    "species_2_to_species_1",
]

ATACSEQ_KEYS = [
    "species_1_peak_file",
    "species_2_peak_file",
]

"""
preprocess_config()
    Input: YAML file for the configuration

    Tasks:
        - Unzips HALPER files
        - Extracts first three columns from ATAC-SEQ + HALPER files
        - Writes new processed config

    Returns:
        Path to processed config
 """
def preprocess_config(config_path: Path) -> Path:

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    cleaned_dir = Path(config["bedtool_preprocess_output_dir"])
    cleaned_dir.mkdir(parents=True, exist_ok=True)


    #Process HALPER files
    for key in HALPER_KEYS:
        if key not in config:
            continue

        halper_path = Path(config[key])

        #If HALPER files are unzipped/missing -> Try .gz extension
        if not halper_path.exists():
            gz_path = halper_path.with_name(halper_path.name + ".gz")
            if gz_path.exists():
                halper_path = gz_path
            else:
                raise FileNotFoundError(f"{halper_path} not found")

        #Ensure HALPER files are unzipped
        halper_path = ensure_unzipped(halper_path)

        # Extract first three columns of BED files
        cleaned_file = cleaned_dir / halper_path.name
        if not cleaned_file.exists():
            extract_bed3(halper_path, cleaned_file)

        config[key + "_cleaned"] = str(cleaned_file)


    # Process ATACSEQ files
    for key in ATACSEQ_KEYS:
        peak_path = Path(config[key])
        ensure_exists(peak_path, "Peak file")

        cleaned_file = cleaned_dir / peak_path.name
        if not cleaned_file.exists():
            extract_bed3(peak_path, cleaned_file)

        config[key + "_cleaned"] = str(cleaned_file)


    # Write new config
    processed_path = config_path.with_suffix(".processed.yaml")
    with open(processed_path, "w") as f:
        yaml.dump(config, f)

    print(f"Processed config written to: {processed_path}")

    return processed_path

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess BED files for bedtools pipeline")
    parser.add_argument("--config", required=True, help="Path to YAML config file")

    args = parser.parse_args()

    processed_config_path = preprocess_config(Path(args.config))
    print(f"Preprocessing done! Processed config: {processed_config_path}")

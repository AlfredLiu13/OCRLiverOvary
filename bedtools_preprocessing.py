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


#Configuration Class
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
            (self.species_1_to_species_2, f"Mapping {self.species_1} -> {self.species_2}"),
            (self.species_2_to_species_1, f"Mapping {self.species_2} -> {self.species_1}"),
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


#Functions
def gunzip_keep(src: Path, dst: Path) -> None:
    """
    Unzip a .gz file without deleting the original.
    
    Inputs:
        src : Path -> Path to .gz file
        dst : Path -> Path to output unzipped file
    
    Errors:
        FileNotFoundError -> If source file doesn't exist
        IOError -> If unzip operation fails
    """
    if not src.exists():
        raise FileNotFoundError(f"Source file not found: {src}")
    
    try:
        with gzip.open(src, "rb") as f_in:
            with open(dst, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.debug(f"Unzipped: {src} -> {dst}")
    except Exception as e:
        logger.error(f"Failed to unzip {src}: {e}")
        raise

def ensure_unzipped(path: Path) -> Path:
    """
    Ensure file is unzipped, returning path to unzipped version.
    
    Inputs:
        path : Path -> Path to file (may be .gz)
    
    Ouputs:
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
    
    Inputs:
        config_path : Path -> Path from config file
        key : str -> Config key name (for logging)
    
    Ouputs:
        Optional[Path] -> Resolved path, or None if not found
    
    Errors:
        FileNotFoundError -> If file and .gz variant don't exist
    """
    config_path = Path(config_path)
    
    if config_path.exists():
        return config_path
    
    #Try .gz variant
    if config_path.suffix != ".gz":
        gz_path = config_path.with_suffix(config_path.suffix + ".gz")
        if gz_path.exists():
            logger.debug(f"{key}: Using .gz variant: {gz_path}")
            return gz_path
    
    raise FileNotFoundError(f"{key} not found: {config_path} (or .gz variant)")


def extract_bed3(input_path: Path, output_path: Path) -> int:
    """
    Extract first 3 columns (BED3 format) from input file. Keeps only chrom, start, end
    
    Inputs:
        input_path : Path -> Path to input BED file
        output_path : Path -> Path to output BED3 file
    
    Ouputs:
        int -> Number of lines in output file
    
    Errors:
        FileNotFoundError -> If input file doesn't exist
        subprocess.CalledProcessError -> If cut command fails
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    try:
        with open(output_path, "w") as out:
            result = subprocess.run(
                ["cut", "-f1-3", str(input_path)],
                stdout=out,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        # Count lines
        with open(output_path) as f:
            line_count = sum(1 for _ in f)
        
        logger.debug(f"Extracted BED3: {input_path} -> {output_path} ({line_count} lines)")
        return line_count
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract BED3 from {input_path}: {e.stderr}")
        raise

def sort_bed(input_path: Path, output_path: Path) -> int:
    """
    Sort a BED file by chromosome and position (required for bedtools).
    
    Inputs:
        input_path : Path -> Path to unsorted BED file
        output_path : Path -> Path to output sorted BED file
    
    Outputs:
        int -> Number of lines in output file
    
    Errors:
        FileNotFoundError -> If input file doesn't exist
        subprocess.CalledProcessError -> If sort command fails
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    try:
        with open(output_path, "w") as out:
            result = subprocess.run(
                ["sort", "-k1,1", "-k2,2n", str(input_path)],
                stdout=out,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        # Count lines
        with open(output_path) as f:
            line_count = sum(1 for _ in f)
        
        logger.debug(f"Sorted BED: {input_path} -> {output_path} ({line_count} lines)")
        return line_count
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to sort BED file {input_path}: {e.stderr}")
        raise

def read_yaml_config(config_path: Path) -> Dict[str, Any]:
    """
    Reads in a YAML configuration file.
    
    Inputs:
        config_path : Path -> Path to YAML config file
    
    Outputs:
        Dict[str, Any] -> Parsed YAML configuration
    
    Errors:
        FileNotFoundError -> If config file doesn't exist
        yaml.YAMLError -> If YAML parsing fails
    """
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


def preprocess_config(config_path: Path) -> Path:
    """
    Preprocess raw configuration file and BED files.
    
    Steps:
        1. Load YAML config
        2. Unzip any .gz files (preserving originals)
        3. Extract BED3 format (first 3 columns)
        4. Write processed config with paths to cleaned files
    
    Inputs:
        config_path : Path -> Path to original YAML config file
    
    Ouputs:
        Path -> Path to processed config file
    
    Errors:
        FileNotFoundError -> If required files don't exist
        yaml.YAMLError -> If config parsing fails
    """
    config_path = Path(config_path)
    logger.info(f"Starting preprocessing: {config_path}")
    
    #Load original config
    raw_config = read_yaml_config(config_path)
    
    #Create output directory
    cleaned_dir = Path(raw_config.get("bedtool_preprocess_output_dir", "temp/cleaned_bed"))
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {cleaned_dir}")
    
    #Process configuration
    processed_config = raw_config.copy()


    #Process ATAC-seq peak files
    logger.info("Processing ATAC-seq peak files...")
    
    for species_key, species_name in [
        ("species_1_peak_file", raw_config.get("species_1", "species1")),
        ("species_2_peak_file", raw_config.get("species_2", "species2")),
    ]:
        if species_key not in raw_config:
            logger.warning(f"Skipping {species_key} (not in config)")
            continue
        
        peak_path = Path(raw_config[species_key])
        logger.info(f"  Processing {species_name} peaks: {peak_path}")
        
        #Resolve path (handles .gz)
        try:
            peak_path = resolve_file_path(peak_path, species_key)
        except FileNotFoundError as e:
            logger.error(f"    Failed: {e}")
            raise
        
        #Ensure unzipped
        peak_path = ensure_unzipped(peak_path)
        
        #Extract BED3
        cleaned_file = cleaned_dir / peak_path.name
        if cleaned_file.exists():
            logger.info(f"    Cleaned file already exists: {cleaned_file}")
        else:
            try:
                line_count = extract_bed3(peak_path, cleaned_file)
                logger.info(f"    Extracted {line_count} lines -> {cleaned_file}")
            except Exception as e:
                logger.error(f"    Failed to extract BED3: {e}")
                raise
        
        #Sort BED file
        sorted_file = cleaned_dir / f"{cleaned_file.stem}.sorted{cleaned_file.suffix}"
        try:
            line_count = sort_bed(cleaned_file, sorted_file)
            logger.info(f"    Sorted {line_count} lines -> {sorted_file}")
            # Replace original with sorted version
            shutil.move(str(sorted_file), str(cleaned_file))
            logger.info(f"    Replaced with sorted version")
        except Exception as e:
            logger.error(f"    Failed to sort BED: {e}")
            raise

        #Update config with cleaned path
        processed_config[species_key + "_cleaned"] = str(cleaned_file)
    


    #Process HALPER mapping files
    logger.info("Processing HALPER mapping files...")
    
    mapping_keys = [
        ("species_1_to_species_2", "species_1 -> species_2"),
        ("species_2_to_species_1", "species_2 -> species_1"),
    ]
    
    for mapping_key, mapping_name in mapping_keys:
        if mapping_key not in raw_config:
            logger.debug(f"  Skipping {mapping_key} (not in config)")
            continue
        
        mapping_path = Path(raw_config[mapping_key])
        logger.info(f"  Processing {mapping_name} mapping: {mapping_path}")
        
        # Resolve path (handles .gz)
        try:
            mapping_path = resolve_file_path(mapping_path, mapping_key)
        except FileNotFoundError as e:
            logger.error(f"    Failed: {e}")
            raise
        
        #Ensure unzipped
        mapping_path = ensure_unzipped(mapping_path)
        
        #Extract BED3
        cleaned_file = cleaned_dir / mapping_path.name
        if cleaned_file.exists():
            logger.info(f"    Cleaned file already exists: {cleaned_file}")
        else:
            try:
                line_count = extract_bed3(mapping_path, cleaned_file)
                logger.info(f"    Extracted {line_count} lines -> {cleaned_file}")
            except Exception as e:
                logger.error(f"    Failed to extract BED3: {e}")
                raise

        #Sort BED file
        sorted_file = cleaned_dir / f"{cleaned_file.stem}.sorted{cleaned_file.suffix}"
        try:
            line_count = sort_bed(cleaned_file, sorted_file)
            logger.info(f"    Sorted {line_count} lines -> {sorted_file}")
            # Replace original with sorted version
            shutil.move(str(sorted_file), str(cleaned_file))
            logger.info(f"    Replaced with sorted version")
        except Exception as e:
            logger.error(f"    Failed to sort BED: {e}")
            raise
        
        #Update config with cleaned path
        processed_config[mapping_key + "_cleaned"] = str(cleaned_file)
    

    #Process TSS files
    logger.info("Processing TSS files...")
    
    tss_keys = [
        ("species_1_tss_file", raw_config.get("species_1", "species1")),
        ("species_2_tss_file", raw_config.get("species_2", "species2")),
    ]
    
    for tss_key, species_name in tss_keys:
        if tss_key not in raw_config:
            logger.warning(f"Skipping {tss_key} (not in config)")
            continue
        
        tss_path = Path(raw_config[tss_key])
        logger.info(f"  Processing {species_name} TSS: {tss_path}")
        
        # Resolve path (handles .gz)
        try:
            tss_path = resolve_file_path(tss_path, tss_key)
        except FileNotFoundError as e:
            logger.error(f"    Failed: {e}")
            raise
        
        #Ensure unzipped
        tss_path = ensure_unzipped(tss_path)
        
        #Extract BED3
        cleaned_file = cleaned_dir / tss_path.name
        
        if cleaned_file.exists():
            logger.info(f"    Cleaned file already exists: {cleaned_file}")
        else:
            try:
                line_count = extract_bed3(tss_path, cleaned_file)
                logger.info(f"    Extracted {line_count} lines -> {cleaned_file}")
            except Exception as e:
                logger.error(f"    Failed to extract BED3: {e}")
                raise

        #Sort BED file
        sorted_file = cleaned_dir / f"{cleaned_file.stem}.sorted{cleaned_file.suffix}"
        try:
            line_count = sort_bed(cleaned_file, sorted_file)
            logger.info(f"    Sorted {line_count} lines -> {sorted_file}")
            # Replace original with sorted version
            shutil.move(str(sorted_file), str(cleaned_file))
            logger.info(f"    Replaced with sorted version")
        except Exception as e:
            logger.error(f"    Failed to sort BED: {e}")
            raise

        #Update config
        processed_config[tss_key + "_cleaned"] = str(cleaned_file)


    #Write processed config    
    processed_path = config_path.with_stem(config_path.stem + ".processed")
    try:
        with open(processed_path, "w") as f:
            yaml.dump(processed_config, f, default_flow_style=False, sort_keys=False)
        logger.info(f"Processed config written: {processed_path}")
    except Exception as e:
        logger.error(f"Failed to write processed config: {e}")
        raise
    
    logger.info(f"Preprocessing complete!")
    return processed_path
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess BED files for OCR classification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python bedtools_preprocessing.py --config config.yaml
  python bedtools_preprocessing.py --config config.yaml --log-level DEBUG
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
        help="Logging level (default: INFO)",
    )
    
    args = parser.parse_args()
    
    #Set log level
    logger.setLevel(getattr(logging, args.log_level))
    
    try:
        processed_config_path = preprocess_config(args.config)
        logger.info(f"✓ Preprocessing successful")
        logger.info(f"✓ Use processed config for classification: {processed_config_path}")
        print(f"\nProcessed config: {processed_config_path}")
    except Exception as e:
        logger.error(f"✗ Preprocessing failed: {e}", exc_info=True)
        exit(1)
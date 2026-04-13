#!/usr/bin/env python3

"""
Pipeline poriton to classify OCRs as promoters/enhancers and compare regulatory states across species.

Usage:
    ./classify_ocr.py --config /path/to/config.yaml
    
SLURM Script:
    #!/bin/bash
    #SBATCH -J classification
    #SBATCH -p RM-shared
    #SBATCH -N 1
    #SBATCH -n 2
    #SBATCH -t 12:00:00
    #SBATCH --mem-per-cpu=2000M 
    #SBATCH -o logs/classification_%j.out 
    #SBATCH -e logs/classification_%j.err
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=alfredl@andrew.cmu.edu
    module load anaconda3
    conda activate my_bio_env
    ./classify_ocr.py --config /path/to/config.yaml
"""

from pathlib import Path
from pybedtools import BedTool #Available in bioconda environment
from main import run_classification
from typing import Dict, Tuple, Optional
import logging
import argparse

#Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

#Helper Functions
def classifyOcrPromotersEnhancers(
    ocr_bed_path: str,
    tss_bed_path: str,
    output_prefix: str,
    promoter_distance: int = 2000  #Might be worth switching to 5k?
) -> Dict[str, str]:
    """
    Classify OCRs as promoters or enhancers based on distance to nearest TSS.
    Promoters defined as OCRs within promoter_distance bp of a TSS.
    Enhancers are OCRs beyond that distance.
    
    Inputs:
        ocr_bed_path : str -> Path to the OCR BED file (peaks).
        tss_bed_path : str -> Path to the TSS BED file (1bp per TSS).
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_promoters.bed
            - {output_prefix}_enhancers.bed
            - {output_prefix}_nearestTSS.bed
        promoter_distance : int, optional -> Distance threshold for promoter classification (default: 2000 bp).
            Consider 5000 bp for more lenient classification.
    
    Outputs:
        Dict[str, str] -> Dictionary with keys 'promoters', 'enhancers', 'nearest' containing paths to the output BED files.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    
    ocr_bed_path = Path(ocr_bed_path)
    tss_bed_path = Path(tss_bed_path)

    #Validate input files
    if not ocr_bed_path.exists():
        raise FileNotFoundError(f"OCR BED file not found: {ocr_bed_path}")
    if not tss_bed_path.exists():
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_path}")
    
    logging.info(f"Processing OCRs: {ocr_bed_path}")
    
    try:
        #Load BED files
        ocr = BedTool(ocr_bed_path)
        tss = BedTool(tss_bed_path)

        #Validate non-empty
        if len(ocr) == 0:
            raise ValueError("OCR BED file is empty")
        if len(tss) == 0:
            raise ValueError("TSS BED file is empty")
        
        logger.info(f"Loaded {len(ocr)} OCRs and {len(tss)} TSSs")

        #Annotate nearest TSS + distance
        ocr_with_genes = ocr.closest(tss, d=True)
        nearest_tss_path = f"{output_prefix}_nearestTSS.bed"
        ocr_with_genes.saveas(nearest_tss_path)
        logger.info(f"Saved annotated OCRs with nearest TSS: {nearest_tss_path}")
    
        #Classify OCRs into promoters by distance
        promoters_path = f"{output_prefix}_promoters.bed"
        ocr_promoters = ocr_with_genes.filter(lambda x: int(x[-1]) <= promoter_distance).saveas(promoters_path)
        logger.info(
            f"Classified {len(ocr_promoters)} promoters (<=  {promoter_distance} bp): "
            f"{promoters_path}"
        )

        #Classify OCRs into enhancers by distance
        enhancers_path = f"{output_prefix}_enhancers.bed"
        ocr_enhancers = ocr_with_genes.filter(lambda x: int(x[-1]) > promoter_distance).saveas(f"{output_prefix}_enhancers.bed")
        logger.info(
            f"Classified {len(ocr_enhancers)} enhancers (> {promoter_distance} bp): "
            f"{enhancers_path}"
        )

        logging.info(f"Finished processing {ocr_bed_path}")


        return {
            "promoters": promoters_path,
            "enhancers": enhancers_path,
            "nearest": nearest_tss_path
        }
    
    except Exception as e:
        logger.error(f"Error processing {ocr_bed_path}: {e}")
        raise

def classifyConservedRegions(
        conserved_bed_path: str,
        tss_bed_path: str,
        output_prefix: str,
        promoter_distance: int = 2000
) -> Tuple[str, str]:
    """
    Classify conserved/mapped regions as promoters or enhancers using TSS distance.
    
    Inputs:
        conserved_bed_path : str -> Path to the conserved/mapped OCR BED file.
        tss_bed_path : str -> Path to the TSS BED file (1bp per TSS).
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_promoters.bed
            - {output_prefix}_enhancers.bed
            - {output_prefix}_TSS.bed (annotated)
        promoter_distance : int, optional -> Distance threshold for promoter classification (default: 2000 bp).
    
    Outputs:
        Tuple[str, str] -> Paths to promoters and enhancers BED files.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    conserved_bed_path = Path(conserved_bed_path)
    tss_bed_path = Path(tss_bed_path)

    if not conserved_bed_path.exists():
        raise FileNotFoundError(f"Conserved BED file not found: {conserved_bed_path}")
    if not tss_bed_path.exists():
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_path}")

    logging.info(f"Processing conserved regions: {conserved_bed_path}")

    try:
        conserved = BedTool(conserved_bed_path)
        tss = BedTool(tss_bed_path)

        if len(conserved) == 0:
            raise ValueError("Conserved BED file is empty")
        if len(tss) == 0:
            raise ValueError("TSS BED file is empty")
        
        logger.info(f"Loaded {len(conserved)} conserved regions and {len(tss)} TSSs")

        #Annotate distance to nearest TSS
        tss_path = f"{output_prefix}_TSS.bed"
        annotated = conserved.closest(tss, d=True)
        annotated.saveas(tss_path)
        logger.info(f"Saved annotated conserved regions: {tss_path}")

        #Split into promoters
        promoters_path = f"{output_prefix}_promoters.bed"
        promoters = annotated.filter(lambda x: int(x[-1]) <= promoter_distance).saveas(promoters_path)
        logger.info(
            f"Classified {len(promoters)} conserved promoters "
            f"(<= {promoter_distance} bp): {promoters_path}"
        )

        #Split into enhancers
        enhancers_path = f"{output_prefix}_enhancers.bed"
        enhancers = annotated.filter(lambda x: int(x[-1]) >  promoter_distance).saveas(enhancers_path)
        logger.info(
            f"Classified {len(enhancers)} conserved enhancers "
            f"(> {promoter_distance} bp): {enhancers_path}"
        )

        return promoters_path, enhancers_path
    
    except Exception as e:
        logger.error(f"Error processing {conserved_bed_path}: {e}")
        raise


def findSharedElements(
        mapped_file_path: str,
        native_file_path: str,
        output_file_path: str
) -> str:
    """
    Identify conserved regions that overlap native regulatory elements in target species.
        - Finds mapped OCRs from one species that intersect with native
        OCRs detected in another species, identifying true conserved regulatory elements.
    
    Inputs:
        mapped_file_path : str -> Path to mapped/conserved OCR BED file (from source species, mapped to target).
        native_file_path : str -> Path to native OCR BED file (detected in target species).
        output_file_path : str -> Path for output BED file containing overlapping regions.
    
    Ouputs:
        str -> Path to the output BED file with shared elements.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    mapped_file_path = Path(mapped_file_path)
    native_file_path = Path(native_file_path)

    if not mapped_file_path.exists():
        raise FileNotFoundError(f"Mapped BED file not found: {mapped_file_path}")
    if not native_file_path.exists():
        raise FileNotFoundError(f"Native BED file not found: {native_file_path}")

    logging.info(f"Finding shared elements between {mapped_file_path} and {native_file_path}")

    try:
        mapped = BedTool(mapped_file_path)
        native = BedTool(native_file_path)

        if len(mapped) == 0:
            logger.warning("Mapped BED file is empty")
        if len(native) == 0:
            logger.warning("Native BED file is empty")        

        logger.info(f"Loaded {len(mapped)} mapped regions and {len(native)} native regions")

        #Find intersection regions
        shared = mapped.intersect(native, u=True) #u=True -> Keeps only regions that intersect
        shared.saveas(output_file_path)

        logger.info(f"Found {len(shared)} shared elements: {output_file_path}")

        return output_file_path

    except Exception as e:
        logger.error(f"Error finding shared elements: {e}")
        raise

#Main wrapper function
def run_classification(config_path: Path) -> None:
    """
    Main pipeline orchestrator. Loads config and runs classification steps.
    
    Expected YAML config structure:
    ```
    species:
      species1:
        ocr_bed: path/to/species1_ocrs.bed
        tss_bed: path/to/species1_tss.bed
        output_prefix: results/species1
      species2:
        ocr_bed: path/to/species2_ocrs.bed
        tss_bed: path/to/species2_tss.bed
        output_prefix: results/species2
    
    mapping:
      species1_to_species2: path/to/species1_ocrs_mapped_to_species2.bed
      species2_to_species1: path/to/species2_ocrs_mapped_to_species1.bed
    
    parameters:
      promoter_distance: 2000
    ```
    
    Inputs:
        config_path : Path -> Path to YAML configuration file.
    """
    try:
        import yaml
    except ImportError:
        logger.error(
            "PyYAML not found. Install with: pip install pyyaml"
        )
        raise
    
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    logger.info(f"Loading configuration from {config_path}")
    
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading config: {e}")
        raise
    
    #Extract parameters
    species_config = config.get("species", {})
    mapping_config = config.get("mapping", {})
    params = config.get("parameters", {})
    
    promoter_distance = params.get("promoter_distance", 2000)
    logger.info(f"Using promoter distance threshold: {promoter_distance} bp")
    
    #Step 1: Classify native OCRs in each species
    classification_results = {}
    for species_name, species_data in species_config.items():
        logger.info(f"\n{'='*70}")
        logger.info(f"STEP 1: Classifying native OCRs for {species_name}")
        logger.info(f"{'='*70}")
        
        result = classifyOcrPromotersEnhancers(
            ocr_bed_path=species_data["ocr_bed"],
            tss_bed_path=species_data["tss_bed"],
            output_prefix=species_data["output_prefix"],
            promoter_distance=promoter_distance,
        )
        classification_results[species_name] = result
    
    #Step 2: Classify mapped/conserved OCRs
    logger.info(f"\n{'='*70}")
    logger.info("STEP 2: Classifying mapped/conserved OCRs")
    logger.info(f"{'='*70}")
    
    mapping_results = {}
    for mapping_name, mapping_data in mapping_config.items():
        # Parse mapping direction (e.g., "species1_to_species2")
        parts = mapping_name.split("_to_")
        if len(parts) != 2:
            logger.warning(f"Skipping malformed mapping name: {mapping_name}")
            continue
        
        source_species, target_species = parts
        logger.info(f"\nProcessing mapping: {source_species} → {target_species}")
        
        #Get target species TSS for annotation
        if target_species not in species_config:
            logger.warning(
                f"Target species '{target_species}' not found in config. Skipping."
            )
            continue
        
        target_tss = species_config[target_species]["tss_bed"]
        output_prefix = f"{species_config[target_species]['output_prefix']}_mapped_from_{source_species}"
        
        #Classify mapped regions
        promoters_path, enhancers_path = classifyConservedRegions(
            conserved_bed_path=mapping_data,
            tss_bed_path=target_tss,
            output_prefix=output_prefix,
            promoter_distance=promoter_distance,
        )
        
        mapping_results[mapping_name] = {
            "promoters": promoters_path,
            "enhancers": enhancers_path,
        }
    
    #Step 3: Find shared regulatory elements across species
    logger.info(f"\n{'='*70}")
    logger.info("STEP 3: Finding shared regulatory elements")
    logger.info(f"{'='*70}")
    
    shared_results = {}
    for mapping_name, mapping_data in mapping_config.items():
        parts = mapping_name.split("_to_")
        if len(parts) != 2:
            continue
        
        source_species, target_species = parts
        
        if source_species not in species_config or target_species not in species_config:
            logger.warning(f"Skipping mapping with missing species config: {mapping_name}")
            continue
        
        logger.info(f"\nFinding shared elements: {mapping_name}")
        
        # Find overlap between mapped and native OCRs
        mapped_ocrs = mapping_data
        native_ocrs = species_config[target_species]["ocr_bed"]
        output_file = (
            f"{species_config[target_species]['output_prefix']}"
            f"_shared_from_{source_species}.bed"
        )
        
        shared_path = findSharedElements(
            mapped_file_path=mapped_ocrs,
            native_file_path=native_ocrs,
            output_file_path=output_file,
        )
        
        shared_results[mapping_name] = shared_path
    
    #Print summary
    logger.info(f"\n{'='*70}")
    logger.info("CLASSIFICATION PIPELINE COMPLETE")
    logger.info(f"{'='*70}")
    logger.info("Classification Results:")
    for species, files in classification_results.items():
        logger.info(f"  {species}:")
        for file_type, path in files.items():
            logger.info(f"    {file_type}: {path}")
    
    if mapping_results:
        logger.info("Mapped/Conserved Classification:")
        for mapping, files in mapping_results.items():
            logger.info(f"  {mapping}:")
            for file_type, path in files.items():
                logger.info(f"    {file_type}: {path}")
    
    if shared_results:
        logger.info("Shared Elements:")
        for mapping, path in shared_results.items():
            logger.info(f"  {mapping}: {path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers and compare across species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ./classify_ocr.py --config config.yaml
  ./classify_ocr.py --config /path/to/config.yaml --log-level DEBUG
        """,
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",)

    args = parser.parse_args()

    #Set log level
    logger.setLevel(getattr(logging, args.log_level))
    
    try:
        run_classification(args.config)
        logger.info("Pipeline finished successfully")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        exit(1)

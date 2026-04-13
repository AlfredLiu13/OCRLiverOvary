from pathlib import Path
from pybedtools import BedTool #Available in bioconda environment
from main import run_classification
import logging
import argparse

"""
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


"""
classifyOcrPromotersEnhancers() --> Classify OCRs as promoters or enhancers and assign nearest TSS.
Paramaeters:
- ocrBedPath (str): Path to the OCR BED file (peaks).
- tssBedPath (str): Path to the TSS BED file (1bp per TSS).
- genomeFile (str): Genome file (chromosome sizes) for pybedtools slop.
- outputPrefix (str): Prefix for output files. Outputs will be saved as:
    - {outputPrefix}_promoters.bed
    - {outputPrefix}_enhancers.bed
    - {outputPrefix}_nearestTSS.bed
 """
def classifyOcrPromotersEnhancers(
    ocrBedPath: str,
    tssBedPath: str,
    outputPrefix: str,
    promoterDistance: int = 2000  #Might be worth switching to 5k?
):
    
    logging.info(f"Processing OCRs: {ocrBedPath}")
    
    # Load BED files
    ocr = BedTool(ocrBedPath)
    tss = BedTool(tssBedPath)

    #Annotate nearest TSS + distance
    ocrWithGenes = ocr.closest(tss, d=True)
    ocrWithGenes.saveas(f"{outputPrefix}_nearestTSS.bed")
    
    #Classify OCRs into promoters by distance
    ocrPromoters = ocrWithGenes.filter(
        lambda x: int(x[-1]) <= promoterDistance
    ).saveas(f"{outputPrefix}_promoters.bed")

    #Classify OCRs into enhancers by distance
    ocrEnhancers = ocrWithGenes.filter(
        lambda x: int(x[-1]) > promoterDistance
    ).saveas(f"{outputPrefix}_enhancers.bed")

    logging.info(f"Finished processing {ocrBedPath}")
    
    return {
        "promoters": f"{outputPrefix}_promoters.bed",
        "enhancers": f"{outputPrefix}_enhancers.bed",
        "nearest": f"{outputPrefix}_nearestTSS.bed"
    }

#classifyConservedRegions() --> Assigns promoter/enhancer labels to mapped (conserved) regions using TSS distance.
def classifyConservedRegions(
        conservedBedPath: str,
        tssBedPath: str,
        outputPrefix: str,
        promoterDistance: int = 2000
        ):

    logging.info(f"Processing conserved regions: {conservedBedPath}")

    conserved = BedTool(conservedBedPath)
    tss = BedTool(tssBedPath)

    # Annotate distance
    annotated = conserved.closest(tss, d=True)
    annotated.saveas(f"{outputPrefix}_TSS.bed")

    # Split promoter/enhancer
    promoters = annotated.filter(lambda x: int(x[-1]) <= promoterDistance).saveas(
        f"{outputPrefix}_promoters.bed"
    )
    enhancers = annotated.filter(lambda x: int(x[-1]) >  promoterDistance).saveas(
        f"{outputPrefix}_enhancers.bed"
    )

    return promoters, enhancers

#findSharedElements() --> Identify conserved regions that overlap native regulatory elements in the target species.
def findSharedElements(
        mappedFilePath: str,
        nativeFilePath: str,
        outputFilePath: str
        ):

    logging.info(f"Finding shared elements")

    mapped = BedTool(mappedFilePath)
    native = BedTool(nativeFilePath)

    shared = mapped.intersect(native, u=True)
    shared.saveas(outputFilePath)

    return outputFilePath

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    run_classification(Path(args.config))

from pybedtools import BedTool #Available in bioconda environment

#Probably better to create a wrapper function:
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
    outputPrefix: str
    promoter_distance: int = 2000  #Might be worth switching to 5k?
):
    
    # Load BED files
    ocr = BedTool(ocrBedPath)
    tss = BedTool(tssBedPath)

    #Annotate nearest TSS + distance
    ocrWithGenes = ocr.closest(tss, d=True)
    ocrWithGenes.saveas(f"{outputPrefix}_nearestTSS.bed")
    
    #Classify OCRs into promoters by distance
    ocrPromoters = ocrWithGenes.filter(
        lambda x: int(x[-1]) <= promoter_distance
    ).saveas(f"{outputPrefix}_promoters.bed")

    #Classify OCRs into enhancers by distance
    ocrEnhancers = ocrWithGenes.filter(
        lambda x: int(x[-1]) > promoter_distance
    ).saveas(f"{outputPrefix}_enhancers.bed")

    print(f"Finished processing {ocrBedPath}")
    print(f"Promoters (≤{promoter_distance}bp): {outputPrefix}_promoters.bed")
    print(f"Enhancers (>{promoter_distance}bp): {outputPrefix}_enhancers.bed")

    return {
        "promoters": f"{outputPrefix}_promoters.bed",
        "enhancers": f"{outputPrefix}_enhancers.bed",
        "nearest": f"{outputPrefix}_nearestTSS.bed"
    }

#classifyConservedRegions() --> Assigns promoter/enhancer labels to mapped (conserved) regions using TSS distance.
def classifyConservedRegions(conserved_bed, tss_bed, output_prefix):

    conserved = BedTool(conserved_bed)
    tss = BedTool(tss_bed)

    # Annotate distance
    annotated = conserved.closest(tss, d=True)
    annotated.saveas(f"{output_prefix}_TSS.bed")

    # Split promoter/enhancer
    promoters = annotated.filter(lambda x: int(x[-1]) <= 5000).saveas(
        f"{output_prefix}_promoters.bed"
    )
    enhancers = annotated.filter(lambda x: int(x[-1]) > 5000).saveas(
        f"{output_prefix}_enhancers.bed"
    )

    return promoters, enhancers

#findSharedElements() --> Identify conserved regions that overlap native regulatory elements in the target species.
def findSharedElements(mapped_file, native_file, output_file):

    mapped = BedTool(mapped_file)
    native = BedTool(native_file)

    shared = mapped.intersect(native, u=True)
    shared.saveas(output_file)

    return output_file

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    run_classification(Path(args.config))

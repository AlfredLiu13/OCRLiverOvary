#main entry point for the pipeline
from pathlib import Path
from bedtools_preprocessing import load_config
from classification import classifyOcrPromotersEnhancers

def run_quality_control():
    # check quality of human and mouse datasets
    

    
    pass


def run_alignment():
    # map mouse peaks to human genome using HALPER
    # output should be orthologous regions (.bed)
    pass


def run_regulatory_comparison():
    # compare open chromatin regions between species
    # identify:
    # - shared (open in both)
    # - human-specific
    # - mouse-specific
    pass


#run_classification() -> Classifies OCR peaks into promoters and enhancers based on distance to TSS (±2kb)
def run_classification(config_path: Path):
    cfg = load_config(config_path, "output_dir")

    # Output prefix naming
    human_prefix = cfg.output_dir / f"{cfg.species_1}_{cfg.tissue}"
    mouse_prefix = cfg.output_dir / f"{cfg.species_2}_{cfg.tissue}"

    # Run for species 1
    classifyOcrPromotersEnhancers(
        ocrBedPath=str(cfg.species_1_peak_file),
        tssBedPath=str(cfg.species_1_tss_file),
        genomeFile=str(cfg.species_1_genome_fasta),
        outputPrefix=str(human_prefix)
    )

    # Run for species 2
    classifyOcrPromotersEnhancers(
        ocrBedPath=str(cfg.species_2_peak_file),
        tssBedPath=str(cfg.species_2_tss_file),
        genomeFile=str(cfg.species_2_genome_fasta),
        outputPrefix=str(mouse_prefix)
    )


def run_motif_analysis():
    # find transcription factor motifs in regions
    # using HOMER
    pass


def run_enrichment_analysis():
    # find biological processes (GO terms)
    # using GREAT / HOMER
    pass


def run_benchmarking():
    # test pipeline using toy or small real dataset
    # check if outputs make sense
    pass


def main():
    # run full pipeline step by step

    run_quality_control()
    run_alignment()
    run_regulatory_comparison()
    run_classification()
    run_motif_analysis()
    run_enrichment_analysis()
    run_benchmarking()

    print("pipeline finished")


if __name__ == "__main__":
    main()

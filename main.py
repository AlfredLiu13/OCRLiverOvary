#main entry point for the pipeline

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


def run_classification():
    # classify regions into promoters and enhancers
    # based on distance to TSS (±2kb)
    pass


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
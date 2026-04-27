# LiverOCRAnalysis


## Overview

This repository contains a bioinformatics pipeline for analyzing and comparing open chromatin regions (OCRs) in the tissue of two species. It accomplishes the following tasks:
- Mapping OCRs between species  
- Identifying shared and species-specific regions  
- Classifying regions into promoters and enhancers  
- Identifying transcription factor binding motifs  
- Performing enrichment analysis  

---

## Structure

The repository is organized based on the main analysis steps:

- `data_qc/` – quality control of datasets  
- `alignment/` – mapping OCRs across species  
- `classification/` – promoter vs enhancer classification + shared and species-specific regions identification
- `motif_analysis/` – transcription factor motif analysis  
- `enrichment_analysis/` – biological process enrichment 
- `results/` – output files  

Other files:
- `main.py`
- `bedtools_preprocessing.py`
- `requirements.txt`

---

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/BioinformaticsDataPracticum2026/LiverOCRAnalysis.git
cd LiverOCRAnalysis
```

### 2.Install Python Dependencies

You can Install LiverOCRAnalysis with:

#### Option A: Using pip
```bash
pip install -r requirements.txt 
```

#### Option B: Using conda 
```bash
conda create -n ocr_analysis python=3.8
conda activate ocr_analysis
pip install -r requirements.txt
conda install -c bioconda bedtools pybedtools
conda install -c conda-forge readr -y
```

### 3. Install External Tools
The pipeline requires the following external tools. Install them according to your system and environment:

#### [BEDTools](https://bedtools.readthedocs.io/en/latest/) (v2.31+)
```bash
# Using conda
conda install -c bioconda bedtools

# Or download from: https://bedtools.readthedocs.io/en/latest/
```

#### [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing)
```bash
git clone https://github.com/pfenninglab/halLiftover-postprocessing.git

# Follow installation instructions in the HALPER repository
```

#### [HOMER](http://homer.ucsd.edu/homer/motif/) (v5.0+)
```bash
# Download and install from: http://homer.ucsd.edu/homer/motif/
cd ~/tools
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer

# Add HOMER to PATH
export PATH=$PATH:~/tools/homer/bin
```

#### [rGREAT](https://jokergoo.github.io/rGREAT/) (v2.0+)
Install in R:
```R
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("rGREAT")
```

---


## Usage

Each step can be run independently or you can run the full pipeline

###  Input Requirements

Before running LiverOCRAnalysis, you should have:

- OCR peak files (BED / narrowPeak)
- TSS annotation files (BED)
- HAL alignment file (for cross-species mapping)

---

### To run individual modules

```bash
#Alignment
python main.py alignment \
  --human-peaks data/human_peaks.narrowPeak \
  --mouse-peaks data/mouse_peaks.narrowPeak \
  --hal-file data/10plusway.hal

#Classification
python main.py classification \
  --config classification/sample_config.yaml

#Motif Analysis
python main.py motif

#Annotate peaks
python main.py annotate \
  --bed-dir results \
  --outdir results/annotated

#GREAT analysis
python main.py great \
  --bed-dir results \
  --outdir results/great \
  --species hg38
  ```

---

### To Run Full LiverOCRAnalysis

```bash
python main.py full \
  --human-peaks data/human_peaks.narrowPeak \
  --mouse-peaks data/mouse_peaks.narrowPeak \
  --hal-file data/10plusway.hal \
  --config classification/sample_config.yaml
```

---

### Flags Description

| Module | Flag | Description |
|-------|------|------------|
| **Alignment** | `--human-peaks` | Human OCR peak file (BED / narrowPeak) |
|  | `--mouse-peaks` | Mouse OCR peak file |
|  | `--hal-file` | HAL alignment file for cross-species mapping |
|  | `--outdir` | Output directory (default: `results/alignment`) |
|  | `--min-len` | Minimum length for mapped regions |
|  | `--protect-dist` | Distance around peak summit to preserve |
|  | `--max-frac` | Maximum allowed size change after mapping |
|  | `--preserve` | Extra peak columns to keep |
|  | `--hal-liftover-path` | Path to `halLiftover` binary |
|  | `--keep-chr-prefix` | Filter chromosomes by prefix |
| **Classification** | `--config` | YAML file defining OCRs, TSS, and mappings |
|  | `--log-level` | Logging verbosity (INFO, DEBUG, etc.) |
|  | `--script` | Path to classification script |
| **Motif** | `--script` | Motif analysis script |
| **Annotate** | `--bed-dir` | Directory containing BED files |
|  | `--outdir` | Output directory (default: `results/annotated`) |
|  | `--genome` | Genome build or FASTA (default: hg38) |
|  | `--homer-bin` | Path to `annotatePeaks.pl` |
|  | `--gtf` / `--gff3` | Custom annotation file |
|  | `--beds` | Specific BED files to process |
| **GREAT** | `--bed-dir` | Directory containing BED files |
|  | `--outdir` | Output directory (default: `results/great`) |
|  | `--species` | Genome assembly (default: hg38) |
|  | `--beds` | Specific BED files to analyze |
|  | `--rscript-bin` | Path to Rscript executable |
|  | `--script` | GREAT batch script |
|  | `--r-script` | GREAT R script |
| **Full** | `--human-peaks`, `--mouse-peaks`, `--hal-file` | Alignment inputs |
|  | `--config` | Classification config |

---

### Outputs

| Module | Output File / Artifact | Location | Format | Description |
|--------|----------------------|----------|--------|-------------|
| **Alignment** | `*.HALPER.narrowPeak.gz` | `results/alignment_results/human_to_mouse/`<br>`results/alignment_results/mouse_to_human/` | gzipped narrowPeak | Cross-species mapped OCR peaks (main HALPER output) |
| **Alignment** | `*.halLiftover.tFile.bed.gz` | Same as above | gzipped BED | Intermediate file of target-genome regions before HALPER filtering |
| **Alignment** | `*.halLiftover.sFile.bed.gz` | Same as above | gzipped BED | Intermediate file of peak summit positions before HALPER filtering |
| **Classification** | `{species}_all_promoters.bed`<br>`{species}_all_enhancers.bed` | `results/classification_results/raw_results/` | BED | All OCRs per species split into promoters (≤2000 bp from TSS) and enhancers (>2000 bp) |
| **Classification** | `shared_promoters.bed`<br>`shared_enhancers.bed` | `results/classification_results/raw_results/` | BED | OCRs with orthologs open in both species, classified as promoters or enhancers |
| **Classification** | `{species}_specific.bed`<br>`{species}_specific_promoters.bed`<br>`{species}_specific_enhancers.bed` | `results/classification_results/raw_results/` | BED | OCRs open in one species but not the other, with promoter/enhancer sub-classifications |
| **Classification** | `summary_table.csv` | `results/classification_results/results_analyzed/` | CSV | Per-category counts and percentages for all classified OCR sets |
| **Classification** | `promoter_enhancer_classification.png`<br>`ortholog_status.png`<br>`conservation_comparison.png`<br>`percentage_comparison.png` | `results/classification_results/results_analyzed/` | PNG | Bar charts and comparison plots of OCR classification and conservation |
| **Enrichment Analysis** | `gobp.csv` | `results/enrichment/great/{dataset}/` | CSV | GO Biological Process enrichment results (p-values, fold enrichment) per OCR category |
| **Enrichment Analysis** | `metadata.txt` | `results/enrichment/great/{dataset}/` | TXT | Region counts and filtering details for the GREAT submission |
| **Enrichment Analysis** | `great_summary.tsv` | `results/enrichment/summary/` | TSV | Top 10 enriched GO BP terms per dataset consolidated across all OCR categories |
| **Enrichment Analysis** | `barplot_{dataset}.png` | `results/enrichment/plots/` | PNG | Bar plots of top enriched GO BP terms per OCR category (one per dataset) |
| **Enrichment Analysis** | `comparison_heatmap.png` | `results/enrichment/plots/` | PNG | Heatmap comparing GO BP enrichment scores across all OCR categories |
| **Motif Analysis** | `homerResults.html` | `results/findmotifs_results/{region_type}/` | HTML | Interactive report of de novo motif discovery results |
| **Motif Analysis** | `knownResults.html`<br>`knownResults.txt` | `results/findmotifs_results/{region_type}/` | HTML / TSV | Known transcription factor motif enrichment report with statistics |
| **Motif Analysis** | `homerMotifs.all.motifs`<br>`homerMotifs.motifs8`<br>`homerMotifs.motifs10`<br>`homerMotifs.motifs12` | `results/findmotifs_results/{region_type}/` | HOMER motif | De novo motif position weight matrices at 8, 10, and 12 bp lengths |
| **Motif Analysis** | `known{N}.motif`<br>`known{N}.logo.svg` | `results/findmotifs_results/{region_type}/knownResults/` | PWM / SVG | Individual known motif PWM files and sequence logo images |

---

## Citation
To cite this repository, please copy the following:

__Hamda Al Hosani, Alfred Liu, Samridhi Makkar, Bhanvi Paliwal (2026).__ _LiverOCRAnalysis_. 03-713: Bioinformatics Data Integration Practicum, Carnegie Mellon Univeristy.

---

## Contact

For questions, issues, or contributions, please reach out to:

| Name | Email |
|------|-------|
| Alfred Liu | alfredl@andrew.cmu.edu |
| Hamda Al Hosani | halhosan@andrew.cmu.edu |
| Samridhi Makkar | smakkar@andrew.cmu.edu |
| Bhanvi Paliwal | bhanvip@andrew.cmu.edu |

---

## License

This project is part of the 03-713: Bioinformatics Data Integration Practicum at Carnegie Mellon University.

---

## Acknowledgments

- BEDTools: Quinlan & Hall (2010)
- HALPER: [Pfenning Lab](https://pfenninglab.org/)
- HOMER: Heinz et al. (2010)
- rGREAT: Gu et al. (2016)

---
**Last Updated:** 2026
**Version:** 1.0

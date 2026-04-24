# LiverOCRAnalysis


## Overview

This repository contains a bioinformatics pipeline for analyzing and comparing open chromatin regions (OCRs) in the tissue of two species. It accomplishes the following tasks:
- Mapping OCRs between species  
- Identifying shared and species-specific regions  
- Classifying regions into promoters and enhancers  
- Identifying transcription factor binding motifs  
- Performing enrichment analysis  

---
## Citation
To cite this repository, please copy the following:

__Hamda Al Hosani, Alfred Liu, Samridhi Makkar, Bhanvi Paliwal (2026).__ _LiverOCRAnalysis_. 03-713: Bioinformatics Data Integration Practicum, Carnegie Mellon Univeristy.

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

---

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

### To Run Full LiverOCRAnalysis

```bash
python main.py full \
  --human-peaks data/human_peaks.narrowPeak \
  --mouse-peaks data/mouse_peaks.narrowPeak \
  --hal-file data/10plusway.hal \
  --config classification/sample_config.yaml
```

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

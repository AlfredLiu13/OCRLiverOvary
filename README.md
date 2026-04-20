# LiverOCRAnalysis


## Overview

This repository contains a pipeline for analyzing and comparing open chromatin regions (OCRs) in the tissue of two species. It accomplishes the following tasks:
- Mapping OCRs between species  
- Identifying shared and species-specific regions  
- Classifying regions into promoters and enhancers  
- Identifying transcription factor binding motifs  
- Performing enrichment analysis  

---
## Citation
To cite this repository, please copy the following:

Hamda Al Hosani, Alfred Liu, Samridhi Makkar, Bhanvi Paliwal (2026). _LiverOCRAnalysis_. 03-713: Bioinformatics Data Integration Practicum, Carnegie Mellon Univeristy.

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

# Dependencies

You can Install LiverOCRAnalysis with:

Option1:
```bash
pip install -r requirements.txt 
```

Option2:
```bash
conda install ...
```

# External tools required:

* [BEDTools](https://bedtools.readthedocs.io/en/latest/) (v2.31)
* [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing)
* [HOMER](http://homer.ucsd.edu/homer/motif/) (v5)
* [rGREAT](https://jokergoo.github.io/rGREAT/) (v2.0)

---

## Usage
```bash    
```


---
## Demo 


---
## Contact
If you have any further questions please reach out to:
- alfredl@andrew.cmu.edu
- halhosan@andrew.cmu.edu
- smakkar@andrew.cmu.edu
- bhanvip@andrew.cmu.edu

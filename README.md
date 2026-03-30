# LiverOCRAnalysis

## Overview

This repository contains a pipeline for analyzing open chromatin regions (OCRs) in human and mouse and comparing regulatory activity across species.

The pipeline includes:
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
- `regulatory_comparison/` – identifying shared and species-specific regions  
- `classification/` – promoter vs enhancer classification  
- `motif_analysis/` – transcription factor motif analysis  
- `enrichment_analysis/` – biological process enrichment  
- `benchmarking/` – testing and validation  
- `results/` – output files  

Other files:
- `main.py` – pipeline structure  
- `requirements.txt` – Python dependencies  

---

## Dependencies

Install Python dependencies with:

```bash
pip install -r requirements.txt 
```

## External tools required:

* BEDTools
* HALPER
* HOMER
* GREAT / rGREAT
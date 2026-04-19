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
- `classification/` – promoter vs enhancer classification + identifying shared and species-specific regions
- `motif_analysis/` – transcription factor motif analysis  
- `enrichment_analysis/` – biological process enrichment  
- `benchmarking/` – testing and validation  
- `results/` – output files  

Other files:
- `main.py` 
- `requirements.txt`  

---

## Dependencies

You can Install LiverOCRAnalysis with:

Option1:
```bash
pip install -r requirements.txt 
```

Option2:
```bash
conda install ...
```

## External tools required:

* BEDTools
* HALPER
* HOMER
* GREAT / rGREAT
* yaml
* argparse
* logging
* typing

## Usage
```bash    
```



## Demo 



## Contact
If you have any further questions please reach out to:
- alfredl@andrew.cmu.edu
- halhosan@andrew.cmu.edu
- smakkar@andrew.cmu.edu
- bhanvip@andrew.cmu.edu

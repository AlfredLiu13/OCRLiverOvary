from pathlib import Path


INPUT_BED_DIR = Path("results/classification_results/raw_results")
GREAT_OUTPUT_DIR = Path("results/enrichment/great")
SUMMARY_OUTPUT_DIR = Path("results/enrichment/summary")

GREAT_R_SCRIPT = Path("enrichment_analysis/great_online.R")


BED_FILES = {
    "human_shared": {
        "file": "human_shared.bed",
        "genome": "hg38",
    },
    "human_specific": {
        "file": "human_specific.bed",
        "genome": "hg38",
    },
    "mouse_shared": {
        "file": "mouse_shared.bed",
        "genome": "mm10",
    },
    "mouse_specific": {
        "file": "mouse_specific.bed",
        "genome": "mm10",
    },
    "shared_promoters": {
        "file": "shared_promoters.bed",
        "genome": "mm10",
    },
    "shared_enhancers": {
        "file": "shared_enhancers.bed",
        "genome": "mm10",
    },
}
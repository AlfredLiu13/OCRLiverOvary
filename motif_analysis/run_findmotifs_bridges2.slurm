#!/bin/bash
#SBATCH --job-name=findmotifs_liver
#SBATCH --account=bio230007p
#SBATCH --partition=RM
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/jet/home/bpaliwal/LiverOCRAnalysis/LiverOCRAnalysis/logs/findmotifs_%j.out
#SBATCH --error=/jet/home/bpaliwal/LiverOCRAnalysis/LiverOCRAnalysis/logs/findmotifs_%j.err

module load homer

BASE=/jet/home/bpaliwal/LiverOCRAnalysis/LiverOCRAnalysis
GENOME=/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/hg38.fa
PREP=/jet/home/bpaliwal/homer_preparsed/hg38

mkdir -p "$PREP"
mkdir -p "$BASE/motif_analysis/findmotifs_results"

findMotifsGenome.pl \
  "$BASE/results/classification_results/raw_results/human_promoters.bed" \
  "$GENOME" \
  "$BASE/motif_analysis/findmotifs_results/human_promoters" \
  -size given -len 8,10,12 -p 8 -mask \
  -preparsedDir "$PREP"

findMotifsGenome.pl \
  "$BASE/results/classification_results/raw_results/human_enhancers.bed" \
  "$GENOME" \
  "$BASE/motif_analysis/findmotifs_results/human_enhancers" \
  -size given -len 8,10,12 -p 8 -mask \
  -preparsedDir "$PREP"

findMotifsGenome.pl \
  "$BASE/results/classification_results/raw_results/human_mapped_from_mouse_promoters.bed" \
  "$GENOME" \
  "$BASE/motif_analysis/findmotifs_results/human_mapped_from_mouse_promoters" \
  -size given -len 8,10,12 -p 8 -mask \
  -preparsedDir "$PREP"

findMotifsGenome.pl \
  "$BASE/results/classification_results/raw_results/human_mapped_from_mouse_enhancers.bed" \
  "$GENOME" \
  "$BASE/motif_analysis/findmotifs_results/human_mapped_from_mouse_enhancers" \
  -size given -len 8,10,12 -p 8 -mask \
  -preparsedDir "$PREP"

findMotifsGenome.pl \
  "$BASE/results/classification_results/raw_results/human_shared_from_mouse.bed" \
  "$GENOME" \
  "$BASE/motif_analysis/findmotifs_results/human_shared_from_mouse" \
  -size given -len 8,10,12 -p 8 -mask \
  -preparsedDir "$PREP"







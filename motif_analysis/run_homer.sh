#!/bin/bash
#SBATCH -J enrichment
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smakkar@andrew.cmu.edu

INPUT_BED="filepath"
OUTPUT_DIR="/jet/home/smakkar/enrichment/results"

export PATH=$PATH:/jet/home/smakkar/enrichment/homer/bin/

mkdir -p $OUTPUT_DIR

#subject to change for the reference genome and the motif size
findMotifsGenome.pl $INPUT_BED hg18 $OUTPUT_DIR -size 200 -mask -p 2
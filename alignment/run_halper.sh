#!/bin/bash
#SBATCH -J halper_map
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH -o logs/halper_%j.out
#SBATCH -e logs/halper_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smakkar@andrew.cmu.edu

# Load modules
module load anaconda3

# Activate HAL conda environment
source activate hal

# Set file paths
HUMAN_PEAKS="/ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"
MOUSE_PEAKS="/ocean/projects/bio230007p/ikaplow/MouseAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"
OUTPUT_DIR="/jet/home/smakkar/atac_align/results"
HAL_FILE="/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal"

# Make output directory
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/human_to_mouse
mkdir -p $OUTPUT_DIR/mouse_to_human
mkdir -p logs

# Unzip peak files
echo "Unzipping human peaks..."
gunzip -c $HUMAN_PEAKS > $OUTPUT_DIR/human_liver.narrowPeak

echo "Unzipping mouse peaks..."
gunzip -c $MOUSE_PEAKS > $OUTPUT_DIR/mouse_liver.narrowPeak

echo "Starting HALPER mapping..."

# Run mappings in parallel
bash /jet/home/smakkar/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
    -b $OUTPUT_DIR/human_liver.narrowPeak \
    -o $OUTPUT_DIR/human_to_mouse \
    -s "Human" \
    -t "Mouse" \
    -c $HAL_FILE &

bash /jet/home/smakkar/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
    -b $OUTPUT_DIR/mouse_liver.narrowPeak \
    -o $OUTPUT_DIR/mouse_to_human \
    -s "Mouse" \
    -t "Human" \
    -c $HAL_FILE &

wait

echo "HALPER mapping completed"

#!/bin/bash
#SBATCH -J halper_map
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -o logs/halper_%j.out
#SBATCH -e logs/halper_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smakkar@andrew.cmu.edu


module load anaconda3

source activate hal

HUMAN_PEAKS="/ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/rep2/SRR13439655_1.trim.merged.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak.gz"
MOUSE_PEAKS="/ocean/projects/bio230007p/ikaplow/MouseAtac/Liver/peak/rep1/SRR8119852_1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak.gz"
OUTPUT_DIR="/jet/home/smakkar/atac_align/results"
HAL_FILE="/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal"

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/human_to_mouse
mkdir -p $OUTPUT_DIR/mouse_to_human

echo "Unzipping human peaks..."
gunzip -c $HUMAN_PEAKS > $OUTPUT_DIR/human_liver.narrowPeak

echo "Unzipping mouse peaks..."
gunzip -c $MOUSE_PEAKS > $OUTPUT_DIR/mouse_liver.narrowPeak

echo "Starting HALPER mapping..."

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

echo "HALPER mapping jobs submitted"


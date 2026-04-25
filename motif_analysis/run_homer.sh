#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smakkar@andrew.cmu.edu

export PATH="$PATH:/jet/home/smakkar/enrichment/homer/bin"

INPUT_DIR="/jet/home/smakkar/enrichment"
BASE_OUT="/jet/home/smakkar/enrichment/results"
# New directory for HOMER's background files
PREPARSED_DIR="/jet/home/smakkar/enrichment/preparsed_genome"

PROMOTER_BED="$INPUT_DIR/human_promoters.bed"
ENHANCER_BED="$INPUT_DIR/human_enhancers.bed"

PROMOTER_OUT="$BASE_OUT/human_promoters"
ENHANCER_OUT="$BASE_OUT/human_enhancers"

GENOME_PATH="/ocean/projects/bio230007p/ikaplow"
GENOME_FILE="$GENOME_PATH/HumanGenomeInfo/hg38.fa"

mkdir -p "$PROMOTER_OUT"
mkdir -p "$ENHANCER_OUT"
mkdir -p "$PREPARSED_DIR"

echo "Job started: $(date)"

echo "Checking input files..."
[ -f "$PROMOTER_BED" ] || { echo "ERROR: Missing $PROMOTER_BED"; exit 1; }
[ -f "$ENHANCER_BED" ] || { echo "ERROR: Missing $ENHANCER_BED"; exit 1; }

module load homer

echo "Launching promoter analysis..."
findMotifsGenome.pl "$PROMOTER_BED" "$GENOME_FILE" "$PROMOTER_OUT" -size 200 -mask -p 2 -preparsedDir "$PREPARSED_DIR" \
  > "$PROMOTER_OUT/run.log" 2>&1 &
PID1=$!

echo "Launching enhancer analysis..."
findMotifsGenome.pl "$ENHANCER_BED" "$GENOME_FILE" "$ENHANCER_OUT" -size 200 -mask -p 2 -preparsedDir "$PREPARSED_DIR" \
  > "$ENHANCER_OUT/run.log" 2>&1 &
PID2=$!

echo "Waiting for promoter job to finish..."
wait $PID1
STATUS1=$?
echo "Promoter job finished with exit code: $STATUS1"

echo "Waiting for enhancer job to finish..."
wait $PID2
STATUS2=$?
echo "Enhancer job finished with exit code: $STATUS2"

if [ $STATUS1 -ne 0 ] || [ $STATUS2 -ne 0 ]; then
    echo "ERROR: One or more jobs failed"
    echo "Check logs:"
    echo "  $PROMOTER_OUT/run.log"
    echo "  $ENHANCER_OUT/run.log"
    exit 1
fi

echo "All motif analyses completed successfully at $(date)"

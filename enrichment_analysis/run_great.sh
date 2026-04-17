# This script runs GREAT analysis on the given classification bed files and generates GREAT output files for each classification.
#!/bin/bash
#input - classification bed files (eg:- human_shared_from_mouse.bed, human_specific.bed, mouse_shared_from_human.bed, mouse_specific.bed, etc)
#output - GREAT output files (eg:- human_shared_from_mouse_great_output.tsv, human_specific_great_output.tsv, mouse_shared_from_human_great_output.tsv, mouse_specific_great_output.tsv, etc)


#Install rGREAT package in R if not already installed
Rscript -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("rGREAT")'

#Bash wrapper to run GREAT analysis on the given classification bed files

#!/bin/bash
# great_rlocal.sh
# Usage: ./great_rlocal.sh input.bed hg38 output_prefix

BED_FILE=$1
SPECIES=$2
OUTPUT=$3

Rscript --vanilla << EOF
library(rGREAT)
library(rtracklayer)
gr <- import("${BED_FILE}")
result <- great(gr, species="${SPECIES}")                              #eg:- hg38 for human, mm10 for mouse
write.table(as.data.frame(result), file="${OUTPUT}.tsv", sep="\t", row.names=FALSE, quote=FALSE)
EOF

echo "Results saved to ${OUTPUT}.tsv"
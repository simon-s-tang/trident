#!/bin/bash
# Create the conda environment (needs conda >= 23 or mamba).
#
#   bash create_environment.sh              # creates `trident`
#   bash create_environment.sh <name>       # under another name
#
# workflow/envs/trident.yaml holds the conda packages; the APA Python
# packages are then installed from GitHub by workflow/envs/trident.post-deploy.sh.
set -euo pipefail
cd "$(dirname "$0")"
eval "$(conda shell.bash hook)"

ENV_NAME=${1:-trident}

conda env create -n "$ENV_NAME" -f workflow/envs/trident.yaml

conda activate "$ENV_NAME"
export PYTHONNOUSERSITE=1
bash workflow/envs/trident.post-deploy.sh

# Check every tool runs
echo "--- tool versions"
snakemake --version
fastp --version 2>&1
fastqc --version
STAR --version
featureCounts -v 2>&1 | grep featureCounts
samtools --version | head -1
htseq-count --version
multiqc --version
rmats.py --version
kallisto version
Rscript -e 'cat("limma", format(packageVersion("limma")), "\n")'
python -c "import apatools, bulkanalysis, genomemanager, kallistomanager, scanpy; print('APA python packages: OK')"

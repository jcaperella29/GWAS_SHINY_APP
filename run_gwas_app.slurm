#!/bin/bash
#SBATCH --job-name=gwas_shiny
#SBATCH --output=gwas_shiny.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --partition=debug

echo "📦 Launching GWAS app..."
singularity exec gwas_app.sif Rscript -e "options(shiny.port=3838, shiny.host='0.0.0.0'); shiny::runApp('/app')"

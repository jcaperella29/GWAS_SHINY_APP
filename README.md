# 🧬 JCAP GWAS Shiny App

A complete, interactive platform for performing QC, GWAS, ML modeling, gene mapping, enrichment, and predictive analysis — powered by **Shiny** and ready for **HPC deployment via Singularity**.

---

## 🔗 Table of Contents

1. [Features](#features)
2. [Screenshots](#screenshots)
3. [Tech Stack](#tech-stack)
4. [Installation](#installation)
   - [Run Locally](#run-locally)
   - [Run with Singularity (HPC)](#run-with-singularity-hpc)
5. [How to Use](#how-to-use)
6. [Downloadable Outputs](#downloadable-outputs)
7. [Common Errors](#common-errors)
8. [Credits & Contact](#credits--contact)

---

## 🚀 Features

- 📤 Upload your own **VCF**, **phenotype**, and **covariate** data
- ✅ Perform **SNP Quality Control** with MAF, HWE, and Call Rate thresholds
- 📊 Explore data via **PCA** and **UMAP** projections
- 🔬 Run **GWAS** with Bonferroni option and adjustable p-value slider
- 🧪 View **Q-Q** and **Manhattan** plots with chromosome/position filters
- 🧬 Map SNPs to nearest **human genes**
- 🧭 Perform **Pathway Enrichment Analysis** with barplots for KEGG, GO, and Reactome
- 🤖 Train/test a **Random Forest model** using significant SNPs
- 📈 Inspect **model predictions**, **AUC**, **sensitivity/specificity**, and **variable importance**
- 📉 Run **power analysis** with effect size and power curves
- 🧾 Export all results as CSV


---

## 🧱 Tech Stack

- **R 4.x** via [`rocker/shiny`](https://hub.docker.com/r/rocker/shiny)
- Shiny + `data.table`, `DT`, `plotly`, `ranger`, `mlr3`, `VariantAnnotation`, `enrichR`, `pwr`, `ROCR`
- Containerized using **Singularity**
- Designed for HPC **cluster-friendly** use

---

## 🔧 Installation

### 📦 Run Locally (with R)

1. Clone this repo:
   in `bash
   git clone https://github.com/jcaperella29/GWAS_SHINY_APP
   cd GWAS_SHINY_APP

2.Launch the app
in R

shiny::runApp()

Make sure required R packages are installed (listed in app.R)

###⚗️ Run with Singularity (HPC)

1.Build container:

in bash
sudo singularity build gwas_app.sif Singularity.def

2.Run inside HPC environment:
in bash

singularity exec --bind $PWD:/workspace gwas_app.sif R -e "shiny::runApp('/srv/shiny-server/app', host='0.0.0.0', port=3838)"
3.Forward port or launch on compute node via sbatch script

script

🧭 How to Use
Upload Data
Upload VCF, phenotype (CSV), and optional covariates (CSV)

Run QC
Click Run QC button after adjusting MAF/HWE/Call Rate thresholds

Explore Data
Visualize sample space with PCA and UMAP

Power Analysis
Hit Run Power to estimate statistical power and see the power curve

Run GWAS
Adjust the p-value slider and choose Bonferroni if needed, then click Run GWAS

Visualize GWAS
Use Q-Q plot and Manhattan plot (adjust chr/position filters as needed)

Map SNPs to Genes
Click Map Genes and view table + plot

Run Enrichment
Perform pathway enrichment and view top results via Barplot

Run Random Forest
Train/test classifier on final hits. View predictions, metrics, and SNP importance

Download Outputs
Use sidebar download buttons to get CSVs of all results

📥 Downloadable Outputs
GWAS full results

Covariate associations

Final phenotype-only hits

SNP-to-gene annotation

Pathway enrichment table

Power analysis summary

Random forest predictions, metrics, and feature importance

❌ Common Errors
Error	Cause	Solution
Not enough samples in one or both groups	Trait imbalance	Ensure both classes are represented
Can't read output 'rfPredictions'	RF model not trained	Click Run Random Forest first
Gene mapping failed	Invalid SNP coordinates or gene DB issue	Check VCF chr formatting (should match hg19)
Enrichment failed	Too few gene symbols	Ensure gene mapping produced valid symbols
Package not found	Singularity not built properly	Rebuild container with correct R packages

🙋‍♂️ Credits & Contact
🔬 App by John Caperella

📬 Contact: @jcaperella29
“Stealing the distorted desires of legacy pipelines.” 🎩

Ready for GWAS evolution.












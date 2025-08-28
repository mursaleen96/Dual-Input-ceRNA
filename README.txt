# ceRNA Pipeline

[![Python](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/)
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Snakemake-based bioinformatics pipeline for identifying and analyzing competing endogenous RNA (ceRNA) networks from RNA-seq data. This pipeline integrates quality control, normalization, feature engineering, machine learning prediction, statistical validation, and network visualization to discover ceRNA interactions.

## 🧬 Overview

The ceRNA pipeline processes RNA-seq count data to identify lncRNA-miRNA-mRNA regulatory networks using:

- **Database Integration**: miRTarBase 2025 and starBase/ENCORI for miRNA-target interactions
- **Machine Learning**: XGBoost classifier for ceRNA triplet prediction  
- **Statistical Validation**: Mediation analysis with Sobel testing
- **Network Analysis**: GraphML/SIF/CSV exports for Cytoscape visualization
- **Interactive Reports**: HTML reports with embedded network visualizations

### Key Features

- ✅ Handles both synthetic and real RNA-seq datasets
- ✅ Automated miRNA name standardization
- ✅ Robust error handling for empty datasets
- ✅ Multiple output formats (GraphML, SIF, CSV)
- ✅ Interactive HTML reports with Plotly networks
- ✅ Cytoscape-compatible network files
- ✅ Excel-friendly CSV exports

## 🚀 Installation

### Prerequisites

- **Operating System**: Linux/macOS (tested on Ubuntu 20.04/WSL2)
- **Python**: 3.9 or higher
- **Conda/Mamba**: For environment management

### Setup

1. **Clone the repository**
git clone https://github.com/your-username/cerna-pipeline.git
cd cerna-pipeline

text

2. **Create conda environment**
conda env create -f environment.yml
conda activate cerna-pipeline

text

3. **Install additional dependencies**
pip install -r requirements.txt

text

## 📊 Database Setup

### Required Database Files

The pipeline requires three main database files. Download them from the following sources:

#### 1. miRTarBase (miRNA-mRNA interactions)

**Source**: [miRTarBase 2025](https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/php/index.php)

- Download the MTI file (Excel or text format)
- Save as: `databases/miRTarBase_MTI.txt`

#### 2. starBase/ENCORI (miRNA-lncRNA interactions)

**Source**: [starBase/ENCORI](https://rnasysu.com/encori/)

- Navigate to "miRNA Target" > "Download"
- Select: miRNA-lncRNA interactions for human
- Save as: `databases/starBase_miRNA_lncRNA.txt`

#### 3. Ensembl Gene Annotations

**Source**: [Ensembl FTP](https://www.ensembl.org/info/data/ftp/index.html)

- Download GTF annotation file
- Convert to CSV format with columns: gene_id, gene_name, gene_biotype
- Save as: `gene_annotation.csv`

### Process Database Files

After downloading the raw database files:

python modules/download_databases.py

text

## 🏃 Quick Start

### Basic Usage

Run with your RNA-seq counts file
python cerna_pipeline_main.py --input your_counts.csv --threads 4

For WSL/Windows users (add latency handling)
python cerna_pipeline_main.py --input your_counts.csv --threads 4 --latency-wait 60

text

### Standardize miRNA Names (if needed)

If your input file has non-standard miRNA names:

Edit the script to specify your input/output files
python scripts/mirna_name_fix.py

text

## 📁 Input Files

### Required Files

| File | Description | Format |
|------|-------------|--------|
| **Raw Counts CSV** | Gene expression matrix (genes×samples) | CSV with gene IDs as row names |
| **miRTarBase_MTI.txt** | miRNA-mRNA interactions | Tab-separated: miRNA, Target |
| **starBase_miRNA_lncRNA.txt** | miRNA-lncRNA interactions | Tab-separated: miRNA, lncRNA |
| **gene_annotation.csv** | Gene annotations from GTF | CSV: gene_id, gene_name, gene_biotype |

### Input Data Format

Your RNA-seq counts file should be structured as:
gene_id,sample1,sample2,sample3,...
ENSG00000000001,245,312,189,...
hsa-miR-21-5p,1523,1876,1234,...
ENSG00000000002,67,89,45,...
...

text

## 📊 Output Files

All results are saved in the `results/` directory:

### Core Outputs

| File | Description | Use Case |
|------|-------------|----------|
| `validated_triplets.csv` | Statistically validated ceRNA triplets | Publication-ready results |
| `cerna_analysis_report.html` | Interactive HTML report | Quick visualization |
| `cerna_network.graphml` | Network in GraphML format | Cytoscape import |
| `cerna_network.sif` | Network in SIF format | Cytoscape import |
| `cerna_network_nodes.csv` | Network nodes (Excel-compatible) | Further analysis |
| `cerna_network_edges.csv` | Network edges (Excel-compatible) | Further analysis |

## ⚙️ Configuration

Edit `config/config.yaml` to customize pipeline parameters:

Normalization and filtering
low_count_threshold: 5 # Minimum counts per gene
sample_frac_threshold: 0.8 # Proportion of samples to keep gene
normalization_method: CPM # CPM normalization method

Statistical thresholds
confidence_threshold: 0.7 # ML confidence threshold
mediation_pval_cutoff: 0.05 # Mediation analysis p-value

Analysis parameters
random_seed: 42 # Reproducibility

text

## 🔧 Troubleshooting

### Common Issues

#### Empty Triplets/No Results
- **Cause**: Gene name mismatch between input and databases
- **Solution**: Use the miRNA standardization script or check gene IDs

#### Network Not Displaying in HTML Report
- **Cause**: Missing Plotly dependencies or JavaScript issues
- **Solution**: Ensure `plotly` is installed; try opening in different browser

#### Filesystem Errors (WSL/Windows)
- **Cause**: File system latency in mounted drives
- **Solution**: Add `--latency-wait 60` to pipeline command
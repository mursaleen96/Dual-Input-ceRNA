# Dual-Input ceRNA Discovery Pipeline

A comprehensive computational pipeline for discovering competing endogenous RNA (ceRNA) interactions using dual RNA-sequencing inputs: RNA-seq and miRNA-seq data.

## Overview

### What is ceRNA?

Competing endogenous RNAs (ceRNAs) are RNA molecules that regulate each other's expression by competing for shared microRNA (miRNA) binding. This creates a regulatory network where:

- **lncRNAs** (long non-coding RNAs) and **mRNAs** compete for the same miRNA binding sites
- **miRNAs** act as mediators in the regulatory relationship
- The interaction forms **ceRNA triplets**: lncRNA-miRNA-mRNA regulatory units

### Pipeline Purpose

This pipeline identifies biologically relevant ceRNA triplets by:

1. **Integrating dual RNA-seq inputs**: Both mRNA/lncRNA expression (RNA-seq) and miRNA expression (miRNA-seq)
2. **Applying computational prediction**: Machine learning models trained on expression correlation and sequence features
3. **Statistical validation**: Mediation analysis to confirm indirect regulatory effects
4. **Network construction**: Building ceRNA interaction networks for downstream analysis

## Pipeline Architecture

The pipeline is implemented as a **Snakemake workflow** with modular components:

```
Input Data (RNA-seq + miRNA-seq)
    ↓
[1] QC & Normalization
    ↓
[2] Database Download ←→ [3] Feature Engineering
    ↓                        ↓
[4] ML Training ←────────── Features
    ↓                        ↓
[5] Triplet Prediction ←── Models
    ↓
[6] Statistical Validation
    ↓
[7] Network Analysis
    ↓
[8] Report Generation
    ↓
Final Results (Networks + Report)
```

## Workflow Stages

### Stage 1: Quality Control & Normalization (`qc_normalization.py`)

**Purpose**: Prepare and normalize dual RNA-seq inputs for analysis

**Key Functions**:
- **Sample alignment**: Ensures RNA-seq and miRNA-seq have matching samples
- **Gene ID standardization**: Handles Ensembl, Entrez, and gene symbol formats
- **Low-expression filtering**: Removes genes with insufficient expression (configurable thresholds)
- **CPM normalization**: Counts per million followed by log2 transformation
- **Gene annotation**: Integrates gene biotype information (mRNA vs lncRNA classification)

**Outputs**:
- `results/norm_counts.csv`: Normalized expression matrix (genes × samples)
- `results/sample_metadata.csv`: Sample information and QC metrics

### Stage 2: Database Download (`download_databases.py`)

**Purpose**: Retrieve known miRNA-target interaction databases

**Databases Downloaded**:
- **miRTarBase**: Experimentally validated miRNA-mRNA interactions
- **LncBase**: Experimentally validated miRNA-lncRNA interactions  
- **Gene Annotation**: Ensembl gene annotations with biotype information

**Outputs**:
- `databases/miRTarBase.txt`: miRNA-mRNA interaction database
- `databases/LncBase.txt`: miRNA-lncRNA interaction database
- `databases/gene_annotation.csv`: Gene annotation mappings

### Stage 3: Feature Engineering (`feature_engineering.py`)

**Purpose**: Generate candidate ceRNA triplets and compute predictive features

**Triplet Generation Process**:
1. For each miRNA in the expression data:
   - Find validated mRNA targets (from miRTarBase)
   - Find validated lncRNA targets (from LncBase)
   - Form all possible lncRNA-miRNA-mRNA combinations

**Features Computed** (per triplet):

**Expression Correlation Features**:
- Pearson correlation: lncRNA ↔ mRNA
- Pearson correlation: lncRNA ↔ miRNA  
- Pearson correlation: mRNA ↔ miRNA
- Partial correlation: lncRNA ↔ mRNA (controlling for miRNA)

**ceRNA-Specific Features**:
- **SPONGE score**: Measures how much the miRNA "sponges" the lncRNA-mRNA relationship
- **Mediation sensitivity**: Quantifies indirect effect strength

**Sequence Features** (placeholders for future extension):
- miRNA response element (MRE) counts
- Seed match binding energy
- Subcellular localization scores

**Quality Filters Applied**:
- Correlation significance thresholds (p-value < 0.01)
- Minimum SPONGE score (> 0.1)
- Expression variance filters

**Output**:
- `results/features.pkl`: Feature matrix for all candidate triplets

### Stage 4: Machine Learning Training (`ml_training.py`)

**Purpose**: Train predictive models to identify high-confidence ceRNA interactions

**Model Architecture**:
- **Algorithm**: XGBoost (Extreme Gradient Boosting)
- **Label Generation**: Top 30% of SPONGE scores as positive examples
- **Feature Selection**: Numeric features only (excludes gene identifiers)
- **Handling Missing Data**: Imputation for robust model training

**Training Process**:
1. Load engineered features from Stage 3
2. Generate binary labels based on SPONGE score distribution
3. Select and preprocess numeric feature columns
4. Train XGBoost classifier with cross-validation
5. Save trained model for prediction

**Output**:
- `results/models.pkl`: Trained machine learning model

### Stage 5: Triplet Prediction (`predict_triplets.py`)

**Purpose**: Apply trained models to predict ceRNA triplet probabilities

**Prediction Process**:
1. Load trained model and feature matrix
2. Apply model to generate prediction probabilities
3. Rank triplets by prediction confidence
4. Filter for high-confidence predictions

**Output**:
- `results/predicted_triplets.csv`: Predicted ceRNA triplets with confidence scores

### Stage 6: Statistical Validation (`statistical_validation.py`)

**Purpose**: Validate predicted triplets using rigorous statistical mediation analysis

**Mediation Analysis Method**:
Uses **Sobel test** to assess indirect effects in the ceRNA regulatory pathway:

```
lncRNA → miRNA → mRNA
   ↓        ↓
Path a   Path b
   ↓        ↓
   Indirect Effect = a × b
```

**Statistical Tests Performed**:
- **Path a**: lncRNA → miRNA regression
- **Path b**: miRNA → mRNA regression (controlling for lncRNA)
- **Sobel test**: Tests significance of mediation effect (a × b)
- **Multiple testing correction**: FDR (False Discovery Rate) adjustment

**Mediation Types Classified**:
- **Complete mediation**: Indirect effect significant, direct effect non-significant
- **Partial mediation**: Both indirect and direct effects significant
- **No mediation**: Indirect effect non-significant

**Quality Control**:
- Variance filters (removes low-variance expressions)
- Sample size requirements
- Outlier detection and handling

**Outputs**:
- `results/validated_triplets.csv`: All triplets with mediation analysis results
- `results/validated_triplets_ranked.csv`: Ranked by statistical significance

### Stage 7: Network Analysis (`network_analysis.py`)

**Purpose**: Construct ceRNA interaction networks from validated triplets

**Network Construction**:
- **Nodes**: Genes (lncRNAs, mRNAs, miRNAs) with expression and biotype attributes
- **Edges**: ceRNA interactions with statistical significance weights
- **Filtering**: Top-ranked interactions only (typically FDR < 0.05)

**Network Metrics Computed**:
- **Degree centrality**: Number of connections per gene
- **Betweenness centrality**: Importance as network bridge
- **Clustering coefficient**: Local network density
- **Hub identification**: Highly connected regulatory genes

**Output Formats**:
- `results/cerna_network.graphml`: NetworkX format for computational analysis
- `results/cerna_network.sif`: Simple interaction format for Cytoscape
- `results/centrality_scores.csv`: Network centrality metrics
- `results/cerna_network_nodes.csv`: Node attributes
- `results/cerna_network_edges.csv`: Edge attributes

### Stage 8: Report Generation (`generate_report.py`)

**Purpose**: Create comprehensive interactive analysis report

**Report Sections**:

1. **Executive Summary**:
   - Total triplets discovered
   - Statistical significance summary
   - Network topology metrics

2. **Quality Control Metrics**:
   - Sample alignment results
   - Expression filtering statistics
   - Database coverage analysis

3. **Statistical Results**:
   - P-value distribution plots
   - Effect size distributions
   - Multiple testing correction results

4. **Network Visualization**:
   - Interactive network plots (Plotly)
   - Centrality analysis plots
   - Hub gene identification

5. **Top Discoveries Table**:
   - Highest-ranked ceRNA triplets
   - Statistical significance metrics
   - Biological annotations

**Output**:
- `results/cerna_analysis_report.html`: Interactive HTML report

## Input Requirements

### Required Files

1. **RNA-seq counts matrix** (`--rna-input`):
   - Format: CSV file with genes as rows, samples as columns
   - Must include both mRNAs and lncRNAs
   - Row names: Gene identifiers (Ensembl/Entrez/Symbols)
   - Values: Raw read counts (not normalized)

2. **miRNA-seq counts matrix** (`--mirna-input`):
   - Format: CSV file with miRNAs as rows, samples as columns  
   - Row names: miRNA identifiers (miRBase format preferred)
   - Values: Raw read counts (not normalized)
   - **Critical**: Sample names must match RNA-seq file

### Configuration File (Optional)

Create `config/config.yaml` for custom parameters:

```yaml
# Expression filtering thresholds
min_count: 1          # Minimum read count
min_samples: 2        # Minimum samples with min_count

# Statistical thresholds
corr_cutoff: 0.5      # Correlation significance threshold
pval_cutoff: 0.01     # P-value significance threshold
sponge_cutoff: 0.1    # SPONGE score minimum threshold

# Gene annotation
gene_anno_path: "databases/gene_annotation.csv"
```

## Usage

### Installation

1. **Install Conda/Mamba** (if not already installed)

2. **Create environment**:
   ```bash
   conda env create -f environment.yml
   conda activate cerna-pipeline
   ```

3. **Install additional dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

### Basic Usage

```bash
python cerna_pipeline_main.py \
    --rna-input /path/to/rna_counts.csv \
    --mirna-input /path/to/mirna_counts.csv \
    --threads 8
```

### Advanced Usage

```bash
python cerna_pipeline_main.py \
    --rna-input /path/to/rna_counts.csv \
    --mirna-input /path/to/mirna_counts.csv \
    --config config/custom_config.yaml \
    --threads 16 \
    --latency-wait 30
```

### Parameters

- `--rna-input`: Path to RNA-seq counts CSV file (**required**)
- `--mirna-input`: Path to miRNA-seq counts CSV file (**required**)  
- `--config`: Path to YAML configuration file (default: `config/config.yaml`)
- `--threads`: Number of CPU cores to use (default: 4)
- `--latency-wait`: File system latency wait time in seconds (default: 10)

## Output Files

### Key Results

1. **`results/validated_triplets_ranked.csv`**:
   - Final ranked ceRNA triplets
   - Columns: lncRNA, miRNA, mRNA, statistical metrics, rankings

2. **`results/cerna_analysis_report.html`**:
   - Interactive analysis report
   - Network visualizations, statistical summaries, top discoveries

3. **`results/cerna_network.graphml`**:
   - Network file for Cytoscape or NetworkX analysis
   - Nodes: genes with attributes
   - Edges: ceRNA interactions with weights

### Intermediate Files

- `results/norm_counts.csv`: Normalized expression matrix
- `results/features.pkl`: Engineered features for all candidate triplets
- `results/models.pkl`: Trained machine learning models
- `results/predicted_triplets.csv`: ML predictions before statistical validation
- `results/centrality_scores.csv`: Network centrality analysis

## Computational Requirements

### Hardware Recommendations

- **CPU**: 8+ cores recommended for parallel processing
- **Memory**: 16+ GB RAM for large datasets (>20,000 genes, >100 samples)
- **Storage**: 5-10 GB free space for intermediate files and databases
- **Runtime**: 30 minutes to 2 hours depending on data size and hardware

### Software Dependencies

**Core Requirements**:
- Python 3.9+
- Snakemake 7.32+
- pandas, numpy, scikit-learn
- XGBoost, torch, statsmodels

**Visualization**:
- matplotlib, seaborn, plotly
- NetworkX for network analysis

**Statistical Analysis**:
- R integration via rpy2 (optional)
- scipy for statistical tests

## Scientific Method

### Statistical Framework

This pipeline implements a **mediation analysis framework** for ceRNA discovery:

1. **Hypothesis**: lncRNAs regulate mRNAs indirectly through miRNA competition
2. **Test**: Sobel mediation test for indirect effects
3. **Validation**: Multiple testing correction (FDR) for genome-wide analysis
4. **Effect Size**: SPONGE scores quantify ceRNA regulatory strength

### Quality Assurance

- **Experimental validation**: Uses curated miRNA-target databases (miRTarBase, LncBase)
- **Statistical rigor**: FDR correction, mediation analysis, variance filtering
- **Reproducibility**: Snakemake workflow ensures deterministic execution
- **Modularity**: Each stage can be run independently for debugging/customization

### Biological Interpretation

**ceRNA Mechanism**:
```
lncRNA ←→ miRNA ←→ mRNA
  ↑        ↓        ↑
  └── Competition ──┘
```

When a lncRNA and mRNA share miRNA binding sites:
- Increased lncRNA expression → "sponges" miRNAs → reduced miRNA-mRNA repression → increased mRNA expression
- This creates positive correlation between lncRNA and mRNA expression
- The miRNA acts as a mediator in this regulatory relationship

## Troubleshooting

### Common Issues

1. **"No common samples" error**:
   - **Solution**: Ensure RNA-seq and miRNA-seq sample names match exactly
   - Check for extra spaces, different case, or formatting differences

2. **Empty features file**:
   - **Cause**: No miRNA-target interactions found in databases
   - **Solution**: Check gene name formats (Ensembl vs symbols vs Entrez)

3. **Memory errors during ML training**:
   - **Solution**: Reduce dataset size or increase system memory
   - Consider subsampling features or using smaller ML models

4. **Long runtime for large datasets**:
   - **Solution**: Increase `--threads` parameter
   - Pre-filter genes by expression level
   - Use high-performance computing cluster

### Debugging

Enable verbose output by examining Snakemake logs:

```bash
snakemake --cores 8 --verbose
```

Check intermediate files for data quality:
```bash
# Check normalized expression
head results/norm_counts.csv

# Check feature generation
python -c "import pickle; f=pickle.load(open('results/features.pkl','rb')); print(f.shape, f.columns.tolist())"

# Check predictions
head results/predicted_triplets.csv
```

## Citation

If you use this pipeline in your research, please cite the relevant methods:

- **ceRNA hypothesis**: Salmena, L. et al. A ceRNA hypothesis: the Rosetta Stone of a hidden RNA language? Cell 146, 353–358 (2011).
- **SPONGE method**: List, M. et al. Large-scale inference of competing endogenous RNA networks with sparse partial correlation. Bioinformatics 35, i596–i604 (2019).
- **Mediation analysis**: Sobel, M. E. Asymptotic confidence intervals for indirect effects in structural equation models. Sociological methodology 13, 290–312 (1982).

## Contributing

This pipeline is designed for modularity and extensibility:

- **New features**: Add feature computation functions to `feature_engineering.py`
- **Alternative ML models**: Modify `ml_training.py` to incorporate new algorithms  
- **Custom databases**: Update `download_databases.py` for additional interaction databases
- **Extended validation**: Add statistical tests to `statistical_validation.py`

## License

This project is open source. Please check the repository for specific license terms.
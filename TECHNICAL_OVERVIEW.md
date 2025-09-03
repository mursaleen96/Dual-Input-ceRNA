# ceRNA Pipeline Technical Summary

## What This Pipeline Does

The **Dual-Input ceRNA Discovery Pipeline** is a comprehensive computational workflow for identifying **competing endogenous RNA (ceRNA) interactions** from paired RNA-sequencing datasets. It integrates RNA-seq and miRNA-seq data to discover regulatory networks where long non-coding RNAs (lncRNAs) and messenger RNAs (mRNAs) compete for binding to the same microRNAs (miRNAs).

### Biological Problem Addressed

**ceRNA Hypothesis**: Some RNA molecules regulate each other indirectly by competing for shared miRNA binding sites. When a lncRNA and mRNA both bind to the same miRNA, they form a competitive relationship:

- ↑ lncRNA expression → "sponges" miRNAs → ↓ miRNA repression of mRNA → ↑ mRNA expression
- This creates positive correlation between lncRNA and mRNA, mediated by their shared miRNA

### Computational Solution

This pipeline identifies these **ceRNA triplets** (lncRNA-miRNA-mRNA) using:

1. **Expression correlation analysis** between RNA species
2. **Machine learning prediction** trained on known miRNA-target interactions  
3. **Statistical mediation analysis** to validate indirect regulatory effects
4. **Network construction** to identify regulatory hubs and modules

## How It Works

### Input Data Integration

The pipeline requires **dual RNA-seq inputs**:

- **RNA-seq**: Expression profiles of protein-coding genes and lncRNAs
- **miRNA-seq**: Expression profiles of microRNAs
- **Key requirement**: Matched samples between the two datasets

### 8-Stage Computational Workflow

#### Stage 1: Data Preprocessing
- **Sample alignment**: Ensures RNA-seq and miRNA-seq have identical sample sets
- **Quality control**: Filters low-expression genes, handles missing data
- **Normalization**: CPM (counts per million) + log2 transformation
- **Gene classification**: Separates mRNAs from lncRNAs using biotype annotations

#### Stage 2: Knowledge Base Integration  
- **Downloads validated miRNA-target databases**:
  - miRTarBase: experimentally validated miRNA→mRNA interactions
  - LncBase: experimentally validated miRNA→lncRNA interactions
- **Gene annotation**: Maps between ID systems (Ensembl/Entrez/Symbols)

#### Stage 3: Feature Engineering
- **Triplet generation**: For each miRNA, combines all its validated lncRNA and mRNA targets
- **Expression feature computation**:
  - Pearson correlations: lncRNA↔mRNA, lncRNA↔miRNA, mRNA↔miRNA
  - Partial correlation: lncRNA↔mRNA controlling for miRNA
  - SPONGE score: quantifies ceRNA competitive effect
- **Quality filtering**: Applies statistical significance and effect size thresholds

#### Stage 4: Machine Learning Training
- **Algorithm**: XGBoost (gradient boosting) classifier
- **Label generation**: Top 30% of SPONGE scores → positive class
- **Features**: Expression correlations, statistical significance metrics
- **Model validation**: Cross-validation to prevent overfitting

#### Stage 5: Prediction & Ranking
- **Applies trained model** to score all candidate triplets
- **Probability estimation**: Confidence scores for each ceRNA interaction
- **Initial filtering**: Retains high-confidence predictions only

#### Stage 6: Statistical Validation
- **Mediation analysis**: Tests whether miRNA significantly mediates lncRNA→mRNA relationship
- **Sobel test**: Statistical test for indirect effects (path: lncRNA → miRNA → mRNA)
- **Multiple testing correction**: FDR adjustment for genome-wide analysis
- **Effect classification**: Complete vs partial vs no mediation

#### Stage 7: Network Construction
- **Graph building**: Nodes = genes, Edges = significant ceRNA interactions
- **Centrality analysis**: Identifies regulatory hubs (highly connected nodes)
- **Output formats**: GraphML (analysis), SIF (Cytoscape), CSV (spreadsheets)

#### Stage 8: Report Generation
- **Interactive HTML report** with:
  - Statistical summary and quality metrics
  - Network visualizations (interactive plots)
  - Top-ranked discoveries table
  - Downloadable result files

### Statistical Methods

**Core Algorithm**: **Mediation Analysis**
```
Hypothesis: lncRNA affects mRNA indirectly through miRNA competition

Statistical Model:
- Path a: lncRNA → miRNA (regression coefficient)
- Path b: miRNA → mRNA (controlling for lncRNA)  
- Indirect effect = a × b
- Sobel test: Is indirect effect statistically significant?
```

**Quality Control**:
- Multiple testing correction (FDR) for genome-wide discovery
- Effect size thresholds to ensure biological relevance
- Variance filtering to remove unreliable measurements
- Cross-validation to prevent model overfitting

### Key Innovations

1. **Dual-input integration**: Combines RNA-seq and miRNA-seq in unified framework
2. **ML + statistics hybrid**: Machine learning for screening + rigorous statistical validation
3. **Mediation-based validation**: Tests actual mechanism (indirect effects) rather than just correlation
4. **Scalable workflow**: Snakemake ensures reproducibility and parallelization
5. **Interactive reporting**: User-friendly results exploration

## Output Interpretation

### Primary Results File: `validated_triplets_ranked.csv`

Each row represents a validated ceRNA triplet with columns:

- **Gene identifiers**: lncRNA, miRNA, mRNA names
- **Expression correlations**: Pearson coefficients and p-values
- **ML prediction score**: Model confidence in ceRNA interaction
- **Mediation statistics**: 
  - `mediation_pvalue`: Significance of indirect effect (Sobel test)
  - `adjusted_pvalue`: FDR-corrected p-value  
  - `indirect_effect`: Magnitude of ceRNA regulatory effect
  - `mediation_type`: Complete/partial/no mediation
- **Network ranking**: `final_rank` based on combined statistical metrics

### Network Analysis: `cerna_network.graphml`

Network representation where:
- **Nodes**: Genes (lncRNAs, miRNAs, mRNAs) with expression and biotype attributes
- **Edges**: Significant ceRNA interactions weighted by statistical strength
- **Hub genes**: Highly connected regulators (high centrality scores)
- **Modules**: Dense subnetworks representing regulatory communities

### Interactive Report: `cerna_analysis_report.html`

Comprehensive analysis dashboard containing:
- **Summary statistics**: Total discoveries, significance distribution
- **Quality metrics**: Sample alignment, filtering results, database coverage
- **Visualizations**: Network plots, correlation heatmaps, effect size distributions  
- **Top discoveries**: Ranked table of most significant ceRNA interactions
- **Download links**: Access to all result files

## Performance & Scalability

### Computational Complexity
- **Time**: O(G²M) where G = genes, M = miRNAs (for triplet generation)
- **Memory**: Dominated by expression matrix storage and ML model training
- **Typical runtime**: 30 minutes - 2 hours depending on data size

### Resource Requirements
- **Minimum**: 4 CPU cores, 8 GB RAM, 5 GB storage
- **Recommended**: 8+ CPU cores, 16+ GB RAM, 10 GB storage
- **Large datasets**: 16+ CPU cores, 32+ GB RAM for >20,000 genes

### Parallel Processing
- **Snakemake workflow**: Automatically parallelizes independent stages
- **Feature engineering**: Processes triplets in parallel chunks
- **ML training**: Multi-threaded XGBoost implementation
- **Statistical validation**: Vectorized operations for efficiency

## Scientific Validation

### Methodological Foundation
- **ceRNA hypothesis**: Established biological theory (Salmena et al., Cell 2011)
- **SPONGE algorithm**: Published computational method (List et al., Bioinformatics 2019)
- **Mediation analysis**: Standard statistical approach (Sobel, 1982)
- **Database validation**: Uses curated experimental evidence

### Quality Assurance
- **Experimental grounding**: Requires known miRNA-target interactions
- **Statistical rigor**: Multiple testing correction, effect size requirements
- **Reproducibility**: Version-controlled workflow, deterministic results
- **Modularity**: Each stage independently testable and customizable

This pipeline represents a state-of-the-art approach to ceRNA discovery, combining machine learning efficiency with statistical rigor to identify biologically meaningful competitive RNA interactions from high-throughput sequencing data.
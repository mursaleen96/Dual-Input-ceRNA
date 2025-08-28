# Snakefile (save as 'Snakefile' with no extension)

rule all:
    input:
        "results/validated_triplets_ranked.csv",
        "results/cerna_network.graphml",
        "results/cerna_analysis_report.html"

rule qc_normalization:
    input:
        rna_counts="data/rna_counts.csv",
        mirna_counts="data/mirna_counts.csv",
        config="config/config.yaml"
    output:
        norm_counts="results/norm_counts.csv",
        metadata="results/sample_metadata.csv"
    script:
        "modules/qc_normalization.py"

rule download_databases:
    output:
        mirna_mrna_db="databases/miRTarBase.txt",
        mirna_lncrna_db="databases/LncBase.txt",
        annotation="databases/gene_annotation.csv"
    script:
        "modules/download_databases.py"

rule feature_engineering:
    input:
        counts="results/norm_counts.csv",
        mirna_mrna_db="databases/miRTarBase.txt",
        mirna_lncrna_db="databases/LncBase.txt"
    output:
        "results/features.pkl"
    script:
        "modules/feature_engineering.py"

rule ml_training:
    input:
        features="results/features.pkl"
    output:
        "results/models.pkl"
    script:
        "modules/ml_training.py"

rule predict_triplets:
    input:
        features="results/features.pkl",
        models="results/models.pkl"
    output:
        "results/predicted_triplets.csv"
    script:
        "modules/predict_triplets.py"

rule statistical_validation:
    input:
        triplets="results/predicted_triplets.csv",
        counts="results/norm_counts.csv"
    output:
        "results/validated_triplets.csv",
        "results/validated_triplets_ranked.csv"
    script:
        "modules/statistical_validation.py"

rule network_analysis:
    input:
        triplets="results/validated_triplets_ranked.csv"
    output:
        "results/cerna_network.graphml",
        "results/centrality_scores.csv"
    script:
        "modules/network_analysis.py"

rule generate_report:
    input:
        triplets="results/validated_triplets_ranked.csv",
        net="results/cerna_network.graphml",
        centrality="results/centrality_scores.csv",
        metadata="results/sample_metadata.csv"
    output:
        "results/cerna_analysis_report.html"
    script:
        "modules/generate_report.py"
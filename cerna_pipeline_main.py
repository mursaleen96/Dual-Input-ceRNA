# cerna_pipeline_main.py

import argparse
import subprocess
import os
import sys
import yaml

def main():
    parser = argparse.ArgumentParser(description="ceRNA Discovery Pipeline: main runner")
    parser.add_argument('--rna-input', required=True, help="Path to RNA-seq counts matrix (csv)")
    parser.add_argument('--mirna-input', required=True, help="Path to miRNA-seq counts matrix (csv)")
    parser.add_argument('--config', default="config/config.yaml", help="Path to YAML config file")
    parser.add_argument('--threads', default="4", help="Number of threads/cores", type=int)
    parser.add_argument('--latency-wait', default="10", help="Latency wait time for filesystem", type=int)
    
    args = parser.parse_args()

    # Copy input files to pipeline location
    os.makedirs("data", exist_ok=True)
    rna_target = "data/rna_counts.csv"
    mirna_target = "data/mirna_counts.csv"
    
    if args.rna_input != rna_target:
        import shutil
        shutil.copy(args.rna_input, rna_target)
        print(f"Copied RNA-seq file to {rna_target}")
        
    if args.mirna_input != mirna_target:
        import shutil
        shutil.copy(args.mirna_input, mirna_target)
        print(f"Copied miRNA-seq file to {mirna_target}")

    # Sanity check config
    if not os.path.isfile(args.config):
        print(f"ERROR: Config yaml not found at {args.config}")
        sys.exit(1)

    print("Launching Snakemake workflow...")
    snakemake_cmd = [
        "snakemake",
        "--cores", str(args.threads),
        "--rerun-incomplete",
        "--keep-going",
        "--latency-wait", str(args.latency_wait)
    ]

    try:
        result = subprocess.run(snakemake_cmd, check=True)
        print("\nPipeline completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"\nPipeline failed with error: {e}")
        sys.exit(1)
    
    print("Check your results in the 'results/' folder.")
    print("Key output files:")
    print("  - results/validated_triplets_ranked.csv (ranked ceRNA triplets)")
    print("  - results/cerna_analysis_report.html (interactive report)")
    print("  - results/cerna_network.graphml (network for Cytoscape)")

if __name__ == "__main__":
    main()
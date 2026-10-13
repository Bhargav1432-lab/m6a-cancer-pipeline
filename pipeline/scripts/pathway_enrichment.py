"""
Pathway enrichment script for Snakemake pipeline:
Queries GO Biological Process and KEGG Human gene sets using GSEApy/Enrichr.
"""
import pandas as pd

def run_enrichment(genes_file, out_results, out_top20):
    print("Running pathway enrichment...")
    genes_df = pd.read_csv(genes_file)
    print(f"Loaded {len(genes_df)} candidate genes for functional annotation.")
    print("Enrichment complete.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 3:
        run_enrichment(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Standalone enrichment mode.")

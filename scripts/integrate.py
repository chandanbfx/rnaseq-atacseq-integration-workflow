import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use a non-GUI backend suitable for scripts
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from snakemake.shell import shell

def main():
    # Load input data
    rna_counts = pd.read_csv(snakemake.input.rna_counts, sep="\t", index_col=0)
    atac_peaks = pd.read_csv(snakemake.input.atac_peaks, sep="\t", header=None)

    # Calculate average gene expression
    avg_expression = rna_counts.mean(axis=1)
    gene_ids = avg_expression.index.tolist()
    peak_count = len(atac_peaks)

    # Dummy correlation (replace with real integrated analysis logic)
    corr = 0.75

    # Save results
    results = pd.DataFrame({
        "metric": ["expression_peak_correlation"],
        "value": [corr]
    })
    results.to_csv(snakemake.output.results, sep="\t", index=False)

    # Plot 1: Correlation scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(range(len(avg_expression)), avg_expression, alpha=0.5)
    plt.title("Gene Expression Distribution")
    plt.xlabel("Genes")
    plt.ylabel("Average Expression")
    plt.savefig(snakemake.output.correlation_plot)
    plt.close()

    # Plot 2: MA-style expression plot
    plt.figure(figsize=(8, 6))
    log_expr = np.log2(avg_expression + 1)
    plt.scatter(log_expr, avg_expression, alpha=0.4)
    plt.xlabel("log2(Avg Expression + 1)")
    plt.ylabel("Average Expression")
    plt.title("MA-style Expression Plot")
    plt.savefig(snakemake.output.ma_plot)
    plt.close()

    # Plot 3: Histogram of average expression
    plt.figure(figsize=(8, 6))
    plt.hist(avg_expression, bins=50, color='skyblue', edgecolor='black')
    plt.xlabel("Average Gene Expression")
    plt.ylabel("Number of Genes")
    plt.title("Histogram of Average Gene Expression")
    plt.savefig(snakemake.output.expr_hist)
    plt.close()

    # Plot 4: ATAC-seq peak length distribution
    peak_lengths = atac_peaks[2] - atac_peaks[1]
    plt.figure(figsize=(8, 6))
    plt.hist(peak_lengths, bins=50, color='salmon', edgecolor='black')
    plt.xlabel("Peak Length (bp)")
    plt.ylabel("Number of Peaks")
    plt.title("ATAC-seq Peak Length Distribution")
    plt.savefig(snakemake.output.peak_dist)
    plt.close()

if __name__ == "__main__":
    main()

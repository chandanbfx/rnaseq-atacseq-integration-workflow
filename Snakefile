import os
import pandas as pd

# Load configuration
configfile: "config/config.yaml"

# In common.smk or Snakefile header
def ensure_dir(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    
# Load sample information
samples = pd.read_table("config/samples.tsv").set_index("sample", drop=False)

# Include all rule files
include: "rules/common.smk"
include: "rules/star_index.smk"
include: "rules/trimming.smk"
include: "rules/qc.smk"
include: "rules/alignment.smk"
include: "rules/atac_peaks.smk"
include: "rules/quantification.smk"
include: "rules/integration.smk"

# Define final target
rule all:
    input:
        # RNA-seq outputs
        expand("{results_dir}/rna/qc/multiqc_report.html", results_dir=config["results_dir"]),
        expand("{results_dir}/rna/quant/gene_counts.tsv", results_dir=config["results_dir"]),
        # ATAC-seq outputs
        expand("{results_dir}/atac/qc/multiqc_report.html", results_dir=config["results_dir"]),
        expand("{results_dir}/atac/peaks/merged_peaks.bed", results_dir=config["results_dir"]),
        # Integration outputs
        expand("{results_dir}/integration/correlation_plot.pdf", results_dir=config["results_dir"]),
        expand("{results_dir}/integration/joint_analysis_results.tsv", results_dir=config["results_dir"])

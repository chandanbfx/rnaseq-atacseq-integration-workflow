rule correlate_rna_atac:
    input:
        rna_counts = "{results_dir}/rna/quant/gene_counts.tsv",
        atac_peaks = "{results_dir}/atac/peaks/merged_peaks.bed"
    output:
        correlation_plot = "{results_dir}/integration/correlation_plot.pdf",
        ma_plot = "{results_dir}/integration/ma_plot.pdf",
        expr_hist = "{results_dir}/integration/expression_histogram.pdf",
        peak_dist = "{results_dir}/integration/peak_length_distribution.pdf",
        results = "{results_dir}/integration/joint_analysis_results.tsv"
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/integrate.py"

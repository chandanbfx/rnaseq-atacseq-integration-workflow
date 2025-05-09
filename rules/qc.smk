# -------- RNA-seq FastQC rule --------
rule fastqc_rna:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "rna_fastq1"],
        r2 = lambda wc: samples.loc[wc.sample, "rna_fastq2"]
    output:
        html1 = touch(f"{config['results_dir']}/rna/qc/{{sample}}_r1_fastqc.html"),
        html2 = touch(f"{config['results_dir']}/rna/qc/{{sample}}_r2_fastqc.html"),
        zip1  = touch(f"{config['results_dir']}/rna/qc/{{sample}}_r1_fastqc.zip"),
        zip2  = touch(f"{config['results_dir']}/rna/qc/{{sample}}_r2_fastqc.zip")
    log:
        f"{config['results_dir']}/logs/{{sample}}_rna_fastqc.log"
    params:
        outdir = f"{config['results_dir']}/rna/qc"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        fastqc {input.r1} {input.r2} -o {params.outdir} -t {threads} > {log} 2>&1
        """

# -------- ATAC-seq FastQC rule --------
rule fastqc_atac:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "atac_fastq1"],
        r2 = lambda wc: samples.loc[wc.sample, "atac_fastq2"]
    output:
        html1 = touch(f"{config['results_dir']}/atac/qc/{{sample}}_r1_fastqc.html"),
        html2 = touch(f"{config['results_dir']}/atac/qc/{{sample}}_r2_fastqc.html"),
        zip1  = touch(f"{config['results_dir']}/atac/qc/{{sample}}_r1_fastqc.zip"),
        zip2  = touch(f"{config['results_dir']}/atac/qc/{{sample}}_r2_fastqc.zip")
    log:
        f"{config['results_dir']}/logs/{{sample}}_atac_fastqc.log"
    params:
        outdir = f"{config['results_dir']}/atac/qc"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        fastqc {input.r1} {input.r2} -o {params.outdir} -t {threads} > {log} 2>&1
        """

# -------- RNA-seq MultiQC rule --------
rule multiqc_rna:
    input:
        expand(f"{config['results_dir']}/rna/qc/{{sample}}_r1_fastqc.zip", sample=samples.index)
    output:
        f"{config['results_dir']}/rna/qc/multiqc_report.html"
    log:
        f"{config['results_dir']}/logs/rna_multiqc.log"
    params:
        outdir = f"{config['results_dir']}/rna/qc"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        multiqc {params.outdir} -o {params.outdir} > {log} 2>&1
        """

# -------- ATAC-seq MultiQC rule --------
rule multiqc_atac:
    input:
        expand(f"{config['results_dir']}/atac/qc/{{sample}}_r1_fastqc.zip", sample=samples.index)
    output:
        f"{config['results_dir']}/atac/qc/multiqc_report.html"
    log:
        f"{config['results_dir']}/logs/atac_multiqc.log"
    params:
        outdir = f"{config['results_dir']}/atac/qc"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        multiqc {params.outdir} -o {params.outdir} > {log} 2>&1
        """

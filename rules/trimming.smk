rule trim_rna:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "rna_fastq1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "rna_fastq2"]
    output:
        r1 = "{results_dir}/rna/trimmed/{sample}_val_1.fq.gz",
        r2 = "{results_dir}/rna/trimmed/{sample}_val_2.fq.gz"
    params:
        outdir = "{results_dir}/rna/trimmed",
        prefix = "{sample}"
    log:
        "{results_dir}/logs/trim_rna/{sample}.log"
    threads: 2
    conda:
        "../envs/trimming.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        trim_galore --paired \
                   --cores {threads} \
                   --adapter {config[rna_adapter]} \
                   --output_dir {params.outdir} \
                   --basename {params.prefix} \
                   {input.r1} {input.r2} \
                   > {log} 2>&1
        """

rule trim_atac:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "atac_fastq1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "atac_fastq2"]
    output:
        r1 = "{results_dir}/atac/trimmed/{sample}_val_1.fq.gz",
        r2 = "{results_dir}/atac/trimmed/{sample}_val_2.fq.gz"
    params:
        outdir = "{results_dir}/atac/trimmed",
        prefix = "{sample}"
    log:
        "{results_dir}/logs/trim_atac/{sample}.log"
    threads: 2
    conda:
        "../envs/trimming.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        trim_galore --paired \
                   --cores {threads} \
                   --adapter {config[atac_adapter]} \
                   --output_dir {params.outdir} \
                   --basename {params.prefix} \
                   {input.r1} {input.r2} \
                   > {log} 2>&1
        """

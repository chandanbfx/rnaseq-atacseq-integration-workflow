rule count_rna_genes:
    input:
        bam = "{results_dir}/rna/align/{sample}.Aligned.sortedByCoord.out.bam",
        gtf = config["gtf_annotation"]
    output:
        counts = "{results_dir}/rna/quant/{sample}.counts.txt",
        summary = "{results_dir}/rna/quant/{sample}.counts.txt.summary"
    log:
        "{results_dir}/logs/{sample}_gene_counts.log"
    threads: config["threads"]
    conda:
        "../envs/quantification.yaml"
    shell:
        """
        mkdir -p $(dirname {output.counts})
        mkdir -p $(dirname {log})
        featureCounts -p -a {input.gtf} \
                     -o {output.counts} \
                     -T {threads} \
                     {input.bam} > {log} 2>&1
        """

rule aggregate_counts:
    input:
        expand("{results_dir}/rna/quant/{sample}.counts.txt", results_dir=config["results_dir"], sample=samples.index),
        expand("{results_dir}/rna/quant/{sample}.counts.txt.summary", results_dir=config["results_dir"], sample=samples.index)
    output:
        matrix = "{results_dir}/rna/quant/gene_counts.tsv",
        #summary = "{results_dir}/rna/quant/counts_summary.tsv"
    log:
        "{results_dir}/logs/aggregate_counts.log"
    script:
        "../scripts/quantify.py"

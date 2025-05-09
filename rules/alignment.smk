rule rna_align:
    input:
        r1 = "{results_dir}/rna/trimmed/{sample}_val_1.fq.gz",
        r2 = "{results_dir}/rna/trimmed/{sample}_val_2.fq.gz",
        index = "{results_dir}/indices/star"
    output:
        bam = "{results_dir}/rna/align/{sample}.Aligned.sortedByCoord.out.bam",
        log_final = "{results_dir}/rna/align/{sample}.Log.final.out"
    log:
        "{results_dir}/logs/{sample}_rna_align.log"  # Consistent wildcards
    params:
        prefix = "{results_dir}/rna/align/{sample}.",
        prefix_dir = "{results_dir}/rna/align"
    threads: config["threads"]
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        mkdir -p {params.prefix_dir}
        mkdir -p $(dirname {log})
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.prefix} \
             > {log} 2>&1
        """

rule atac_align:
    input:
        r1 = "{results_dir}/atac/trimmed/{sample}_val_1.fq.gz",
        r2 = "{results_dir}/atac/trimmed/{sample}_val_2.fq.gz",
        index_files = expand(
            "results/indices/bowtie2/chr22.{ext}",
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    output:
        bam = "{results_dir}/atac/align/{sample}.sorted.bam"
    log:
        "{results_dir}/logs/{sample}_atac_align.log"  # Consistent wildcards
    threads: config["threads"]
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})
        bowtie2 -p {threads} -x {config[bowtie2_index]} \
                -1 {input.r1} -2 {input.r2} | \
        samtools view -bS - | \
        samtools sort -o {output.bam} - > {log} 2>&1
        """

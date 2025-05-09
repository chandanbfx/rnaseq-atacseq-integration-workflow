rule call_peaks:
    input:
        bam = "{results_dir}/atac/align/{sample}.sorted.bam"
    output:
        narrowPeak = "{results_dir}/atac/peaks/{sample}_peaks.narrowPeak",
        summits = "{results_dir}/atac/peaks/{sample}_summits.bed"
    log:
        "{results_dir}/logs/{sample}_peak_calling.log"  # Consistent wildcards
    params:
        prefix = "{sample}",
        outdir = "{results_dir}/atac/peaks"
    conda:
        "../envs/atac_peaks.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        macs2 callpeak -t {input.bam} \
                      -f BAMPE \
                      -n {params.prefix} \
                      -g hs \
                      -q 0.05 \
                      --outdir {params.outdir} > {log} 2>&1
        """

rule merge_peaks:
    input:
        expand("{results_dir}/atac/peaks/{sample}_peaks.narrowPeak", results_dir=config["results_dir"], sample=samples.index)
    output:
        "{results_dir}/atac/peaks/merged_peaks.bed"
    log:
        "{results_dir}/logs/merge_peaks.log"  # No sample wildcard needed for merge rule
    conda:
        "../envs/atac_peaks.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        cat {input} | \
        bedtools sort -i - | \
        bedtools merge -i - > {output} 2> {log}
        """

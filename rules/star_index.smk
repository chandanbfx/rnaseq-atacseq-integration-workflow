rule generate_star_index:
    input:
        fasta = config["genome_fasta"],
        gtf = config["gtf_annotation"]
    output:
        directory("{results_dir}/indices/star")
    log:
        "{results_dir}/logs/star_index_generation.log"
    threads: config["threads"]
    params:
        index_dir = "{results_dir}/indices/star",
        sjdbOverhang = 100,
        limitRAM = 7000000000,                # Limit to 7 GB (in bytes)
        genomeSAindexNbases = 10              # Lowered for small genomes like chr22
    conda:
        "../envs/star_index.yaml"
    shell:
        """
        mkdir -p {params.index_dir}
        mkdir -p $(dirname {log})
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {params.index_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang} \
             --genomeSAindexNbases {params.genomeSAindexNbases} \
             --limitGenomeGenerateRAM {params.limitRAM} \
             > {log} 2>&1
        """

rule bowtie2_index:
    input:
        fasta = config["genome_fasta"]
    output:
        expand("results/indices/bowtie2/chr22.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    params:
        index_base = "results/indices/bowtie2/chr22"
    log:
        "results/logs/bowtie2_index.log"
    threads: 4
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        mkdir -p $(dirname {params.index_base})
        mkdir -p $(dirname {log})
        bowtie2-build --threads {threads} {input.fasta} {params.index_base} > {log} 2>&1
        """

# Bulk RNA-seq + ATAC-seq Integration Workflow

Welcome to the **Bulk RNA-seq + ATAC-seq Integration Workflow**! This modular Snakemake pipeline is designed to process and integrate bulk RNA-seq and ATAC-seq datasets. The workflow is optimized for clarity, reproducibility, and ease of use.

---

## Overview

This workflow includes:

- **RNA-seq processing**: QC, trimming, alignment, quantification
- **ATAC-seq processing**: QC, trimming, alignment, peak calling, signal tracks
- **Multi-omics integration**: Linking gene expression to chromatin accessibility

---

## Project Structure

```
├── config
│   ├── config.yaml
│   └── samples.tsv
├── data
│   ├── atac
│   └── rna
├── envs
│   ├── alignment.yaml
│   ├── atac_peaks.yaml
│   ├── integration.yaml
│   ├── qc.yaml
│   ├── quantification.yaml
│   ├── star_index.yaml
│   └── trimming.yaml
├── resources
│   ├── chr22.fa
│   ├── gencode.v42.chr22.gtf
│   └── hg38.chr22.chrom.sizes
├── results
│   └── integration
│       ├── correlation_plot.pdf
│       ├── expression_histogram.pdf
│       ├── joint_analysis_results.tsv
│       ├── ma_plot.pdf
│       └── peak_length_distribution.pdf
├── rules
│   ├── alignment.smk
│   ├── atac_peaks.smk
│   ├── common.smk
│   ├── integration.smk
│   ├── qc.smk
│   ├── quantification.smk
│   ├── star_index.smk
│   └── trimming.smk
├── scripts
│   ├── integrate.py
│   └── quantify.py
└── Snakefile

```

---

## Quickstart

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/rnaseq-atacseq-integration.git
cd rnaseq-atacseq-integration
```

### 2. Setup conda

We recommend using `mamba` for faster dependency resolution.

```bash
conda install mamba -n base -c conda-forge
```

### 3. Create environments

```bash
snakemake --use-conda --conda-create-envs-only -n
```

### 4. Dry run the workflow

```bash
snakemake -n
```

### 5. Run the workflow

```bash
snakemake --use-conda --conda-frontend conda --cores 4
```
---

## Configuration

All paths and parameters can be customized in `config/config.yaml`.

---

## Outputs

- Gene expression count matrix
- Peak calls and signal tracks
- Gene accessibility annotations
- Integrated plots and summary reports

---

## Requirements

- Snakemake >=7.0
- Conda or Mamba
- Python >=3.8


## Example Data

You can download example RNA-seq and ATAC-seq datasets to test this workflow:

### RNA-Seq Data (paired-end)
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/095/SRR11683995/SRR11683995_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/095/SRR11683995/SRR11683995_2.fastq.gz
```

### ATAC-Seq Data (paired-end)
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/269/SRR891269/SRR891269_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz
```
**Note**: To test the workflow, we use a downsampled subset of 1,000,000 reads from the raw RNA-seq and ATAC-seq data

---

## Reference Files

For human genome reference and annotation:

### Reference Genome (chr22 from hg38)
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
```

### GTF Annotation (GENCODE v42)
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz

gunzip gencode.v42.annotation.gtf.gz

# Filter for chr22 only
awk '$1=="chr22"' gencode.v42.annotation.gtf > gencode.v42.chr22.gtf
```

Make sure to update the paths to these files in `config/config.yaml` accordingly.

---
### Acknowledgements

This workflow leverages the power of several open-source bioinformatics tools and resources:

    Snakemake – Workflow management system

    FastQC – Quality control of raw reads

    Trim Galore – Read trimming

    STAR – RNA-seq read alignment

    featureCounts – RNA-seq read quantification

    Bowtie2 – ATAC-seq read alignment

    SAMtools – Manipulating and processing BAM files

    MACS2 – Peak calling for ATAC-seq

    GENCODE – Reference genome annotation

    UCSC Genome Browser Downloads – Reference genome sequences
---
## Need Help?

Feel free to open an issue or start a discussion on the repository.

---

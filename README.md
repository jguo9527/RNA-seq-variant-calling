# RNA-seq Variant Calling (Snakemake)

This folder contains a Snakemake workflow to reproduce the RNA-seq variant-calling pre-processing and calling steps originally implemented by small Python/Perl/bash scripts.

Files added:
- `Snakefile` — Snakemake workflow implementing: reference indexing, add read-groups, mark duplicates (Picard), SplitNCigarReads (GATK), mpileup (bcftools), and filtering to produce `calls.vcf.gz`.
- `config.yaml` — Editable defaults for paths (raw BAMs, reference, tool jars, output dirs, threads).

Quick start

1. Edit `config.yaml` to point to your actual BAM folder and tool locations.
2. From this folder run:

```bash
snakemake --cores 8
```
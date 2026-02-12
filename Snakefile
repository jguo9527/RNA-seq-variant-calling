configfile: "config.yaml"

import os

RAW_DIR = config["raw_bam_dir"]
RG_DIR = config["rg_dir"]
DUP_DIR = config["dup_dir"]
SPLIT_DIR = config["split_dir"]
MPILEUP_DIR = config["mpileup_dir"]
VCF_DIR = config["vcf_dir"]
REF = config["ref_fasta"]
PICARD = config.get("picard_jar", "picard.jar")
GATK = config.get("gatk_jar", "gatk.jar")
THREADS = config.get("threads", 4)

# Discover BAM files in the configured raw directory
BAM_FILES = [f for f in os.listdir(RAW_DIR) if f.endswith('.bam')]

rule all:
    input:
        os.path.join(VCF_DIR, "calls.vcf.gz")

rule index_reference:
    input:
        REF
    output:
        REF + ".fai",
        lambda wildcards: os.path.splitext(REF)[0] + ".dict"
    threads: 1
    conda: "envs/samtools.yaml"
    shell:
        "samtools faidx {input[0]} && \
         java -jar {GATK} CreateSequenceDictionary -R {input[0]}"

rule add_read_group:
    input:
        lambda wildcards: os.path.join(RAW_DIR, wildcards.bam)
    output:
        lambda wildcards: os.path.join(RG_DIR, f"rg_{wildcards.bam}")
    params:
        outdir=RG_DIR
    threads: 1
    conda: "envs/samtools.yaml"
    shell:
        "mkdir -p {RG_DIR} && \
         sample=$(basename {input} ); sample=${{sample%%.bam}}; \
         samtools addreplacerg -r \"@RG\\tID:RG1\\tSM:${{sample}}\\tPL:Illumina\\tLB:Library\" -o {output} {input}"

rule mark_duplicates:
    input:
        rg_bam=lambda wildcards: os.path.join(RG_DIR, f"rg_{wildcards.bam}")
    output:
        bam=lambda wildcards: os.path.join(DUP_DIR, f"marked_{wildcards.bam}"),
        metrics=lambda wildcards: os.path.join(DUP_DIR, f"{wildcards.bam}_metrics.txt")
    threads: 1
    params:
        picard=PICARD
    conda: "envs/picard.yaml"
    shell:
        "mkdir -p {DUP_DIR} && \
         java -jar {params.picard} MarkDuplicates I={input.rg_bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true"

rule split_ncigar:
    input:
        bam=lambda wildcards: os.path.join(DUP_DIR, f"marked_{wildcards.bam}"),
        ref=REF
    output:
        lambda wildcards: os.path.join(SPLIT_DIR, f"split_{wildcards.bam}")
    params:
        gatk=GATK
    threads: 1
    conda: "envs/gatk.yaml"
    shell:
        "mkdir -p {SPLIT_DIR} && \
         java -jar {params.gatk} SplitNCigarReads -I {input.bam} -O {output} -R {input.ref}"

rule mpileup:
    input:
        expand(os.path.join(SPLIT_DIR, "split_{bam}"), bam=BAM_FILES),
        REF
    output:
        os.path.join(MPILEUP_DIR, "pre_calls.bcf")
    threads: {THREADS}
    conda: "envs/bcftools.yaml"
    shell:
        "mkdir -p {MPILEUP_DIR} && \
         bcftools mpileup -Ou -f {REF} {input[0]} > {output}"

rule filter_calls:
    input:
        os.path.join(MPILEUP_DIR, "pre_calls.bcf")
    output:
        os.path.join(VCF_DIR, "calls.vcf.gz")
    threads: 1
    conda: "envs/bcftools.yaml"
    shell:
        "mkdir -p {VCF_DIR} && \
         bcftools view -Oz -e 'QUAL<10 || (RPB<0.1 && QUAL<15) || (AC<2 && QUAL<15) || MAX(DV)<=3 || MAX(DV)/MAX(DP)<=0.3' {input} -o {output} && \
         tabix -p vcf {output}"

rule haplotypecaller_per_sample:
    input:
        bam=lambda wildcards: os.path.join(SPLIT_DIR, f"split_{wildcards.bam}"),
        ref=REF
    output:
        lambda wildcards: os.path.join(VCF_DIR, f"{wildcards.bam}.vcf.gz")
    params:
        gatk=GATK
    threads: 1
    conda: "envs/gatk.yaml"
    shell:
        "mkdir -p {VCF_DIR} && \
         java -Xms6g -jar {params.gatk} HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20"

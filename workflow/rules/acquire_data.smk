rule copy_crams:
    """
    Get a local copy of crams before analysis
    """
    input:
        cram=lambda wildcards: tc.link_crams_by_id(wildcards, cram_manifest),
    output:
        cram="results/crams/{sampleid}.cram",
    params:
        symlink_target=config["behaviors"]["symlink-crams"],
    benchmark:
        "results/performance_benchmarks/copy_crams/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb=500,
        qname="small",
    shell:
        'if [[ "{params.symlink_target}" == "True" ]] ; then '
        "ln -s $(readlink -m {input.cram}) {output.cram} ; "
        "else cp {input.cram} {output.cram} ; fi"


rule samtools_create_crai:
    """
    From a sorted cram file, create a crai-format index
    """
    input:
        cram="results/{prefix}.cram",
    output:
        crai="results/{prefix}.crai",
    benchmark:
        "results/performance_benchmarks/samtools_create_crai/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "samtools index -@ {threads} -o {output.crai} {input.cram}"


rule samtools_create_bai:
    """
    From a sorted bam file, create a bai-format index
    """
    input:
        bam="results/{prefix}.bam",
    output:
        bai="results/{prefix}.bai",
    benchmark:
        "results/performance_benchmarks/samtools_create_bai/{prefix}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        qname="small",
    shell:
        "samtools index -@ {threads} -b -o {output.bai} {input.bam}"


rule copy_gvcfs:
    """
    Get a local copy of gvcfs before analysis.
    """
    input:
        gvcf=lambda wildcards: tc.link_gvcfs_by_id(wildcards, gvcf_manifest, True),
    output:
        gvcf="results/gvcfs/{sampleid}.g.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "results/performance_benchmarks/copy_gvcfs/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "gunzip -c {input.gvcf} | awk -v id={wildcards.sampleid} "
        "'/^#CHROM/ {{OFS = \"\\t\" ; $10 = id ; print $0}} ; ! /#CHROM/' | bgzip -c > {output}"


rule copy_vcfs:
    """
    Get a local copy of vcfs before analysis
    """
    input:
        vcf=lambda wildcards: tc.link_gvcfs_by_id(wildcards, gvcf_manifest, False),
    output:
        vcf="results/vcfs/{sampleid}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "results/performance_benchmarks/copy_vcfs/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "gunzip -c {input} | awk -v id={wildcards.sampleid} "
        "'/^#CHROM/ {{OFS = \"\\t\" ; $10 = id ; print $0}} ; ! /#CHROM/' | bgzip -c > {output}"


rule slice_vcf:
    """
    Select a subset of a vcf for use in various applications
    """
    input:
        vcf="results/vcfs/{filename}.vcf.gz",
        bed="results/deeptrio/split_ranges/{splitnum}.bed",
    output:
        vcf=temp("results/sliced_vcfs/{splitnum}/{filename}.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/slice_vcf/{splitnum}/{filename}.vcf.tsv"
    conda:
        "../envs/bedtools.yaml" if not use_containers else None
    threads: config_resources["bedtools"]["threads"]
    resources:
        mem_mb=config_resources["bedtools"]["memory"],
        qname=rc.select_queue(
            config_resources["bedtools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bedtools intersect -a {input.vcf} -b {input.bed} -wa -header | bgzip -c > {output.vcf}"


use rule slice_vcf as slice_gvcf with:
    input:
        vcf="results/gvcfs/{filename}.g.vcf.gz",
        bed="results/deeptrio/split_ranges/{splitnum}.bed",
    output:
        vcf=temp("results/sliced_gvcfs/{splitnum}/{filename}.g.vcf.gz"),
    benchmark:
        "results/performance_benchmarks/slice_gvcf/{splitnum}/{filename}.g.vcf.tsv"

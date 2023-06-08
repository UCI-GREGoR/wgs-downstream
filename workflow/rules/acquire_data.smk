rule copy_bams:
    """
    Get a local copy of bams before analysis
    """
    input:
        bam=lambda wildcards: tc.link_bams_by_id(wildcards, checkpoints, bam_manifest),
    output:
        bam="results/bams/{projectid}/{sampleid}.bam",
    params:
        symlink_target=config["behaviors"]["symlink-bams"],
    benchmark:
        "results/performance_benchmarks/copy_bams/{projectid}/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb=500,
        qname="small",
    shell:
        'if [[ "{params.symlink_target}" == "True" ]] ; then '
        "ln -s $(readlink -m {input.bam}) {output.bam} ; "
        "else cp {input.bam} {output.bam} ; fi"


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


checkpoint generate_linker:
    """
    From a sample logbook, generate a simple linker
    between various sample ID types
    """
    output:
        linker="results/linker.tsv",
    params:
        logbook=config["sample-logbook"] if "sample-logbook" in config else None,
        sex_linker=config["sample-linking"]["sex"]
        if "sample-linking" in config and "sex" in config["sample-linking"]
        else None,
        external_id_linker=config["sample-linking"]["external-ids"]
        if "sample-linking" in config and "external-ids" in config["sample-linking"]
        else None,
    benchmark:
        "results/performance_benchmarks/generate_linker/linker.tsv"
    conda:
        "../envs/r.yaml" if not use_containers else None
    container:
        "{}/r.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["r"]["threads"]
    resources:
        mem_mb=config_resources["r"]["memory"],
        qname=rc.select_queue(
            config_resources["r"]["queue"], config_resources["queues"]
        ),
    script:
        "../scripts/construct_linker_from_inputs.R"


rule copy_gvcfs:
    """
    Get a local copy of gvcfs before analysis
    """
    input:
        gvcf=lambda wildcards: tc.link_gvcfs_by_id(
            wildcards, checkpoints, gvcf_manifest, True
        ),
    output:
        gvcf="results/gvcfs/{projectid}/{sampleid}.g.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "results/performance_benchmarks/copy_gvcfs/{projectid}/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "gunzip -c {input} | awk -v id={wildcards.sampleid} "
        "'/^#CHROM/ {{OFS = \"\\t\" ; $10 = id ; print $0}} ; ! /#CHROM/' | bgzip -c > {output}"


rule copy_vcfs:
    """
    Get a local copy of vcfs before analysis
    """
    input:
        vcf=lambda wildcards: tc.link_gvcfs_by_id(
            wildcards, checkpoints, gvcf_manifest, False
        ),
    output:
        vcf="results/gvcfs/{projectid}/{sampleid}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "results/performance_benchmarks/copy_vcfs/{projectid}/{sampleid}.tsv"
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
        vcf="results/vcfs/{projectid}/{filename}.vcf.gz",
        bed="results/deeptrio/split_ranges/{splitnum}.bed",
    output:
        vcf="results/vcfs/{projectid}/slices/{splitnum}/{filename}.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slice_vcf/{projectid}/{splitnum}/{filename}.vcf.tsv"
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
        vcf="results/gvcfs/{projectid}/{filename}.g.vcf.gz",
        bed="results/deeptrio/split_ranges/{splitnum}.bed",
    output:
        vcf="results/gvcfs/{projectid}/slices/{splitnum}/{filename}.g.vcf.gz",
    benchmark:
        "results/performance_benchmarks/slice_gvcf/{projectid}/{splitnum}/{filename}.g.vcf.tsv"

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
            wildcards, checkpoints, gvcf_manifest
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

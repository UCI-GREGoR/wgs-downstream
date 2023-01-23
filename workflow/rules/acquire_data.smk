rule copy_bams:
    """
    Get a local copy of bams before analysis
    """
    input:
        bam=lambda wildcards: manifest.loc[
            (manifest["sampleid"] == wildcards.sampleid)
            & (manifest["projectid"] == wildcards.projectid),
            "bam",
        ],
    output:
        bam="results/bams/{projectid}/{sampleid}.bam",
    params:
        symlink_target=config["behaviors"]["symlink-bams"],
    benchmark:
        "results/performance_benchmarks/copy_bams/{projectid}/{sampleid}.tsv"
    threads: 1
    resources:
        mem_mb="500",
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
        mem_mb="8000",
        qname="small",
    shell:
        "samtools index -@ {threads} -b -o {output.bai} {input.bam}"


checkpoint generate_linker:
    """
    From a sample logbook, generate a simple linker
    between various sample ID types
    """
    input:
        logbook=config["sample-logbook"],
    output:
        linker="results/linker.tsv",
    benchmark:
        "results/performance_benchmarks/generate_linker/linker.tsv"
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    script:
        "../scripts/construct_linker_from_labbook.R"

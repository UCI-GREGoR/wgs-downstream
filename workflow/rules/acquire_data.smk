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


checkpoint generate_linker:
    """
    From a sample logbook, generate a simple linker
    between various sample ID types
    """
    input:
        logbook=config["sample-logbook"],
    output:
        linker="results/export/linker.tsv",
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

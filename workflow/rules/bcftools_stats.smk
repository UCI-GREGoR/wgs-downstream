rule bcftools_stats:
    """
    Use bcftools stats to generate metrics
    for a vcf.
    """
    input:
        vcf="{prefix}.vcf.gz",
        tbi="{prefix}.vcf.gz.tbi",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
    output:
        stats="{prefix}.vcf.stats.txt",
    conda:
        "../envs/bcftools.yaml"
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools stats --threads {threads} --af-bins 0.01,0.05,0.1,1 -F {input.fasta} -s- {input.vcf} > {output.stats}"


rule plot_vcfstats:
    """
    Use plot-vcfstats to emit a wide array of results plots and files
    from bcftools stats output.
    """
    input:
        "{prefix}.vcf.stats.txt",
    output:
        "{prefix}/substitutions.0.png",
    params:
        outdir="{prefix}",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "plot-vcfstats -P -p {params.outdir} {input}"

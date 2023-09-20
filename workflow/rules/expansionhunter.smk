localrules:
    expansionhunter_rename_bai,


rule expansionhunter_rename_bai:
    """
    The path to the bam index is not exposed as a controllable
    parameter by expansionhunter. As such, copy the standard
    bai index to the modified extension .bam.bai.

    Due to ambiguity with other bai generation rules,
    this specifies both the bam and the bai with the
    original naming, though technically the bam is not required
    for this operation.
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bai",
    output:
        bai=temp("results/bams/{projectid}/{sampleid}.bam.bai"),
    shell:
        "cp {input.bai} {output.bai}"


rule expansionhunter_run:
    """
    Run ExpansionHunter using user-configured run settings.
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bam.bai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        linker="results/linker.tsv",
    output:
        "results/expansionhunter/{projectid}/{sampleid}.output.vcf",
        "results/expansionhunter/{projectid}/{sampleid}.output.json",
    params:
        variant_catalog="$CONDA_PREFIX/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json",
        output_prefix="results/expansionhunter/{projectid}/{sampleid}.output",
        region_extension_length=config["expansionhunter"]["region-extension-length"],
        sex=lambda wildcards: tc.get_sample_sex(wildcards, checkpoints),
        aligner=config["expansionhunter"]["aligner"],
        analysis_mode="seeking",
    conda:
        "../envs/expansionhunter.yaml"
    threads: 2
    resources:
        mem=2000,
        qname="small",
    shell:
        "ExpansionHunter --reads {input.bam} --reference {input.fasta} "
        "--variant-catalog {params.variant_catalog} "
        "--output-prefix {params.output_prefix} "
        "--region-extension-length {params.region_extension_length} "
        "--sex {params.sex} "
        "--aligner {params.aligner} "
        "--analysis-mode {params.analysis_mode} "
        "-n {threads}"


rule expansionhunter_filter_vcfs:
    """
    Filter expansionhunter raw output on PASS.
    """
    input:
        vcf="results/expansionhunter/{projectid}/{sampleid}.output.vcf",
    output:
        output="results/expansionhunter/{projectid}/{sampleid}.filtered.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "bcftools view --threads {threads} -f 'PASS' -Oz -o {output} {input.vcf}"


rule expansionhunter_combine_vcfs:
    """
    Merge all samples' filtered expansionhunter results together into one vcf.
    """
    input:
        vcf=lambda wildcards: tc.select_expansionhunter_subjects(
            wildcards,
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/expansionhunter",
            "filtered.vcf.gz",
            config["expansionhunter"]["excluded-samples"],
        ),
        tbi=lambda wildcards: tc.select_expansionhunter_subjects(
            wildcards,
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/expansionhunter",
            "filtered.vcf.gz.tbi",
            config["expansionhunter"]["excluded-samples"],
        ),
    output:
        "results/expansionhunter/expansionhunter_all_subjects.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=16000,
        qname="small",
    shell:
        "bcftools merge --threads {threads} -Oz -o {output} {input.vcf}"


rule expansionhunter_create_report:
    """
    Create report describing repeat length output found in expansionhunter runs.
    """
    input:
        results="results/expansionhunter/expansionhunter_all_subjects.vcf.gz",
        pubtable="resources/10.1038.nrg2828.tsv",
    output:
        "results/expansionhunter/expansionhunter_report.html",
    params:
        output_tsv="results/expansionhunter/expansionhunter_report.tsv",
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        qname="small",
    script:
        "../scripts/expansionhunter.Rmd"

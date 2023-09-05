rule telomerecat_bam_to_telbam:
    """
    Prune WGS bam to telomere read bams with telomerecat.
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bai",
    output:
        telbam="results/telomerecat/{projectid}/{sampleid}_telbam.bam",
    params:
        outdir="results/telomerecat/{projectid}",
    conda:
        "../envs/telomerecat.yaml"
    threads: 4
    resources:
        mem_mb=64000,
        qname="large",
    shell:
        "telomerecat bam2telbam -p {threads} --outbam_dir {params.outdir} {input.bam}"


rule telomerecat_telbam_to_csv:
    """
    Using the pruned telomere read bam, estimate telomere length with telomerecat.
    """
    input:
        bam="results/telomerecat/{projectid}/{sampleid}_telbam.bam",
    output:
        csv="results/telomerecat/{projectid}/{sampleid}_length.csv",
    params:
        max_handles=2,
    conda:
        "../envs/telomerecat.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        qname="large",
    shell:
        "telomerecat telbam2length -p {threads} --output {output.csv} {input.bam}"


localrules:
    telomerecat_run_all,


rule telomerecat_run_all:
    """
    Use dummy tracker file to trigger telomerecat runs for all valid bams.
    """
    input:
        linker="results/linker.tsv",
        csvs=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            checkpoints,
            bam_manifest["projectid"].to_list(),
            bam_manifest["sampleid"].to_list(),
            "results/telomerecat/",
            "_length.csv",
        ),
    output:
        "results/telomerecat/tracker.tsv",
    shell:
        "touch {output}"

rule run_multiqc_cross_flowcell:
    """
    Run multiqc across all input bams

    Currently limited to things that rely on pairwise data, so:
    - somalier
    """
    input:
        bam=manifest["internal-bam-paths"],
        multiqc_config=config["multiqc-config"],
    output:
        html="results/multiqc/multiqc.cross-flowcell.html",
        data_zip="results/multiqc/multiqc.cross-flowcell_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_cross_flowcell/out.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}/{projectid}",
                    toolname=["somalier"],
                    projectid=list(set(manifest["projectid"])),
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "-m somalier "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Cross-flowcell Data' "
        "-n {output.html}"

rule run_multiqc_cross_flowcell:
    """
    Run multiqc across all input bams

    Currently limited to things that rely on pairwise data, so:
    - somalier
    """
    input:
        somalier_relate="results/somalier/relate/somalier.html",
        somalier_ancestry="results/somalier/ancestry/results.somalier-ancestry.tsv",
        multiqc_config=tc.annotate_remote_file(config["multiqc-config"]),
    output:
        html="results/multiqc/multiqc.cross-flowcell.html",
        data_zip="results/multiqc/multiqc.cross-flowcell_data.zip",
    benchmark:
        "results/performance_benchmarks/run_multiqc_cross_flowcell/out.tsv"
    params:
        target_dirs=list(
            set(
                expand(
                    "results/{toolname}",
                    toolname=["somalier/ancestry", "somalier/relate"],
                )
            )
        ),
    conda:
        "../envs/multiqc.yaml"
    threads: config_resources["default"]["threads"]
    resources:
        mem_mb=config_resources["default"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["default"]["queue"], config_resources["queues"]
        ),
    shell:
        "multiqc {params.target_dirs} "
        "--config {input.multiqc_config} "
        "-m somalier "
        "--interactive "
        "--profile-runtime --zip-data-dir "
        "-f -i 'MultiQC for Cross-flowcell Data' "
        "-n {output.html}"

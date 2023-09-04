localrules:
    cyrius_create_manifest,
    cyrius_clone_repo,


rule cyrius_clone_repo:
    """
    Clone a local copy of Cyrius because
    there is no conda package for it, probably
    due to licensing.
    """
    output:
        directory("results/cyrius/repo"),
    params:
        github="git@github.com:Illumina/Cyrius.git",
        version=config["cyrius"]["version"],
    shell:
        "git clone {params.github} {output} && "
        "cd {output} && "
        "git checkout {params.version}"


rule cyrius_create_manifest:
    """
    Generate a simple plaintext file with *absolute*
    paths to input bams or crams.
    """
    input:
        bams=lambda wildcards: tc.select_cyrius_subjects(
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams",
            "bam",
        ),
    output:
        tsv="results/cyrius/input_manifest.tsv",
    run:
        with open(output.tsv, "w") as f:
            f.writelines(["{}\n".format(x) for x in input.bams])


rule cyrius_run:
    """
    Execute cyrius star caller for CYP2D6 alleles.
    """
    input:
        bams=lambda wildcards: tc.select_cyrius_subjects(
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams",
            "bam",
        ),
        bais=lambda wildcards: tc.select_cyrius_subjects(
            checkpoints,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams",
            "bai",
        ),
        tsv="results/cyrius/input_manifest.tsv",
        repo="results/cyrius/repo",
    output:
        tsv="results/cyrius/results.tsv",
        json="results/cyrius/results.json",
    params:
        outdir="results/cyrius",
        prefix="results",
    benchmark:
        "results/performance_benchmarks/cyrius_run/out.tsv"
    conda:
        "../envs/cyrius.yaml"
    threads: 4
    resources:
        mem_mb=16000,
        qname="large",
    shell:
        "python3 {input.repo}/star_caller.py --manifest {input.tsv} --genome 38 --prefix {params.prefix} --outDir {params.outdir} --threads {threads}"

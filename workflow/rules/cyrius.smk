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
        github="https://github.com/Illumina/Cyrius.git",
        version=config["cyrius"]["version"],
    shell:
        "git clone {params.github} {output} && "
        "cd {output} && "
        "git checkout {params.version}"


rule cyrius_create_manifest:
    """
    Generate a simple plaintext file with *absolute*
    paths to input bams or crams.

    Apparently, absolute in their docs actually means
    absolute or relative.
    """
    input:
        cram="results/crams/{sampleid}.cram",
    output:
        tsv=temp("results/cyrius/{sampleid}.input_manifest.tsv"),
    shell:
        "echo '{input}' > {output}"


rule cyrius_run:
    """
    Execute cyrius star caller for CYP2D6 alleles.
    """
    input:
        cram="results/crams/{sampleid}.cram",
        crai="results/crams/{sampleid}.crai",
        tsv="results/cyrius/{sampleid}.input_manifest.tsv",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        repo="results/cyrius/repo",
    output:
        tsv=temp("results/cyrius/{sampleid}.tsv"),
        json=temp("results/cyrius/{sampleid}.json"),
    params:
        outdir="results/cyrius",
        prefix="{sampleid}",
    benchmark:
        "results/performance_benchmarks/cyrius_run/{sampleid}.tsv"
    conda:
        "../envs/cyrius.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        qname="large",
    shell:
        "python3 {input.repo}/star_caller.py --manifest {input.tsv} --genome 38 --reference {input.fasta} "
        "--prefix {params.prefix} --outDir {params.outdir} --threads {threads}"


localrules:
    cyrius_combine_results,


rule cyrius_combine_results:
    """
    Concatenate cyrius outputs across all subjects.
    """
    input:
        tsvs=lambda wildcards: tc.select_cyrius_subjects(
            cram_manifest["sampleid"],
            config["cyrius"]["excluded-samples"],
            "results/cyrius",
            "tsv",
        ),
    output:
        "results/cyrius/cyrius_all_subjects.tsv",
    shell:
        "cat {input.tsvs} | awk 'NR == 1 || NR % 2 == 0' > {output}"


rule create_cyrius_report:
    """
    create report describing genotype distribution found in cyrius runs
    """
    input:
        results="results/cyrius/cyrius_all_subjects.tsv",
        pubtable="resources/10.1038.s41397-020-00205-5.tsv",
        somalier="results/somalier/ancestry/results.somalier-ancestry.tsv",
    output:
        "results/cyrius/cyrius_report.html",
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        time="1:00:00",
        mem=2000,
    script:
        "../scripts/cyrius.Rmd"

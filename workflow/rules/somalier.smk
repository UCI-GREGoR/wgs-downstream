rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        cram="results/crams/{sampleid}.cram",
        crai="results/crams/{sampleid}.crai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        sites_vcf="reference_data/somalier/{}/ref.sites.vcf.gz".format(reference_build),
    output:
        "results/somalier/extract/{sampleid}.somalier",
    benchmark:
        "results/performance_benchmarks/somalier_extract/{sampleid}.tsv"
    params:
        extract_dir="results/somalier/extract",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        qname="small",
    shell:
        "somalier extract -d {params.extract_dir} "
        "--sites {input.sites_vcf} "
        "-f {input.fasta} --sample-prefix {wildcards.sampleid}_ {input.cram} ; "
        'if [[ "$(samtools samples {input.cram} | cut -f 1)" != "{wildcards.sampleid}" ]] ; then '
        "mv {params.extract_dir}/$(samtools samples {input.cram} | cut -f 1).somalier {output} ; "
        "fi"


checkpoint somalier_relate:
    """
    Compute relatedness metrics on preprocessed alignment data with somalier.
    """
    input:
        somalier=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            cram_manifest["sampleid"].to_list(),
            "results/somalier/extract",
            ".somalier",
        ),
        ped="results/somalier/somalier.ped",
    output:
        html="results/somalier/relate/somalier.html",
        pairs="results/somalier/relate/somalier.pairs.tsv",
        samples="results/somalier/relate/somalier.samples.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_relate/out.tsv"
    params:
        outprefix="results/somalier/relate/somalier",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        "somalier relate --ped {input.ped} -o {params.outprefix} {input.somalier}"


checkpoint somalier_split_by_family:
    """
    Try to improve DAG solving performance by creating subsets of the pairs
    file that just contain results from a particular family
    """
    input:
        tsv="results/somalier/relate/somalier.pairs.tsv",
    output:
        pairs="results/somalier/relate/by_family/{family_id}.pairs.tsv",
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        'grep -E "#sample|-{wildcards.family_id}-" {input.tsv} > {output.pairs}'


rule somalier_build_pedfile:
    """
    Generate a pedfile for sex check for somalier.

    In the post 0.3.0 world, this is going to try to suck in self-reported sex
    information from the newly-generated sample ID linker. However, this information
    is only unreliably reported upstream, and as such this sexcheck will still be
    very incomplete. This is flagged to eventually be replaced with queries to retool,
    once that is actually implemented and operational.
    """
    input:
        sex_manifest=config["sample-sex"],
    output:
        ped="results/somalier/somalier.ped",
        problems="results/somalier/discordant_annotations.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_build_pedfile/somalier.tsv"
    params:
        sampleids=lambda wildcards: cram_manifest["sampleid"].to_list(),
        valid_subjectids=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            cram_manifest["sampleid"].to_list(),
            "",
            "",
        ),
        use_somalier_ids=True,
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    script:
        "../scripts/construct_somalier_pedfile.py"


rule somalier_get_reference_files:
    """
    Download somalier-provided extracted files
    """
    input:
        "reference_data/somalier/{}/kg.reference.data.tar.gz".format(reference_build),
    output:
        directory("results/somalier/references"),
    benchmark:
        "results/performance_benchmarks/somalier_get_reference_files/out.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "mkdir -p {output} && "
        "tar --directory={output} -z -x -v -f {input} && "
        "mv {output}/1kg-somalier/*somalier {output} && "
        "rmdir {output}/1kg-somalier"


rule somalier_ancestry:
    """
    Run ancestry inference with somalier
    """
    input:
        reference_labels="reference_data/somalier/{}/kg.labels.tsv".format(
            reference_build
        ),
        somalier_reference="results/somalier/references",
        somalier_experimental=lambda wildcards: tc.get_valid_subjectids(
            wildcards,
            cram_manifest["sampleid"].to_list(),
            "results/somalier/extract",
            ".somalier",
        ),
    output:
        "results/somalier/ancestry/results.somalier-ancestry.tsv",
        "results/somalier/ancestry/results.somalier-ancestry.html",
    params:
        outprefix="results/somalier/ancestry/results.",
    benchmark:
        "results/performance_benchmarks/somalier_ancestry/out.tsv"
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        "somalier ancestry --labels {input.reference_labels} -o {params.outprefix} "
        "{input.somalier_reference}/*somalier ++ {input.somalier_experimental}"


rule somalier_plot_pca:
    """
    The somalier default plot is nice for being interactive but doesn't
    make particularly attractive static plots.
    """
    input:
        somalier_ancestry="results/somalier/ancestry/results.somalier-ancestry.tsv",
    output:
        plotname="results/somalier/ancestry/results.somalier-ancestry.pcplot.png",
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        qname="small",
    script:
        "../scripts/plot_pca.R"

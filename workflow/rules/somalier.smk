rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        bam="results/bams/{projectid}/{sampleid}.bam",
        bai="results/bams/{projectid}/{sampleid}.bai",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        sites_vcf="reference_data/somalier/{}/ref.sites.vcf.gz".format(reference_build),
    output:
        "results/somalier/extract/{projectid}/{sampleid}.somalier",
    benchmark:
        "results/performance_benchmarks/somalier_extract/{projectid}/{sampleid}.tsv"
    params:
        extract_dir="results/somalier/extract/{projectid}",
    conda:
        "../envs/somalier.yaml"
    threads: 1
    resources:
        mem_mb="2000",
        qname="small",
    shell:
        "somalier extract -d {params.extract_dir} "
        "--sites {input.sites_vcf} "
        "-f {input.fasta} --sample-prefix {wildcards.sampleid}_ {input.bam} && "
        "mv {params.extract_dir}/$(samtools samples {input.bam} | cut -f 1).somalier {output}"


def get_valid_pmgrcs(wildcards, projectids, sampleids, prefix, suffix):
    valid_samples = []
    valid_projects = []
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    for projectid, sampleid in zip(projectids, sampleids):
        df_matches = df.loc[(df["ru"] == projectid) & (df["sq"] == sampleid), "pmgrc"]
        if len(df_matches) == 1:
            valid_samples.append(df_matches.to_list()[0])
            valid_projects.append(projectid)
    res = [
        "{}{}/{}{}".format(prefix, projectid, sampleid, suffix)
        for projectid, sampleid in zip(valid_projects, valid_samples)
    ]
    return res


rule somalier_relate:
    """
    Compute relatedness metrics on preprocessed alignment data with somalier.
    """
    input:
        somalier=lambda wildcards: get_valid_pmgrcs(
            wildcards,
            manifest["projectid"].to_list(),
            manifest["sampleid"].to_list(),
            "results/somalier/extract/",
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
        mem_mb="4000",
        qname="small",
    shell:
        "somalier relate --ped {input.ped} -o {params.outprefix} -i {input.somalier}"


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
        linker="results/linker.tsv",
    output:
        ped="results/somalier/somalier.ped",
        problems="results/somalier/discordant_annotations.tsv",
    benchmark:
        "results/performance_benchmarks/somalier_build_pedfile/somalier.tsv"
    params:
        projectids=lambda wildcards: manifest["projectid"].to_list(),
        subjectids=lambda wildcards: manifest["sampleid"].to_list(),
        valid_pmgrcids=lambda wildcards: get_valid_pmgrcs(
            wildcards,
            manifest["projectid"].to_list(),
            manifest["sampleid"].to_list(),
            "",
            "",
        ),
    threads: 1
    resources:
        mem_mb="1000",
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
        mem_mb="1000",
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
        somalier_experimental=lambda wildcards: get_valid_pmgrcs(
            wildcards,
            manifest["projectid"].to_list(),
            manifest["sampleid"].to_list(),
            "results/somalier/extract/",
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
        mem_mb="4000",
        qname="small",
    shell:
        "somalier ancestry --labels {input.reference_labels} -o {params.outprefix} "
        "{input.somalier_reference}/*somalier ++ {input.somalier_experimental}"

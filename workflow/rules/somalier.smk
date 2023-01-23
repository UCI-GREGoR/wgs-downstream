def link_bams_by_id(wildcards, suffix):
    outfn = str(checkpoints.generate_linker.get().output[0])
    linker = pd.read_csv(outfn, sep="\t")
    projectid = linker.loc[linker["pmgrc"] == wildcards.sampleid, "ru"]
    sampleid = linker.loc[linker["pmgrc"] == wildcards.sampleid, "sq"]
    return "results/bams/{}/{}.{}".format(projectid, sampleid, suffix)


rule somalier_extract:
    """
    Run somalier extract on a single bam.
    """
    input:
        bam=lambda wildcards: link_bams_by_id(wildcards, "bam"),
        bai=lambda wildcards: link_bams_by_id(wildcards, "bai"),
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
        mem_mb="2000",
        qname="small",
    shell:
        "somalier extract -d {params.extract_dir} "
        "--sites {input.sites_vcf} "
        "-f {input.fasta} "
        "{input.bam}"


def get_valid_pmgrcs(wildcards, projectids, sampleids):
    outfn = str(checkpoints.generate_linker.get().output[0])
    linker = pd.read_csv(outfn, sep="\t")
    valid_ids = []
    for ru, sq, pmgrc in zip(linker["ru"], linker["sq"], linker["pmgrc"]):
        if (ru in projectids) and (sq in sampleids):
            valid_ids.append(pmgrc)
    return valid_ids


rule somalier_relate:
    """
    Compute relatedness metrics on preprocessed alignment data with somalier.
    """
    input:
        somalier=lambda wildcards: expand(
            "results/somalier/extract/{sampleid}.somalier",
            sampleid=get_valid_pmgrcs(
                wildcards, manifest["projectid"], manifest["sampleid"]
            ),
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
        "somalier relate --ped {input.ped} -o {params.outprefix} {input.somalier}"


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
        pmgrcids=lambda wildcards: get_valid_pmgrcs(
            wildcards, manifest["projectid"], manifest["sampleid"]
        ),
    threads: 1
    resources:
        mem_mb="1000",
        qname="small",
    script:
        "../scripts/construct_somalier_pedfile.py"

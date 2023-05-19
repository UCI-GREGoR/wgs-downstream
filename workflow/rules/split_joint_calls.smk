checkpoint get_sample_list_from_vcf:
    """
    Extract actual sample list from joint called vcf
    """
    input:
        "results/glnexus/{subset}/merged_callset.filtered.regions.vcf.gz",
    output:
        temp("results/glnexus/{subset}/joint_called_samples.tsv"),
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        qname="small",
    shell:
        "bcftools query -l {input} > {output}"


def get_family_clusters(checkpoints):
    """
    From the joint call set, compute the
    set of family clusters present in the
    filtered vcf
    """
    subject_ids = ""
    with checkpoints.get_sample_list_from_vcf.get().output[0].open() as f:
        subject_ids = f.readlines()
    family_ids = []
    for subject_id in subject_ids:
        split_id = subject_id.split("-")
        if len(split_id) == 4:
            family_ids.append(split_id[2])
    return family_ids


def get_subjects_by_family(checkpoints, family_id):
    """
    From the joint call set, find the set of
    subjects that belong to a specified family
    cluster
    """
    subject_ids = ""
    with checkpoints.get_sample_list_from_vcf.get().output[0].open() as f:
        subject_ids = f.readlines()
    family_subjects = []
    for subject_id in subject_ids:
        split_id = subject_id.split("-")
        if len(split_id) == 4:
            if split_id[2] == family_id:
                family_subjects.append(subject_id)
    return family_subjects


localrules:
    create_family_list,
    aggregate_split_joint_calls,


rule create_family_list:
    """
    For a family cluster, list the set
    of subjects in the joint call that
    belong to that cluster and write to
    file for use with bcftools
    """
    input:
        "results/glnexus/{family_cluster}/joint_called_samples.tsv",
    output:
        temp("results/split_joint_calls/{family_cluster}.tsv"),
    shell:
        'awk -v var={wildcards.family_cluster} \'{{pat = "-"var"-[0-9]$" ; if ($1 ~ pat) {{print $1}}}}\' {input} > {output}'


rule split_joint_calls_by_family:
    """
    Create per-sample vcfs with filtering logic
    applied to family clusters to try to balance
    information, redundancy, and quality
    """
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.vcf.gz",
        tsv="results/split_joint_calls/{family_cluster}.tsv",
    output:
        vcf=temp(
            "results/split_joint_calls/PMGRC-{subject_id}-{family_cluster}-{relationship}_family-variants.vcf.gz"
        ),
    params:
        subject_id="PMGRC-{subject_id}-{family_cluster}-{relationship}",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        "bcftools view --threads {threads} -S {input.tsv} -O u {input.vcf} | "
        "bcftools view --threads {threads} -c 1 -O u | "
        "bcftools view --threads {threads} -s {params.subject_id} -O z -o {output}"


rule split_joint_calls_single_sample:
    """
    Provide an alternative filtering path for samples with no family identifier
    """
    input:
        vcf="results/glnexus/NA{giab_code}/merged_callset.filtered.regions.vcf.gz",
    output:
        vcf=temp("results/split_joint_calls/NA{giab_code}_family-variants.vcf.gz"),
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        qname="small",
    shell:
        "bcftools view --threads {threads} -s NA{wildcards.giab_code} -O u {input.vcf} | "
        "bcftools view --threads {threads} -c 1 -O u -o {output.vcf}"


rule remove_all_problematic_regions:
    """
    At downstream user request, be even more
    stringent than usual and remove anything
    that NIST annotates as a challenging region
    for sequencing. The primary target of this
    filter is highly repetitive content.
    """
    input:
        vcf="results/split_joint_calls/{subject_id}_family-variants.vcf.gz",
        regions="resources/GRCh38_alllowmapandsegdupregions.bed",
    output:
        temp("results/split_joint_calls/{subject_id}_nist-filters.vcf.gz"),
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.vcf} -b {input.regions} -wa -v -header | bgzip -c > {output}"


rule apply_single_sample_qc:
    """
    Mirroring behavior from the upstream
    pipeline's export preparation,
    apply sample-specific allelic balance
    filters, but in this case set the genotypes
    to missing when the filter fails instead
    of excluding the genotype entirely.
    """
    input:
        "results/split_joint_calls/{subject_id}_nist-filters.vcf.gz",
    output:
        "results/split_joint_calls/{subject_id}_{lsid}_{sqid}.snv.vcf.gz",
    params:
        exportid="{subject_id}_{lsid}_{sqid}",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools +setGT -O v {input} -- -t q -n . "
        '-e \'((FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.2 & FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.8 & GT != "1/1") | '
        ' (GT == "0/0") | '
        ' (FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.05 & GT = "1/1"))\' | '
        'bcftools reheader -s <(echo -e "{wildcards.subject_id}\\t{params.exportid}") | '
        "bcftools sort -O v | "
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


def compute_expected_single_samples(manifest, checkpoints):
    valid_samples = {}
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    for projectid, sampleid in zip(manifest["projectid"], manifest["sampleid"]):
        df_matches = df.loc[(df["ru"] == projectid) & (df["sq"] == sampleid), :]
        if len(df_matches["pmgrc"]) == 1:
            ## ad hoc handler for rerun samples: pick the later project id
            if df_matches["pmgrc"].to_list()[0] in valid_samples.keys():
                if projectid > valid_samples[df_matches["pmgrc"].to_list()[0]]:
                    valid_samples[df_matches["pmgrc"].to_list()[0]] = (
                        projectid,
                        "{}_{}_{}".format(
                            df_matches["pmgrc"].to_list()[0],
                            df_matches["ls"].to_list()[0],
                            df_matches["sq"].to_list()[0],
                        ),
                    )
            else:
                valid_samples[df_matches["pmgrc"].to_list()[0]] = (
                    projectid,
                    "{}_{}_{}".format(
                        df_matches["pmgrc"].to_list()[0],
                        df_matches["ls"].to_list()[0],
                        df_matches["sq"].to_list()[0],
                    ),
                )
    res = [
        "results/split_joint_calls/{}.snv.vcf.gz".format(x[1][1])
        for x in valid_samples.items()
    ]
    return res


rule aggregate_split_joint_calls:
    """
    Determine the full set of expected output vcfs
    from a single sample split of a joint call run
    """
    input:
        lambda wildcards: compute_expected_single_samples(gvcf_manifest, checkpoints),
    output:
        temp("results/split_joint_calls/.joint_calls_split"),
    shell:
        "touch {output}"

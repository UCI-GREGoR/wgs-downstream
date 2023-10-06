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
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools query -l {input} > {output}"


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
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
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
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
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
    threads: config_resources["bedtools"]["threads"]
    resources:
        mem_mb=config_resources["bedtools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bedtools"]["queue"], config_resources["queues"]
        ),
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
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=lambda wildcards: rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools +setGT -O v {input} -- -t q -n . "
        '-e \'((FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.2 & FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.8 & GT != "1/1") | '
        ' (GT == "0/0") | '
        ' (FORMAT/AD[0:0] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.05 & GT = "1/1"))\' | '
        'bcftools reheader -s <(echo -e "{wildcards.subject_id}\\t{params.exportid}") | '
        "bcftools sort -O v | "
        "sed 's|\\t1/0:|\\t0/1:|' | bgzip -c > {output}"


rule aggregate_split_joint_calls:
    """
    Determine the full set of expected output vcfs
    from a single sample split of a joint call run
    """
    input:
        lambda wildcards: tc.compute_expected_single_samples(gvcf_manifest),
    output:
        temp("results/split_joint_calls/.joint_calls_split"),
    shell:
        "touch {output}"

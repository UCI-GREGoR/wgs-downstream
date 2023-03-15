localrules:
    glnexus_create_gvcf_list,


rule glnexus_create_gvcf_list:
    """
    Emit target set of gvcf filenames to file in preparation
    for glnexus
    """
    input:
        gvcfs=gvcf_manifest["gvcf"],
    output:
        tsv="results/glnexus/gvcf_list.tsv",
    run:
        with open(output.tsv, "w") as f:
            f.writelines([x + "\n" for x in input.gvcfs])


def get_calling_range_by_chrom(wildcards: Wildcards, ranges: str):
    """
    Report a particular calling range as specified by chromosome
    """
    all_ranges = []
    with open(ranges, "r") as f:
        all_ranges = [x.rstrip() for x in f.readlines()]
    for fn in all_ranges:
        if re.search("{}.bed".format(wildcards.chrom), fn) is not None:
            return fn


def get_all_calling_ranges(ranges: str):
    """
    Report all calling ranges from user-specified range list file
    """
    all_ranges = []
    with open(ranges, "r") as f:
        all_ranges = [x.rstrip() for x in f.readlines()]
    res = []
    for fn in all_ranges:
        current_match = re.match("^.*/(.*)\\.bed$", fn)
        if current_match is None:
            raise ValueError("unable to find match in {}".format(fn))
        else:
            res.append(current_match[1])
    return res


rule glnexus_joint_calling:
    """
    Given gvcfs, create a joint called dataset.
    """
    input:
        tsv="results/glnexus/gvcf_list.tsv",
        calling_ranges=lambda wildcards: get_calling_range_by_chrom(
            wildcards, config["glnexus"]["calling-ranges"]
        ),
    output:
        bcf=temp("results/glnexus/merged_callset_{chrom}.bcf"),
        tmp=temp(directory("results/glnexus/tmp_{chrom}")),
    params:
        gvcf_manifest=gvcf_manifest,
        memlimit="16",
    conda:
        "../envs/glnexus.yaml"
    threads: 4
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "cat {input.tsv} | "
        "LD_PRELOAD=$(jemalloc-config --libdir)/libjemalloc.so.$(jemalloc-config --revision) "
        "xargs glnexus_cli --config DeepVariant "
        "--bed {input.calling_ranges} -a -m {params.memlimit} "
        "-t {threads} --dir {output.tmp} && touch {output}"


rule prepare_joint_calling_output:
    """
    glnexus emits unconditional bcfs; use bcftools to create gzipped vcf
    """
    input:
        bcf=lambda wildcards: expand(
            "results/glnexus/merged_callset_{chrom}.bcf",
            chrom=get_all_calling_ranges(config["glnexus"]["calling-ranges"]),
        ),
    output:
        vcf="results/glnexus/merged_callset.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    resources:
        mem_mb="4000",
        qname="small",
    shell:
        "bcftools concat {input.bcf} --threads {threads} -O z -o {output.vcf}"

import os
import re
from warnings import warn

import pandas as pd
from snakemake.io import AnnotatedString, Namedlist, Wildcards
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider()
HTTP = HTTPRemoteProvider()


def annotate_remote_file(fn: str):
    """
    try our best to automagically wrap files with remote handlers.
    the FTP remote handler, as of early 2023, is a bit of a mess
    that breaks everything, so it's not included here. note that
    rules using the S3 remote with environment SSO caching should
    be set to high priority, lest the SSO session time out while
    the rule is queued.
    """
    if fn.startswith("https://") or fn.startswith("http://"):
        return HTTP.remote(fn)
    if fn.startswith("s3://"):
        return S3.remote(fn)
    return fn


def format_reference_build(ref: str) -> str:
    """
    Given a config key for a reference genome, try to turn it
    into a string that Moon will acknowledge

    So here's the thing: Moon is a very strange tool that implicitly
    enforces a bunch of restrictions on otherwise valid vcf and then
    doesn't bother to tell you about it. One of those restrictions is
    for its support of GRCh38: it tries to determine reference genome
    with the following logic:

    - if it finds a vcf header entry that contains anywhere in the line
      '##reference=', it takes the rest of the line after those characters
      and searches exactly for 'hg19' or 'GRCh37'; and then if that fails
      it searches for 'hg38' or 'GRCh38'
    - if the above fails for any reason, it then sniffs the entire header
      first for 'GRCh38' or 'hg38'
    - if the above fails, it returns 37

    Note that the '##reference' search is not anchored, and the 'GRCh' searches
    are case-sensitive
    """
    if ref.lower() == "grch38" or ref.lower() == "hg38":
        return "GRCh38"
    if ref.lower() == "grch37" or ref.lower() == "hg19":
        return "GRCh37"
    warn(
        "The user-provided genome build configuration key '"
        + ref
        + "' does not match any of the simple patterns expected by Moon; "
        "as such, this will cause mismatch problems if you try to load "
        "a variant VCF into Moon without first modifying the header to "
        "contain 'GRCh38' or 'hg38' somewhere in the header. If your "
        "intended build is GRCh37/hg19, Moon should default to this anyway."
    )
    return ref


def map_reference_file(wildcards: Namedlist, config: dict) -> str | AnnotatedString:
    """
    Use wildcard information to figure out what configured
    reference file is needed, and then wrap that file in a
    remote provider structure as required.
    """
    queries = wildcards.reference_file.split("/")
    queries[len(queries) - 1] = (
        queries[len(queries) - 1].removeprefix("ref.").replace(".", "-")
    )
    current_lvl = config
    for query in queries:
        current_lvl = current_lvl[query]
    ## The intention for this function was to distinguish between S3 file paths and others,
    ## and return wrapped objects related to the remote provider service when appropriate.
    ## There have been periodic issues with the remote provider interface, and the FTP one
    ## is completely unusable with conda env creation due to timeouts, so I'm giving up
    return current_lvl


def get_fid(sampleid):
    """
    From a subject ID, extract the effective family cluster identifier.
    For PGMRC IDs, this will be PMGRC-{numeric sample id}-{numeric family id}-{relationship code}.
    For other IDs, at this point, just treat the ID itself as its own cluster.
    """
    split_id = sampleid.split("-")
    if len(split_id) == 1:
        return split_id[0]
    else:
        return split_id[2]


def get_valid_pmgrcs(wildcards, checkpoints, projectids, sampleids, prefix, suffix):
    valid_samples = {}
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    for projectid, sampleid in zip(projectids, sampleids):
        df_matches = df.loc[(df["ru"] == projectid) & (df["sq"] == sampleid), "pmgrc"]
        if len(df_matches) == 1:
            ## ad hoc handler for rerun samples: pick the later project id
            if df_matches.to_list()[0] in valid_samples.keys():
                if projectid > valid_samples[df_matches.to_list()[0]]:
                    valid_samples[df_matches.to_list()[0]] = projectid
            else:
                valid_samples[df_matches.to_list()[0]] = projectid
    res = []
    for sampleid, projectid in valid_samples.items():
        if "subset" in wildcards:
            fid = get_fid(sampleid)
            if wildcards.subset == "all" or wildcards.subset == fid:
                res.append("{}{}/{}{}".format(prefix, projectid, sampleid, suffix))
        else:
            res.append("{}{}/{}{}".format(prefix, projectid, sampleid, suffix))
    return res


def link_bams_by_id(wildcards, checkpoints, bam_manifest):
    sampleid = ""
    projectid = wildcards.projectid
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    df_sampleid = df.loc[
        (df["pmgrc"] == wildcards.sampleid) & (df["ru"] == projectid), "sq"
    ]
    if len(df_sampleid) == 1:
        sampleid = df_sampleid.to_list()[0]
    else:
        raise ValueError(
            "cannot find pmgrc id in manifest: {}".format(wildcards.sampleid)
        )
    res = bam_manifest.loc[
        (bam_manifest["sampleid"] == sampleid)
        & (bam_manifest["projectid"] == projectid),
        "bam",
    ]
    return [annotate_remote_file(x) for x in res]


def link_gvcfs_by_id(wildcards, checkpoints, gvcf_manifest):
    sampleid = ""
    projectid = wildcards.projectid
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    df_sampleid = df.loc[
        (df["pmgrc"] == wildcards.sampleid) & (df["ru"] == projectid), "sq"
    ]
    if len(df_sampleid) == 1:
        sampleid = df_sampleid.to_list()[0]
    else:
        raise ValueError(
            "cannot find pmgrc id in manifest: {}".format(wildcards.sampleid)
        )
    res = gvcf_manifest.loc[
        (gvcf_manifest["sampleid"] == sampleid)
        & (gvcf_manifest["projectid"] == projectid),
        "gvcf",
    ]
    return [annotate_remote_file(x) for x in res]


def select_expansionhunter_denovo_subjects(
    wildcards, checkpoints, projectids, sampleids, prefix, suffix
) -> list:
    """
    Access checkpoint output to determine the subset of input manifest
    subjects that also have available phenotype information in the input
    data model
    """
    fn = str(checkpoints.expansionhunter_denovo_create_manifest.get().output[0])
    df = pd.read_table(fn, sep="\t", header=None, names=["sample", "status", "json"])
    outfn = str(checkpoints.generate_linker.get().output[0])
    linker = pd.read_table(outfn, sep="\t")
    res = []
    for projectid, sampleid in zip(projectids, sampleids):
        linker_sampleid = linker.loc[
            (linker["ru"] == projectid) & (linker["sq"] == sampleid), "pmgrc"
        ]
        if len(linker_sampleid) == 1:
            pmgrcid = linker_sampleid.to_list()[0]
            if pmgrcid in df["sample"].to_list():
                res.append("{}/{}/{}.{}".format(prefix, projectid, pmgrcid, suffix))
    return res


def get_calling_range_by_chrom(wildcards: Wildcards, ranges: str):
    """
    Report a particular calling range as specified by chromosome
    """
    all_ranges = []
    with open(ranges, "r") as f:
        all_ranges = [x.rstrip() for x in f.readlines()]
    for fn in all_ranges:
        if re.search("{}.bed".format(wildcards.chrom), fn) is not None:
            return annotate_remote_file(fn)


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


def get_probands_with_structure(checkpoints):
    """
    From the joint call set, get proband-flagged subjects
    with at least one parent also present in the dataset
    """
    subject_ids = ""
    with checkpoints.get_sample_list_from_vcf.get().output[0].open() as f:
        subject_ids = f.readlines()
    parents = {}
    children = {}
    results = []
    for subject_id in subject_ids:
        split_id = subject_id.split("-")
        if len(split_id) == 4:
            if split_id[3] == "0":
                children[split_id[2]] = subject_id
            elif split_id[3] == "1" or split_id[3] == "2":
                parents[split_id[2] + "-" + split_id[3]] = subject_id
    for cluster, child in children:
        if (
            "{}-1".format(cluster) in parents.keys()
            or "{}-2".format(cluster) in parents.keys()
        ):
            results.append(child)
    return results


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

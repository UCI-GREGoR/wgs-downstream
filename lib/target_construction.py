import os
import pathlib
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


def get_valid_subjectids(wildcards, checkpoints, projectids, sampleids, prefix, suffix):
    valid_samples = {}
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    for projectid, sampleid in zip(projectids, sampleids):
        df_matches = df.loc[(df["ru"] == projectid) & (df["sq"] == sampleid), "subject"]
        if len(df_matches) == 0:
            ## the "project ID" is not straightforwardly conveyed when specifying new form IDs;
            ## as a temporary workaround, allow unique matches when linker project and sq IDs are NA
            df_matches = df.loc[
                (df["sq"].isna()) & (df["ru"].isna()) & (df["subject"] == sampleid),
                "subject",
            ]
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
        (df["subject"] == wildcards.sampleid) & (df["ru"] == projectid), "sq"
    ]
    if len(df_sampleid) == 0:
        ## the "project ID" is not straightforwardly conveyed when specifying new form IDs;
        ## as a temporary workaround, allow unique matches when linker project and sq IDs are NA
        df_sampleid = df.loc[
            (df["sq"].isna())
            & (df["ru"].isna())
            & (df["subject"] == wildcards.sampleid),
            "subject",
        ]
    if len(df_sampleid) == 1:
        sampleid = df_sampleid.to_list()[0]
    else:
        raise ValueError(
            "cannot find subject id in manifest: {}, {}".format(
                wildcards.sampleid, wildcards.projectid
            )
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
        (df["subject"] == wildcards.sampleid) & (df["ru"] == projectid), "sq"
    ]
    if len(df_sampleid) == 0:
        ## the "project ID" is not straightforwardly conveyed when specifying new form IDs;
        ## as a temporary workaround, allow unique matches when linker project and sq IDs are NA
        df_sampleid = df.loc[
            (df["sq"].isna())
            & (df["ru"].isna())
            & (df["subject"] == wildcards.sampleid),
            "subject",
        ]
    if len(df_sampleid) == 1:
        sampleid = df_sampleid.to_list()[0]
    else:
        raise ValueError(
            "cannot find subject id in manifest: {}".format(wildcards.sampleid)
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
            (linker["ru"] == projectid) & (linker["sq"] == sampleid), "subject"
        ]
        if len(linker_sampleid) == 1:
            subjectid = linker_sampleid.to_list()[0]
            if subjectid in df["sample"].to_list():
                res.append("{}/{}/{}.{}".format(prefix, projectid, subjectid, suffix))
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
    with checkpoints.get_sample_list_from_vcf.get(subset="all").output[0].open() as f:
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


def get_subjects_by_family(
    wildcards,
    checkpoints,
    family_id,
    relationship_code,
    projectids,
    sampleids,
    prefix,
    suffix,
):
    """
    From the joint call set, find the set of
    subjects that belong to a specified family
    cluster
    """
    subject_ids = get_valid_subjectids(
        wildcards, checkpoints, projectids, sampleids, prefix, suffix
    )
    family_subjects = []
    for subject_id in subject_ids:
        split_id = pathlib.PurePosixPath(subject_id).name.split("-")
        if len(split_id) == 4:
            if split_id[2] == str(family_id) and re.search(
                r"^{}[\._]".format(relationship_code), split_id[3]
            ):
                family_subjects.append(subject_id)
    return family_subjects


def compute_expected_single_samples(manifest, checkpoints):
    valid_samples = {}
    outfn = str(checkpoints.generate_linker.get().output[0])
    df = pd.read_table(outfn, sep="\t")
    for projectid, sampleid in zip(manifest["projectid"], manifest["sampleid"]):
        df_matches = df.loc[(df["ru"] == projectid) & (df["sq"] == sampleid), :]
        if len(df_matches["subject"]) == 1:
            ## ad hoc handler for rerun samples: pick the later project id
            if df_matches["subject"].to_list()[0] in valid_samples.keys():
                if projectid > valid_samples[df_matches["subject"].to_list()[0]]:
                    valid_samples[df_matches["subject"].to_list()[0]] = (
                        projectid,
                        "{}_{}_{}".format(
                            df_matches["subject"].to_list()[0],
                            df_matches["ls"].to_list()[0],
                            df_matches["sq"].to_list()[0],
                        ),
                    )
            else:
                valid_samples[df_matches["subject"].to_list()[0]] = (
                    projectid,
                    "{}_{}_{}".format(
                        df_matches["subject"].to_list()[0],
                        df_matches["ls"].to_list()[0],
                        df_matches["sq"].to_list()[0],
                    ),
                )
    res = [
        "results/split_joint_calls/{}.snv.vcf.gz".format(x[1][1])
        for x in valid_samples.items()
    ]
    return res


def caller_interval_file_count(config: dict) -> int:
    """
    Get simple line count of a file; this is intended to
    count a tiny text file containing <100 interval filenames.
    """
    fn = config["deeptrio"][config["genome-build"]]["calling-ranges"]
    x = 0
    with open(fn, "r") as f:
        x = len(f.readlines())
    return x


def determine_trio_structure(
    wildcards, checkpoints, config, manifest, sampleid, chrcode
) -> str:
    """
    Use somalier relatedness/sexcheck output and expected pedigree to determine
    which version of deeptrio should be run for a particular chromosome.

    somalier IDs have format "PMGRC..._SQ" for whatever reason. this creates theoretical
    ambiguities

    somalier.samples.tsv column 25 "Y_n" should be 0 for female samples, >0 for male

    somalier.pairs.tsv column 3 "relatedness" should be >= 0.3535534

    chrcode is the index of the line in the tsv file defining the calling regions
    """
    somalier_samples = str(checkpoints.somalier_relate.get().output["samples"])
    somalier_pairs = str(checkpoints.somalier_relate.get().output["pairs"])

    subject_ids = get_valid_subjectids(
        wildcards, checkpoints, manifest["projectid"], manifest["sampleid"], "", ""
    )
    family_id = sampleid.split("-")[3]
    mother_id = ""
    father_id = ""
    for subject_id in subject_ids:
        split_id = pathlib.PurePosixPath(subject_id).name.split("-")
        if len(split_id) == 4:
            if split_id[2] == str(family_id) and re.search(r"^1[\._]", split_id[3]):
                father_id = subject_id
            if split_id[2] == str(family_id) and re.search(r"^2[\._]", split_id[3]):
                mother_id = subject_id

    # check for parent relatedness in somalier output
    somalier_pairs_df = pd.read_table(somalier_pairs, sep="\t")
    somalier_pairs_df[["sample_a", "suffix_a"]] = somalier_pairs_df[
        "#sample_a"
    ].str.split("_", n=1, expand=True)
    somalier_pairs_df[["sample_b", "suffix_b"]] = somalier_pairs_df[
        "sample_b"
    ].str.split("_", n=1, expand=True)
    somalier_pairs_df = somalier_pairs_df[
        (somalier_pairs_df["sample_a"] == sampleid)
        | (somalier_pairs_df["sample_b"] == sampleid)
    ]
    mother_relatedness = somalier_pairs_df[
        (somalier_pairs_df["sample_a"] == mother_id)
        | (somalier_pairs_df["sample_b"] == mother_id)
    ]
    if mother_relatedness["relatedness"][0] < 0.3535534:
        mother_id = ""
    father_relatedness = somalier_pairs_df[
        (somalier_pairs_df["sample_a"] == father_id)
        | (somalier_pairs_df["sample_b"] == father_id)
    ]
    if father_relatedness["relatedness"][0] < 0.3535534:
        father_id = ""
    if mother_id == "" and father_id == "":
        raise ValueError("cannot solve family structure for {}".format(sampleid))

    # check for proband sex in somalier output
    somalier_samples_df = pd.read_table(somalier_samples, sep="\t")
    somalier_samples_df[["sample", "suffix"]] = somalier_samples_df[
        "sample_id"
    ].str.split("_", n=1, expand=True)
    somalier_samples_df = somalier_samples_df.set_index("sample")
    sample_is_male = somalier_samples_df.loc[sampleid, "Y_n"][0] < 1

    if mother_id == "":
        return "father_only"
    if father_id == "":
        return "mother_only"

    # use chrcode as line index to probe region file and determine whether it's an X/Y/par nonsense
    fn = config["deeptrio"][config["genome-build"]]["calling-ranges"]
    lines = []
    with open(fn, "r") as f:
        lines = f.readlines()
    interval = lines[chrcode].rstrip()
    if "chrX" in interval:
        if "non-PAR" in interval and sample_is_male:
            return "mother_only"
    if "chrY" in interval:
        if "non-PAR" in interval and sample_is_male:
            return "father_only"
    return "full_trio"

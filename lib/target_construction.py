import os
from warnings import warn

import pandas as pd
from snakemake.io import AnnotatedString, Namedlist
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

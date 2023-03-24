import os

import pandas as pd
from snakemake.io import AnnotatedString, Namedlist
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider()
GS = GSRemoteProvider()
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
    if fn.startswith("gs://"):
        return GS.remote(fn)
    if fn.startswith("s3://"):
        return S3.remote(fn)
    return fn


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

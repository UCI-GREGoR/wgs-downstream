rule create_sequence_dictionary:
    """
    For a reference fasta, create a sequence dictionary (.dict extension)

    Because GATK is ridiculous, make two copies of the dict, so it's always available
    no matter how the downstream tool thinks it should be named.
    """
    input:
        "{prefix}fasta",
    output:
        standard="{prefix}fasta.dict",
        modified="{prefix}dict",
    benchmark:
        "results/performance_benchmarks/create_sequence_dictionary/{prefix}fasta.tsv"
    params:
        tmpdir=tempDir,
        java_args="-Djava.io.tmpdir={} -XX:CompressedClassSpaceSize=200m -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx2000m".format(
            tempDir
        ),
    conda:
        "../envs/gatk4.yaml"
    threads: 1
    resources:
        mem_mb="10000",
        qname="small",
        tmpdir=tempDir,
    shell:
        "mkdir -p {params.tmpdir} && "
        'gatk --java-options "{params.java_args}" CreateSequenceDictionary '
        "-REFERENCE {input} "
        "-OUTPUT {output.standard} "
        "--TMP_DIR {params.tmpdir} && "
        "cp {output.standard} {output.modified}"

localrules:
    caller_scatter_tasks,


rule caller_scatter_tasks:
    """
    For scatter-gather: take input set of range annotations
    and prepare for individual caller runs
    """
    input:
        lambda wildcards: config["deeptrio"][reference_build]["calling-ranges"],
    output:
        "results/{caller}/split_ranges/{splitnum}.ssv",
    benchmark:
        "results/performance_benchmarks/{caller}/split_ranges/{splitnum}.tsv"
    shell:
        "cat $(awk 'NR == {wildcards.splitnum}' {input}) | tr '\\n' ' ' > {output}"


rule deeptrio_make_examples_full_trio:
    """
    Run deeptrio make_examples in a hybrid
    embarrassingly parallel fashion. This rule
    will apply for: autosomes, X-PAR, Y-PAR,
    X-nonPAR for female proband
    """
    input:
        child_bam="results/bams/{projectid0}/PMGRC-{childid}-{childid}-0.bam",
        parent1_bam=lambda wildcards: tc.get_subjects_by_family(
            wildcards,
            checkpoints,
            wildcards.childid,
            1,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams/",
            ".bam",
        ),
        parent2_bam=lambda wildcards: tc.get_subjects_by_family(
            wildcards,
            checkpoints,
            wildcards.childid,
            2,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams/",
            ".bam",
        ),
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deeptrio/split_ranges/{splitnum}.ssv",
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    output:
        examples=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/full_trio/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                relation=["child", "parent1", "parent2"],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        gvcfs=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/full_trio/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
                relation=["child", "parent1", "parent2"],
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        jsons=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/full_trio/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz.example_info.json",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_make_examples/full_trio/{projectid0}/PMGRC-{childid}-{childid}-0.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/full_trio/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        gvcf_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/full_trio/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        tmpdir="/tmp",
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["make_examples_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
        tmpdir="/tmp",
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "mkdir -p {params.tmpdir} && '
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir {params.tmpdir} "
        "make_examples --mode calling "
        "--ref {input.fasta} "
        "--reads {input.child_bam} --reads_parent1 {input.parent1_bam} --reads_parent2 {input.parent2_bam} "
        '--regions \\"$(cat {input.intervals})\\" '
        "--examples {params.shard_string} --channels insert_size "
        "--gvcf {params.gvcf_string} "
        '--task {{}}"'


rule deeptrio_make_examples_mother_only:
    """
    Run deeptrio make_examples in a hybrid
    embarrassingly parallel fashion. This rule
    will apply for: X-PAR for male proband, or any
    proband lacking a paternal sample for whatever reason
    """
    input:
        child_bam="results/bams/{projectid0}/PMGRC-{childid}-{childid}-0.bam",
        parent2_bam=lambda wildcards: tc.get_subjects_by_family(
            wildcards,
            checkpoints,
            wildcards.childid,
            2,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams/",
            ".bam",
        ),
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deeptrio/split_ranges/{splitnum}.ssv",
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    output:
        examples=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/mother_only/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                relation=["child", "parent2"],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        gvcfs=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/mother_only/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
                relation=["child", "parent2"],
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        jsons=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/mother_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz.example_info.json",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_make_examples/mother_only/{projectid0}/PMGRC-{childid}-{childid}-0.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/mother_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        gvcf_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/mother_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        tmpdir="/tmp",
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["make_examples_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
        tmpdir="/tmp",
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "mkdir -p {params.tmpdir} && '
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir {params.tmpdir} "
        "make_examples --mode calling "
        "--ref {input.fasta} "
        "--reads {input.child_bam} --reads_parent2 {input.parent2_bam} "
        '--regions \\"$(cat {input.intervals})\\" '
        "--examples {params.shard_string} --channels insert_size "
        "--gvcf {params.gvcf_string} "
        '--task {{}}"'


rule deeptrio_make_examples_father_only:
    """
    Run deeptrio make_examples in a hybrid
    embarrassingly parallel fashion. This rule
    will apply for: Y-PAR for male proband, or any
    proband lacking a maternal sample for whatever reason
    """
    input:
        child_bam="results/bams/{projectid0}/PMGRC-{childid}-{childid}-0.bam",
        parent1_bam=lambda wildcards: tc.get_subjects_by_family(
            wildcards,
            checkpoints,
            wildcards.childid,
            1,
            bam_manifest["projectid"],
            bam_manifest["sampleid"],
            "results/bams/",
            ".bam",
        ),
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        intervals="results/deeptrio/split_ranges/{splitnum}.ssv",
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    output:
        examples=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/father_only/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                relation=["child", "parent1"],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        gvcfs=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/father_only/PMGRC-{{childid}}-{{childid}}-0_{relation}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
                relation=["child", "parent1"],
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
        jsons=temp(
            expand(
                "results/deeptrio/{{projectid0}}/make_examples/father_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz.example_info.json",
                shardnum=[
                    str(i).rjust(5, "0")
                    for i in range(config_resources["deeptrio"]["threads"])
                ],
                shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
            )
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_make_examples/father_only/{projectid0}/PMGRC-{childid}-{childid}-0.{splitnum}.tsv"
    params:
        shard_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/father_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        gvcf_string=expand(
            "results/deeptrio/{{projectid0}}/make_examples/father_only/PMGRC-{{childid}}-{{childid}}-0.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        tmpdir="/tmp",
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["make_examples_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
        tmpdir="/tmp",
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "mkdir -p {params.tmpdir} && '
        "seq 0 $(({threads}-1)) | parallel -j{threads} --tmpdir {params.tmpdir} "
        "make_examples --mode calling "
        "--ref {input.fasta} "
        "--reads {input.child_bam} --reads_parent1 {input.parent1_bam} "
        '--regions \\"$(cat {input.intervals})\\" '
        "--examples {params.shard_string} --channels insert_size "
        "--gvcf {params.gvcf_string} "
        '--task {{}}"'


rule deeptrio_call_variants:
    """
    Run deeptrio call_variants in an
    embarrassingly parallel fashion.
    """
    input:
        gz=lambda wildcards: expand(
            "results/deeptrio/{{projectid}}/make_examples/{trio_structure}/{{sampleid}}_{{relation}}.{{splitnum}}.tfrecord-{shardnum}-of-{shardmax}.gz",
            trio_structure=tc.determine_trio_structure(
                wildcards,
                checkpoints,
                config,
                bam_manifest,
                wildcards.sampleid,
                wildcards.splitnum,
            ),
            shardnum=[
                str(i).rjust(5, "0")
                for i in range(config_resources["deeptrio"]["threads"])
            ],
            shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
        ),
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    output:
        gz=temp(
            "results/deeptrio/{projectid}/call_variants/{sampleid}_{relation,child|parent1|parent2}.{splitnum}.tfrecord.gz",
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_call_variants/{projectid}/{sampleid}.{splitnum}.{relation}.tsv"
    params:
        shard_string=lambda wildcards: expand(
            "results/deeptrio/{{projectid}}/make_examples/{trio_structure}/{{sampleid}}_{{relation}}.{{splitnum}}.tfrecord@{shardmax}.gz",
            trio_structure=tc.determine_trio_structure(
                wildcards,
                checkpoints,
                config,
                bam_manifest,
                wildcards.sampleid,
                wildcards.splitnum,
            ),
            shardmax=config_resources["deeptrio"]["threads"],
        ),
        docker_model=lambda wildcards: "/opt/models/deeptrio/wgs/{}/model.ckpt".format(
            "child" if wildcards.relation == "child" else "parent"
        ),
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: config_resources["deeptrio"]["threads"]
    resources:
        mem_mb=config_resources["deeptrio"]["call_variants_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "'
        "call_variants "
        "--outfile {output.gz} "
        "--examples {params.shard_string} "
        '--checkpoint \\"{params.docker_model}\\""'


rule deeptrio_postprocess_variants:
    """
    Run deeptrio postprocess_variants in an
    embarrassingly parallel fashion.
    """
    input:
        gz="results/deeptrio/{projectid}/call_variants/{sampleid}_{relation}.{splitnum}.tfrecord.gz",
        fasta="reference_data/bwa/{}/ref.fasta".format(reference_build),
        gvcf=lambda wildcards: expand(
            "results/deeptrio/{{projectid}}/make_examples/{trio_structure}/{{sampleid}}_{{relation}}.{{splitnum}}.gvcf.tfrecord-{shardnum}-of-{shardmax}.gz",
            trio_structure=tc.determine_trio_structure(
                wildcards,
                checkpoints,
                config,
                bam_manifest,
                wildcards.sampleid,
                wildcards.splitnum,
            ),
            shardnum=[
                str(i).rjust(5, "0")
                for i in range(config_resources["deeptrio"]["threads"])
            ],
            shardmax=str(config_resources["deeptrio"]["threads"]).rjust(5, "0"),
        ),
        fai="reference_data/bwa/{}/ref.fasta.fai".format(reference_build),
        sif="results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    output:
        vcf=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}_{relation,child|parent1|parent2}.{splitnum}.vcf.gz",
        ),
        gvcf=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}_{relation,child|parent1|parent2}.{splitnum}.g.vcf.gz",
        ),
        tbi=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}_{relation,child|parent1|parent2}.{splitnum}.vcf.gz.tbi",
        ),
        html=temp(
            "results/deeptrio/{projectid}/postprocess_variants/{sampleid}_{relation,child|parent1|parent2}.{splitnum}.visual_report.html"
        ),
    params:
        gvcf_string=lambda wildcards: expand(
            "results/deeptrio/{{projectid}}/make_examples/{trio_structure}/{{sampleid}}_{{relation}}.{{splitnum}}.gvcf.tfrecord@{shardmax}.gz",
            trio_structure=tc.determine_trio_structure(
                wildcards,
                checkpoints,
                config,
                bam_manifest,
                wildcards.sampleid,
                wildcards.splitnum,
            ),
            shardmax=config_resources["deeptrio"]["threads"],
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_postprocess_variants/{projectid}/{sampleid}.{relation}.{splitnum}.tsv"
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: 1
    resources:
        mem_mb=config_resources["deeptrio"]["postprocess_variants_memory"],
        qname=rc.select_queue(
            config_resources["deeptrio"]["queue"], config_resources["queues"]
        ),
    shell:
        'apptainer exec -B /usr/lib/locale/:/usr/lib/locale/ {input.sif} sh -c "'
        "postprocess_variants "
        "--ref {input.fasta} "
        "--infile {input.gz} "
        "--nonvariant_site_tfrecord_path {params.gvcf_string} "
        "--gvcf_outfile {output.gvcf} "
        '--outfile {output.vcf}"'


rule deeptrio_combine_regions:
    """
    Combine per-region deeptrio vcfs
    into a single mega vcf.
    """
    input:
        expand(
            "results/deeptrio/{{projectid}}/postprocess_variants/{{sampleid}}_{{relation}}.{splitnum}.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        "results/deeptrio/{projectid}/{sampleid}_{relation,child|parent1|parent2}.sorted.vcf.gz",
    benchmark:
        "results/performance_benchmarks/deeptrio_combine_regions/{projectid}/{sampleid}_{relation}.tsv"
    conda:
        "../envs/bcftools.yaml" if not use_containers else None
    container:
        "{}/bcftools.sif".format(apptainer_images) if use_containers else None
    threads: config_resources["bcftools"]["threads"]
    resources:
        mem_mb=config_resources["bcftools"]["memory"],
        qname=rc.select_queue(
            config_resources["bcftools"]["queue"], config_resources["queues"]
        ),
    shell:
        "bcftools concat --threads {threads} -O z -o {output} {input}"


use rule deeptrio_combine_regions as deeptrio_combine_gvcfs with:
    input:
        expand(
            "results/deeptrio/{{projectid}}/postprocess_variants/{{sampleid}}_{{relation}}.{splitnum}.g.vcf.gz",
            splitnum=[i + 1 for i in range(caller_num_intervals)],
        ),
    output:
        temp(
            "results/deeptrio/{projectid}/{sampleid}_{relation,child|parent1|parent2}.sorted.g.vcf.gz"
        ),
    benchmark:
        "results/performance_benchmarks/deeptrio_combine_gvcfs/{projectid}/{sampleid}_{relation}.tsv"


rule deeptrio_rename_vcf_outputs:
    """
    Rename the deeptrio-internal "parent1/parent2/child" labels to more reasonable
    things that the rest of the logic understands
    """
    input:
        lambda wildcards: "results/deeptrio/"
        + "{}.sorted.{{suffix}}".format(
            tc.get_subjects_by_family(
                wildcards,
                checkpoints,
                wildcards.probandid,
                0,
                bam_manifest["projectid"],
                bam_manifest["sampleid"],
                "",
                "_parent{}".format(wildcards.relcode)
                if wildcards.sampleid != wildcards.probandid
                else "_child",
            )[0]
        ),
    output:
        "results/deeptrio/{projectid}/PMGRC-{sampleid}-{probandid}-{relcode,[0-9]}.sorted.{suffix}",
    shell:
        "cp {input} {output}"


rule rtg_create_sdf:
    """
    Convert a fasta to an sdf format *folder* for rtg tools' particularities
    """
    input:
        fasta="reference_data/bwa/{genome}/ref.fasta",
    output:
        directory("results/{genome}/ref.fasta.sdf"),
    benchmark:
        "results/performance_benchmarks/create_sdf/{genome}.tsv"
    conda:
        "../envs/vcfeval.yaml"
    threads: 1
    resources:
        qname="small",
        mem_mb=16000,
    shell:
        "rtg RTG_MEM=12G format -f fasta -o {output} {input}"


use rule somalier_build_pedfile as rtg_create_cluster_pedigree with:
    output:
        ped="results/deeptrio/{subset}.ped",
        problems=temp("results/deeptrio/{subset}.discordant_annotations.tsv"),
    benchmark:
        "results/performance_benchmarks/rtg_create_cluster_pedigree/{subset}.tsv"


rule rtg_annotate_vcf:
    input:
        vcf="results/glnexus/{family_cluster}/merged_callset.filtered.regions.vcf.gz",
        ped="results/deeptrio/{family_cluster}.ped",
        sdf="results/{}/ref.fasta.sdf".format(reference_build),
    output:
        vcf="results/deeptrio/{family_cluster}.annotated.vcf.gz",
    conda:
        "../envs/vcfeval.yaml"
    threads: 1
    resources:
        qname="small",
        mem_mb=16000,
    shell:
        "rtg RTG_MEM=12G mendelian -i {input.vcf} -o {output.vcf} --pedigree {input.ped} -t {input.sdf}"


rule aggregate_deeptrio_output:
    """
    Dispatch "deeptrio" runs for any proband who also has
    at least one parent present
    """
    input:
        lambda wildcards: tc.get_probands_with_structure(checkpoints),
    output:
        temp("results/deeptrio/.deeptrio_calls_split"),
    shell:
        "touch {output}"

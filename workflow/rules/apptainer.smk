rule apptainer_deeptrio:
    """
    Because snakemake's internal container support doesn't
    mesh with the deeptrio image, pull the image and create
    a sif image with apptainer.
    """
    output:
        "results/apptainer_images/deepvariant_{}.sif".format(
            config["deeptrio"]["docker-version"]
        ),
    params:
        outdir="results/apptainer_images",
        image_version=config["deeptrio"]["docker-version"],
    conda:
        "../envs/apptainer.yaml" if not use_containers else None
    threads: 1
    resources:
        mem_mb=2000,
        qname="small",
    shell:
        "apptainer pull --dir {params.outdir} docker://google/deepvariant:{params.image_version}"

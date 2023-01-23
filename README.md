# Snakemake workflow: WGS Downstream Analyses

This workflow is intended to be the R&D space for PMGRC WGS downstream (post-alignment or calling) analyses. Target analyses are as follows:

- Cross-flowcell QC with [somalier](https://github.com/brentp/somalier)
- [ExpansionHunter Denovo](https://github.com/Illumina/ExpansionHunterDenovo)
- Joint calling with [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- others TBD

New global targets should be added in `workflow/Snakefile`. Content in `workflow/Snakefile` and the snakefiles in `workflow/rules` should be specifically _rules_; python infrastructure should be composed as subroutines under `lib/` and constructed in such a manner as to be testable with [pytest](https://docs.pytest.org/en/7.2.x/). Rules can call embedded scripts (in python or R/Rmd) from `workflow/scripts`; again, these should be constructed to be testable with pytest or [testthat](https://testthat.r-lib.org/).

## Authors

* Lightning Auriga (@lightning.auriga)

## Usage

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone git@gitlab.com:lightning.auriga1/wgs-pipeline.git
```

Note that this requires local git ssh key configuration; see [here](https://docs.gitlab.com/ee/user/ssh.html) for instructions as required.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `manifest.tsv` to specify your sample setup.

The following settings are recognized in `config/config.yaml`.

- `manifest`: relative path to run manifest
- `multiqc-config`: relative path to configuration settings for cross-flowcell multiQC report
- `genome-build`: requested genome reference build to use for this analysis run. this should match the tags used in the reference data blocks below.
- `behaviors`: user-configurable modifiers to how the pipeline will run
  - `symlink-bams`: whether to copy (no) or symlink (yes) input bams into workspace. symlinking is faster and more memory-efficient, but
    less reproducible, as the upstream files may vanish leaving no way to regenerate your analysis from scratch.
- `parameters`: tool-specific parameters. note that this section is a work in progress, somewhat more than the rest
- `references`: human genome reference data applicable to multiple tools

The following columns are expected in the run manifest, by default at `config/manifest.tsv`:
- `projectid`: run ID, or other desired grouping of sequencing samples. this will be a subdirectory under individual tools in `results/`
- `sampleid`: sequencing ID for sample
- `bam`: absolute or relative path to sample alignment (post-BQSR) bam file

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --profile sge-profile --cluster-config config/cluster.yaml --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Cluster profiles

Snakemake interfaces with job schedulers via _cluster profiles_. For running jobs on SGE, you can use
the cookiecutter template [here](https://github.com/Snakemake-Profiles/sge).



### Step 5: Investigate results

#### Read Quality Control

Quality control data from fastqc and fastp are aggregated in a multiqc report at `results/multiqc/{projectid}/multiqc.fastq.html`.
This version of the quality control report is split by lane within the report, if per-lane fastqs have been provided
and annotated in the manifest.

#### Alignment

Quality control data from the above read QC tools as well as somalier, verifybamid, alignstats,  and assorted gatk4/picard analysis tools
are aggregated in a multiqc report at `results/multiqc/{projectid}/multiqc.alignment.html`. This version of the quality report is currently
under active modification, and should only be used for considering alignment QC results (as opposed to pre-alignment read QC) until further notice.

#### Variant Calling

All variant calling in this pipeline is conducted per-subject, ignoring batch data. SNV calls from the user-configured tool
(e.g. DeepVariant, Octopus) are aggregated in `results/{toolname}/{projectid}/{sampleid}.sorted.vcf.gz`. SV calls from
ensemble calling based on user-configured tools and exclusion criteria are aggregated in `results/final/{projectid}/{sampleid}.sv.vcf.gz`.
Note that these paths and filenames are subject to change before stabilization of the workflow.

#### (DeepVariant only) GVCFs for Batch Calling

If the user has selected DeepVariant for SNV calling, gvcf files per-subject will be collected at
`results/deepvariant/{projectid}/{sampleid}.g.vcf.gz`. These files are not used in the pipeline itself,
and represent the raw output of `deepvariant postprocess_variants`. At some point, these will likely receive
further processing in anticipation of use with e.g. GLnexus in a different pipeline. Also at some point,
I'll probably add a flag for disabling the production of gvcf output, but that's not urgent.

gvcf output is not supported by Octopus and so is not possible in this pipeline.

#### Optional: emit methods description and software version data

Some users may be interested in a specific breakdown of workflow methods, relevant software versions pulled from
conda, and the effects of certain important user configuration settings. This information can be generated upon
completion of a workflow with the command `snakemake -j1 summarize_methods`. This will create an additional output
file, `results/reports/methods_summary.md`, that contains the best description of the workflow as was actually run.
Note that this information focuses on methodology as opposed to results, and only requires the relevant conda
environments be installed; so if you want to predict what the workflow will do before actually running it,
complete user configuration, run `snakemake -j1 --use-conda --conda-create-envs-only`, and then run the `summarize_methods`
target to generate a markdown description of what _would_ happen if the pipeline were deployed.

#### Data Release

When the pipeline is run in `release` mode, postprocessed output files for each flowcell will be emitted under
`results/export/{flowcell_id}`. The following files will be present, with modifications as annotated:

- aligned reads, represented as lossless crams
  - the source reference fasta used for these files is both in the cram header as a `@CO` and also
    annotated in the output methods html
- `crai` index files for the above cram files
- SNV vcfs, bgzip-compressed, with the following modifications (derived from [Pedersen _et al._](https://doi.org/10.1038/s41525-021-00227-3):
  - only FILTER=PASS variants
  - multiallelics split to biallelics
  - GQ >= 20
  - DP >= 10
  - for heterozygotes, allele balance on `[0.2, 0.8]`
  - for homozygous alts, allele balance `<= 0.04`
  - variants intersected with configured exclusion regions removed
- tabix indices for the above vcfs
- md5 sums for all above files
- a plaintext `manifest.tsv` containing a list of the above data files
- an immutable `methods_summary.html` representing the rendered version of the above methods description for the version of the pipeline that created the release
- SVs TBD

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@gitlab.com:lightning.auriga1/wgs-pipeline.git` or `git remote add -f upstream https://github.com/snakemake-workflows/wgs-pipeline.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/default workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/default config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository. This project follows git flow; feature branches off of dev are welcome.

1. [Clone](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html) the fork to your local system, to a different place than where you ran your analysis.
2. Check out a branch off of dev:
```
git fetch
git checkout dev
git checkout -b your-new-branch
```
3. Make whatever changes best please you to your feature branch.
4. Commit and push your changes to your branch.
5. Create a [merge request](https://docs.gitlab.com/ee/user/project/merge_requests/) against dev.

## Testing

Testing infrastructure for embedded python and R scripts is installed under `lib/` and `workflow/scripts/`. Additional testing
coverage for the Snakemake infrastructure itself should be added once the workflow is more mature ([see here](https://github.com/lightning-auriga/snakemake-unit-tests)).

### Python testing with `pytest`
The testing under `lib/` is currently functional. Partial testing exists for the builtin scripts under `workflow/scripts`: the new utilities
for this implementation are tested, but some code inherited from the legacy pipeline(s) is not yet covered. To run the tests, do the following (from top level):

```bash
mamba install pytest-cov
pytest --cov=lib --cov=workflow/scripts lib workflow/scripts
```


### R testing with `testthat`
The testing under `workflow/scripts` is currently functional. The tests can be run with the utility script `run_tests.R`:

```bash
Rscript ./run_tests.R
```

To execute the above command, the environment must have an instance of R with appropriate libraries:

```bash
mamba install -c bioconda -c conda-forge "r-base>=4" r-testthat r-covr r-r.utils r-desctools
```

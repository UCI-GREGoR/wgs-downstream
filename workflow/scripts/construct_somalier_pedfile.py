import pandas as pd


def convert_sex_representation(sample_sex: str) -> int:
    """
    Take a self-reported sex representation written out as
    an English word, and convert to a plink-style integer
    representation. 0 -> Unknown, 2 -> Female, 1 -> Male
    """
    if sample_sex.lower() == "female":
        return 2
    elif sample_sex.lower() == "male":
        return 1
    else:
        return 0


def add_problem(problems: dict, sampleid: str, problem: str) -> dict:
    """
    Add annotation of problematic ID/annotation configuration to
    aggregator for later reporting
    """
    if sampleid in problems.keys():
        problems[sampleid] = problems[sampleid] + "; " + problem
    else:
        problems[sampleid] = problem
    return problems


def run_construct_somalier_pedfile(
    linker: str,
    projectids: list,
    sampleids: list,
    valid_pmgrcids: list,
    outfn: str,
    problemfn: str,
) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder 0 for unknown.
    """

    ## strip ruids
    valid_pmgrcids = [x.split("/")[1] for x in valid_pmgrcids]

    ## load linker information formatted from the lab logbook
    linker_data = pd.read_csv(linker, sep="\t")

    ## Grab pmgrcids from linker
    aligned_pmgrcids = []
    for projectid, sampleid in zip(projectids, sampleids):
        res = linker_data.loc[
            (linker_data["ru"] == projectid) & (linker_data["sq"] == sampleid), "pmgrc"
        ]
        if len(res) == 0:
            aligned_pmgrcids.append(projectid + "_" + sampleid)
        else:
            aligned_pmgrcids.append(res.to_list()[0])

    ## Do not tolerate duplicates. This shouldn't actually happen, but hey.
    ids = pd.DataFrame(
        data={"ruid": projectids, "sampleid": sampleids, "pmgrcid": aligned_pmgrcids}
    ).drop_duplicates()

    ## track detected problems with ID consistency
    problems = {}

    ## iterate across manifest subjects and try to find matching values
    self_reported_sex = []
    family_id = []
    mat_id = []
    pat_id = []
    parent_data = {}
    for pmgrcid, ruid, sqid in zip(
        linker_data["pmgrc"], linker_data["ru"], linker_data["sq"]
    ):
        if not (pmgrcid in valid_pmgrcids):
            parent_data[str(ruid) + "_" + str(sqid)] = str(ruid) + "_" + str(sqid)
        else:
            parsed_sample_id = pmgrcid.split("-")
            parent_data[
                "{}-{}".format(parsed_sample_id[2], parsed_sample_id[3])
            ] = pmgrcid

    for ruid, sampleid, pmgrc_id in zip(ids["ruid"], ids["sampleid"], ids["pmgrcid"]):
        sample_sex = linker_data.loc[
            (linker_data["ru"] == ruid) & (linker_data["sq"] == sampleid), "sex"
        ]
        if len(sample_sex) == 1:
            parsed_sample_id = pmgrc_id.split("-")
            invalid_family_structure = False
            if (
                parsed_sample_id[1] == parsed_sample_id[2]
                and parsed_sample_id[3] != "0"
            ):
                problems = add_problem(
                    problems, sampleid, "subject with proband ID is not flagged -0"
                )
                invalid_family_structure = True
            sex_representation = convert_sex_representation(sample_sex.to_list()[0])
            invalid_sex_configuration = False
            if sex_representation != 0:
                if parsed_sample_id[3] == "1" and sex_representation != 1:
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband father without self-reported male sex",
                    )
                    invalid_sex_configuration = True
                elif parsed_sample_id[3] == "2" and sex_representation != 2:
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband mother without self-reported female sex",
                    )
                    invalid_sex_configuration = True
                if not invalid_sex_configuration:
                    self_reported_sex.append(sex_representation)
                else:
                    self_reported_sex.append(0)
            elif parsed_sample_id[3] == "1":
                self_reported_sex.append(convert_sex_representation("Male"))
            elif parsed_sample_id[3] == "2":
                self_reported_sex.append(convert_sex_representation("Female"))
            else:
                self_reported_sex.append(0)
            if (
                "{}-1".format(parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                pat_id.append(parent_data["{}-1".format(parsed_sample_id[1])])
            else:
                pat_id.append("0")
            if (
                "{}-2".format(parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                mat_id.append(parent_data["{}-2".format(parsed_sample_id[1])])
            else:
                mat_id.append("0")
            family_id.append(parsed_sample_id[2])
        else:
            self_reported_sex.append(0)
            mat_id.append(0)
            pat_id.append(0)
            family_id.append(sampleid)
            problems = add_problem(
                problems, sampleid, "self-reported sex missing from annotations"
            )

    x = pd.DataFrame(
        data={
            "FID": family_id,
            "Sample": ids["pmgrcid"],
            "Pat": pat_id,
            "Mat": mat_id,
            "Sex": self_reported_sex,
            "Pheno": ["-9" for x in ids["sampleid"]],
        }
    )
    x.to_csv(outfn, sep="\t", index=False, header=False)
    for problem_key in problems.keys():
        problems[problem_key] = [problems[problem_key]]
    problems = pd.DataFrame(data=problems)
    problems.to_csv(problemfn, sep="\t", index=False, header=True)


run_construct_somalier_pedfile(
    snakemake.input[0],  # noqa: F821
    snakemake.params["projectids"],  # noqa: F821
    snakemake.params["subjectids"],  # noqa: F821
    snakemake.params["valid_pmgrcids"],  # noqa: F821
    snakemake.output["ped"],  # noqa: F821
    snakemake.output["problems"],  # noqa: F821
)

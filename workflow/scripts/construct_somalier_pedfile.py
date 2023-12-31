import re

import numpy as np
import pandas as pd


def convert_sex_representation(sample_sex: str) -> str:
    """
    Take a self-reported sex representation written out as
    an English word, and convert to a plink-style integer (as string)
    representation. 0 -> Unknown, 2 -> Female, 1 -> Male
    """
    if sample_sex.lower() == "female":
        return "2"
    elif sample_sex.lower() == "male":
        return "1"
    else:
        return "0"


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


def construct_final_id(sampleid: str, use_somalier_id: bool) -> str:
    """
    Adjust sample ID format for use in either somalier or, really, everything else
    """
    if use_somalier_id:
        res = sampleid
    else:
        id_parse = re.match("^(.*)_SQ[0-9]{4}$", sampleid)
        res = id_parse[1] if id_parse else sampleid
    ## catch corner case malformation of missing parent
    if res == "0_":
        res = "0"
    return res


def add_affected_status(df: pd.DataFrame, affected_status: str) -> pd.DataFrame:
    """
    Suck in user annotations of sample affected status and add to results data
    frame column "Pheno"
    """
    affected = pd.read_table(affected_status, sep="\t").set_index("participant_id")
    df = df.join([affected])
    df["Pheno"] = df["affected_status"].map(
        {"Affected": "2", "Unaffected": "1", np.nan: "1"}
    )
    return df[df.columns[range(6)]]


def run_construct_somalier_pedfile(
    sex_manifest_fn: str,
    affected_status: str,
    sampleids: list,
    valid_subjectids: list,
    use_somalier_ids: bool,
    outfn: str,
    problemfn: str,
) -> None:
    """
    Create a somalier format pedfile for linking sample id to sex annotation.

    Currently, this populates the sex entry for all subjects with placeholder 0 for unknown.
    """

    ## track detected problems with ID consistency
    problems = {}

    ## iterate across manifest subjects and try to find matching values
    self_reported_sex = []
    family_id = []
    mat_id = []
    pat_id = []
    parent_data = {}
    for sampleid in sampleids:
        if not (sampleid in valid_subjectids):
            parent_data[sampleid] = sampleid
        else:
            parsed_sample_id = sampleid.split("-")
            parent_data[
                "{}-{}".format(parsed_sample_id[2], parsed_sample_id[3])
            ] = sampleid

    sex_manifest = pd.read_csv(sex_manifest_fn, sep="\t").set_index(
        "sampleid", drop=False
    )

    for sampleid in sampleids:
        if sampleid in sex_manifest.index:
            sample_sex = sex_manifest.loc[sampleid, "sex"]
            parsed_sample_id = sampleid.split("-")
            is_verifiable = len(parsed_sample_id) == 4
            invalid_family_structure = False
            if (
                is_verifiable
                and parsed_sample_id[1] == parsed_sample_id[2]
                and parsed_sample_id[3] != "0"
            ):
                problems = add_problem(
                    problems, sampleid, "subject with proband ID is not flagged -0"
                )
                invalid_family_structure = True
            sex_representation = convert_sex_representation(sample_sex)
            invalid_sex_configuration = False
            if sex_representation != "0":
                if (
                    is_verifiable
                    and parsed_sample_id[3] == "1"
                    and sex_representation != "1"
                ):
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband father without self-reported male sex",
                    )
                    invalid_sex_configuration = True
                elif (
                    is_verifiable
                    and parsed_sample_id[3] == "2"
                    and sex_representation != "2"
                ):
                    problems = add_problem(
                        problems,
                        sampleid,
                        "subject is proband mother without self-reported female sex",
                    )
                    invalid_sex_configuration = True
                if not invalid_sex_configuration:
                    self_reported_sex.append(sex_representation)
                else:
                    self_reported_sex.append("0")
            elif is_verifiable and parsed_sample_id[3] == "1":
                self_reported_sex.append(convert_sex_representation("Male"))
            elif is_verifiable and parsed_sample_id[3] == "2":
                self_reported_sex.append(convert_sex_representation("Female"))
            else:
                self_reported_sex.append("0")
            if (
                is_verifiable
                and "{}-1".format(parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                pat_id.append(parent_data["{}-1".format(parsed_sample_id[1])])
            else:
                pat_id.append("0")
            if (
                is_verifiable
                and "{}-2".format(parsed_sample_id[1]) in parent_data
                and not invalid_family_structure
            ):
                mat_id.append(parent_data["{}-2".format(parsed_sample_id[1])])
            else:
                mat_id.append("0")
            if is_verifiable:
                family_id.append(parsed_sample_id[2])
            else:
                family_id.append(sampleid)
        else:
            self_reported_sex.append("0")
            mat_id.append("0")
            pat_id.append("0")
            family_id.append(sampleid)
            problems = add_problem(
                problems, sampleid, "self-reported sex missing from annotations"
            )

    x = (
        pd.DataFrame(
            data={
                "FID": family_id,
                "Sample": [
                    construct_final_id(sampleid, use_somalier_ids)
                    for sampleid in sampleids
                ],
                "Pat": [construct_final_id(str(y), use_somalier_ids) for y in pat_id],
                "Mat": [construct_final_id(str(y), use_somalier_ids) for y in mat_id],
                "Sex": self_reported_sex,
                "Pheno": ["1" for x in sampleids],
            }
        )
        .sort_values(by=["FID", "Sample"])
        .set_index("Sample", drop=False)
    )
    if affected_status is not None:
        x = add_affected_status(x, affected_status)
    x.to_csv(outfn, sep="\t", index=False, header=False)

    problems = {"sampleid": problems.keys(), "problem": problems.values()}
    problems = pd.DataFrame(data=problems)
    problems.to_csv(problemfn, sep="\t", index=False, header=True)


run_construct_somalier_pedfile(
    snakemake.input["sex_manifest"],  # noqa: F821
    snakemake.input["affected_status"]  # noqa: F821
    if "affected_status" in snakemake.input.keys()  # noqa: F821
    else None,
    snakemake.params["sampleids"],  # noqa: F821
    snakemake.params["valid_subjectids"],  # noqa: F821
    snakemake.params["use_somalier_ids"],  # noqa: F821
    snakemake.output["ped"],  # noqa: F821
    snakemake.output["problems"],  # noqa: F821
)

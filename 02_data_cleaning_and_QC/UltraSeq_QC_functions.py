#!/usr/bin/env python
# coding: utf-8
"""QC metrics for UltraSeq / immunoediting analysis (cohort + sample-specific)."""

import math

import numpy as np
import pandas as pd


def Find_Controls(input_gRNA_df, input_pattern):
    """Return gRNAs whose Targeted_gene_name matches the regex pattern."""
    return input_gRNA_df.loc[
        input_gRNA_df["Targeted_gene_name"].str.contains(input_pattern, na=False, regex=True),
        "gRNA",
    ].unique()


def Generate_ref_input_df(input_df, input_sample_list, input_cell_cutoff):
    return input_df[
        (input_df["Cell_number"] > input_cell_cutoff) & (input_df["Sample_ID"].isin(input_sample_list))
    ]


def LN_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    temp_var = log_vector.var()
    if len(log_vector) == 1:
        temp_var = 0
    return math.exp(temp_mean + 0.5 * temp_var)


def Geometric_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    return math.exp(temp_mean)


def Cal_Tumor_Size_simple(x, input_percentile):
    d = {}
    temp_vect = x["Cell_number"]
    if isinstance(temp_vect, int):
        temp_vect = [temp_vect]
    d["LN_mean"] = LN_Mean(temp_vect)
    d["Geo_mean"] = Geometric_Mean(temp_vect)
    Percentile_list = list(np.percentile(temp_vect, input_percentile))
    for c, y in enumerate(input_percentile):
        temp_name = str(y) + "_percentile"
        d[temp_name] = Percentile_list[c]
    d["TTN"] = len(temp_vect)
    d["TTB"] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys()))


def Generate_Normalized_Metrics(input_df1, input_df2, trait_list):
    """Normalize input_df1 by input_df2 (per-gRNA ratio of cohort metrics)."""
    temp1 = input_df1.set_index("gRNA")
    temp2 = input_df2.set_index("gRNA").loc[temp1.index]
    temp_output_df = pd.DataFrame({"gRNA": temp1.index.values})
    for temp_cname in trait_list:
        temp_cname_new = temp_cname + "_normalized"
        temp_output_df[temp_cname_new] = np.array(temp1[temp_cname].to_list()) / np.array(
            temp2[temp_cname].to_list()
        )
    return temp_output_df


def Add_Corhort_Specific_Relative_Metrics(input_df, input_control_list):
    temp_sub = input_df[input_df["gRNA"].isin(input_control_list)]
    for temp_cname in input_df.drop(columns=["gRNA"], inplace=False).columns:
        temp_name = temp_cname + "_relative"
        input_df[temp_name] = input_df[temp_cname] / temp_sub[temp_cname].median()


def Calculate_Relative_Normalized_Metrics(input_df1, input_df2, percentile_list, input_control_gRNA_list):
    temp_df = input_df1.groupby(["gRNA"], as_index=False).apply(
        lambda g: Cal_Tumor_Size_simple(g, percentile_list), include_groups=False
    )
    temp_df2 = input_df2.groupby(["gRNA"], as_index=False).apply(
        lambda g: Cal_Tumor_Size_simple(g, percentile_list), include_groups=False
    )
    temp_out = Generate_Normalized_Metrics(temp_df, temp_df2, ["TTN", "TTB"])
    temp_df = temp_df.merge(temp_out, on="gRNA")
    Add_Corhort_Specific_Relative_Metrics(temp_df, input_control_gRNA_list)
    temp_df["Type"] = temp_df.apply(
        lambda x: "Inert" if (x["gRNA"] in input_control_gRNA_list) else "Experiment", axis=1
    )
    temp_df = temp_df.merge(
        input_df1[["gRNA", "Targeted_gene_name", "Identity", "Numbered_gene_name"]].drop_duplicates(),
        how="inner",
        on="gRNA",
    )
    return temp_df


def Shannon_Index(input_vector):
    temp_fraction_list = np.array(input_vector) / sum(input_vector)
    return -sum(temp_fraction_list * np.log(temp_fraction_list))


def Simpson_Index(input_vector):
    temp_fraction_list = np.array(input_vector) / sum(input_vector)
    return 1 / sum(temp_fraction_list**2)


def Cal_Tumor_Size(x, input_percentile):
    d = {}
    temp_vect = x["Cell_number"]
    if isinstance(temp_vect, int):
        temp_vect = [temp_vect]
    d["LN_mean"] = LN_Mean(temp_vect)
    d["Geo_mean"] = Geometric_Mean(temp_vect)
    Percentile_list = list(np.percentile(temp_vect, input_percentile))
    for c, y in enumerate(input_percentile):
        temp_name = str(y) + "_percentile"
        d[temp_name] = Percentile_list[c]
    d["TTB"] = sum(temp_vect)
    d["TTN"] = len(temp_vect)
    d["Shannon_diversity"] = Shannon_Index(temp_vect)
    d["Simpson_diversity"] = Simpson_Index(temp_vect)
    return pd.Series(d, index=list(d.keys()))


def Generate_Normalized_Metrics_by_KT(input_df1, input_df2):
    """Normalize sample-level metrics using KT mice (mean profile over gRNAs)."""
    temp1 = input_df1.set_index("gRNA")
    kt_mean = input_df2.groupby("gRNA")[["TTB", "TTN"]].mean()
    temp2 = (kt_mean / kt_mean.sum()).loc[temp1.index]
    temp_output_df = pd.DataFrame({"gRNA": temp1.index.values, "Sample_ID": temp1.Sample_ID.values})
    dict1 = {sid: temp1[temp1.Sample_ID == sid]["TTN"].sum() for sid in temp1.Sample_ID.unique()}
    dict2 = {sid: temp1[temp1.Sample_ID == sid]["TTB"].sum() for sid in temp1.Sample_ID.unique()}
    temp_list = list(temp1.apply(lambda x: x["TTN"] / dict1.get(x["Sample_ID"]), axis=1))
    temp_output_df["TTN_normalized"] = temp_list / np.array(temp2["TTN"].to_list())
    temp_list = list(temp1.apply(lambda x: x["TTB"] / dict2.get(x["Sample_ID"]), axis=1))
    temp_output_df["TTB_normalized"] = temp_list / np.array(temp2["TTB"].to_list())
    return temp_output_df


def Add_Sample_Specific_Relative_Metrics(input_df, input_control_list, input_percentile):
    temp_sub = input_df[input_df["gRNA"].isin(input_control_list)]
    temp_df = temp_sub.groupby("Sample_ID").apply(
        lambda g: Cal_Relative_Metrics(g, input_percentile), include_groups=False
    )
    for temp_cnames in temp_df.columns:
        temp1 = temp_cnames.replace("_median", "")
        temp2 = temp_cnames.replace("_median", "_relative")
        input_df[temp2] = (
            np.array(input_df[temp1]) / temp_df.loc[input_df["Sample_ID"].to_list(), temp_cnames].to_list()
        )


def Cal_Relative_Metrics(x, input_percentile):
    d = {}
    d["LN_mean_median"] = np.median(x["LN_mean"])
    d["Geo_mean_median"] = np.median(x["Geo_mean"])
    for c, y in enumerate(input_percentile):
        temp_name = str(y) + "_percentile"
        temp_name1 = str(y) + "_percentile_median"
        d[temp_name1] = np.median(x[temp_name])
    d["TTB_normalized_median"] = np.median(x["TTB_normalized"])
    d["TTN_normalized_median"] = np.median(x["TTN_normalized"])
    d["Shannon_diversity_median"] = np.median(x["Shannon_diversity"])
    d["Simpson_diversity_median"] = np.median(x["Simpson_diversity"])
    return pd.Series(d, index=list(d.keys()))


def Generate_Sample_Specific_Metrics(input_df, input_cell_number_cutoff, input_control, input_percentil_list):
    temp_input = input_df[input_df["Cell_number"] > input_cell_number_cutoff]
    temp_df1 = temp_input[temp_input["Identity"] == "gRNA"].groupby(
        [
            "Sample_ID",
            "gRNA",
            "Mouse_genotype",
            "Sex",
            "Targeted_gene_name",
            "Numbered_gene_name",
            "Time_after_tumor_initiation",
            "Total_lung_weight",
            "Virus_titer",
        ],
        as_index=False,
    ).apply(lambda g: Cal_Tumor_Size(g, input_percentil_list), include_groups=False)
    temp_out = Generate_Normalized_Metrics_by_KT(temp_df1, temp_df1[temp_df1.Mouse_genotype == "KT"])
    temp_out = temp_df1.merge(temp_out, on=["gRNA", "Sample_ID"])
    Add_Sample_Specific_Relative_Metrics(temp_out, input_control, input_percentil_list)
    temp_out["Type"] = temp_out.apply(
        lambda x: "Negative_control" if (x["gRNA"] in input_control) else "Experiment", axis=1
    )
    return temp_out


def Generate_sample_specific_treatment_effect(
    input_df, treatment_trait, ref_treatment, other_treatment_list, focal_metrics, experimental_info_name
):
    temp_final = pd.pivot_table(
        input_df, values=focal_metrics, index=experimental_info_name, columns="Vector_type"
    ).reset_index()
    for temp_metric in focal_metrics:
        for y in other_treatment_list:
            temp_new_name = temp_metric + "_" + y + "_standardized"
            temp_final[temp_new_name] = temp_final[temp_metric][y] / temp_final[temp_metric][ref_treatment]
    return temp_final

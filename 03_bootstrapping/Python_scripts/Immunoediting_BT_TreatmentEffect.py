#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import copy
import argparse
from scipy.stats import rankdata
# from Ultra_Seq_Boostrapping import *

## Functions
def Generate_treatment_ratio_effect_df(input_ref_df, input_treatment_df,trait_of_interest):
    ref_trait = ['Type','Targeted_gene_name','Numbered_gene_name','Bootstrap_id','gRNA']
    trait_of_interest = [x for x in trait_of_interest if x in input_ref_df.columns]
    ref_trait = list(set(input_ref_df.columns)&set(ref_trait))
    # I used how ='inner', so it is not accurate for sample with few tumor that don't appear in some bootstrapping cycles.
    temp_combined = pd.merge(input_ref_df[trait_of_interest+ref_trait],
                             input_treatment_df[trait_of_interest+ref_trait],
                            on=ref_trait, how ='inner', suffixes=('_ref', '_treatment'))
    
    # temp_combined = temp_combined.fillna(0)
    temp_new_name_list = []
    for x in trait_of_interest:
        temp = x+'_fold'
        temp_new_name_list.append(temp)
        temp_combined[temp] = temp_combined.apply(lambda y: y[x+'_treatment']/y[x + '_ref'],axis=1)

    temp_final_df = temp_combined[ref_trait+temp_new_name_list]
    return(temp_final_df)


def Generate_treatment_dif_effect_df(input_ref_df, input_treatment_df,trait_of_interest):
    ref_trait = ['Type','Targeted_gene_name','Numbered_gene_name','Bootstrap_id','gRNA']
    trait_of_interest = [x for x in trait_of_interest if x in input_ref_df.columns]
    ref_trait = list(set(input_ref_df.columns)&set(ref_trait))
    temp_combined = pd.merge(input_ref_df[trait_of_interest+ref_trait],
                             input_treatment_df[trait_of_interest+ref_trait],
                            on=ref_trait, how ='outer', suffixes=('_ref', '_treatment'))
    temp_combined = temp_combined.fillna(0)
    temp_new_name_list = []
    for x in trait_of_interest:
        temp = x+'_dif'
        temp_new_name_list.append(temp)
        temp_combined[temp] = temp_combined.apply(lambda y: y[x+'_treatment']-y[x + '_ref'],axis=1)

    temp_final_df = temp_combined[ref_trait+temp_new_name_list]
    return(temp_final_df)

def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]





def Generate_Final_Summary_Dataframe(input_df,trait_of_interest,ref_cutoff,group_trait='gRNA'):
    temp_summary = input_df[input_df['Bootstrap_id']!='Real'].groupby(group_trait,as_index = False).apply(Cal_Bootstrapping_Summary,(trait_of_interest),(ref_cutoff))
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] =='Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on = group_trait)
    # Dictionary to hold new columns
    new_columns = {}
    # Loop over traits of interest
    for temp_trait in trait_of_interest:
        temp_name0 = temp_trait + '_fraction_greater_than_ref'
        temp_name1 = temp_trait + '_pvalue'
        temp_name2 = temp_name1 + '_FDR'
        temp_name3 = temp_name1 + '_twoside'
        temp_name4 = temp_name1 + '_twoside_FDR'

        # Calculate the p-value, FDR, two-sided p-value, and FDR for two-sided
        new_columns[temp_name1] = temp_output_df.apply(lambda x: min(x[temp_name0], 1 - x[temp_name0]), axis=1)
        new_columns[temp_name2] = fdr(new_columns[temp_name1])
        new_columns[temp_name3] = new_columns[temp_name1] * 2
        new_columns[temp_name4] = fdr(new_columns[temp_name3])

    # Convert new columns into a DataFrame and concatenate with the original DataFrame
    new_columns_df = pd.DataFrame(new_columns)

    # Concatenate the new columns to the output DataFrame in one go
    temp_output_df = pd.concat([temp_output_df, new_columns_df], axis=1)
    
    # Return the updated DataFrame
    return(temp_output_df)


def Cal_Bootstrapping_Summary(x,trait_of_interest,ref_cutoff):
    d = {}
    for temp_trait in trait_of_interest:
        temp0 = temp_trait + '_95P'
        temp1 = temp_trait + '_5P'
        temp2 = temp_trait +'_fraction_greater_than_ref' # t_test pvalue column name
        temp3 = temp_trait +'_bootstrap_median'
        temp4 = temp_trait +'_bootstrap_mean'
        temp5 = temp_trait + '_97.5P'
        temp6 = temp_trait + '_2.5P'
        d[temp0] = x[temp_trait].quantile(0.95)
        d[temp1] = x[temp_trait].quantile(0.05)
        d[temp2] = sum(x[temp_trait]>ref_cutoff)/len(x[temp_trait])
        d[temp3] = x[temp_trait].mean()
        d[temp4] = x[temp_trait].median()
        d[temp5] = x[temp_trait].quantile(0.975)
        d[temp6] = x[temp_trait].quantile(0.025)
    return pd.Series(d, index=list(d.keys())) 
# -----


# -----

def main():
    parser = argparse.ArgumentParser(description='Calculate treatment effects using bootstrapping data')
    parser.add_argument("--a0", required=True, help="Path to reference data CSV")
    parser.add_argument("--a1", required=True, help="Path to treatment data CSV")
    parser.add_argument("--a2", required=True, help="Name of the group trait")
    parser.add_argument("--o1", required=True, help="Output path for summary data")

    # data input
    args = parser.parse_args()
    
    ref_address  = args.a0
    treatment_address  = args.a1
    output_address = args.o1
    group_trait = args.a2

    ref_df = pd.read_csv(ref_address) # ref data frame
    treatment_df = pd.read_csv(treatment_address) # treatment data frame

    # trait list 
    # these are the traits that are value 
    temp_trait_list1 = [x for x in ref_df.columns if ('relative' in x)&('P_' not in x)]+['P_0_percentile_relative','P_30_percentile_relative','P_50_percentile_relative','P_70_percentile_relative']
    temp_trait_list2 = [x.split('_relative')[0] for x in temp_trait_list1]    
    temp_trait_list1 = temp_trait_list1 + temp_trait_list2
    temp_ratio = Generate_treatment_ratio_effect_df(ref_df, treatment_df,temp_trait_list1)
    # Generate summary statistics
    temp_trait_new = [x for x in temp_ratio.columns if ('fold' in x)]
    Final_ratio_gRNA_summary_df = Generate_Final_Summary_Dataframe(temp_ratio,temp_trait_new,1,group_trait)
    Final_df = Final_ratio_gRNA_summary_df
    
    if 'gene' in group_trait:
        Final_df.to_csv(f"{output_address}_gene.csv", index=False)
    else:
        Final_df.to_csv(f"{output_address}_{group_trait}.csv", index=False)
    print(f"All steps finished") 

if __name__ == "__main__":
    main() 

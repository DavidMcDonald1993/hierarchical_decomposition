import os

import numpy as np 
import pandas as pd 

from scipy.stats import kruskal

def main():


    core_of_the_core = set(["pi3k", "pip3", "gab1"])

    proliferation_grow_tfs = set(["elk1", "creb", "ap1", "cmyc", 
        "p70s6_2", "hsp27"])
    apoptosis_tfs = set(["pro_apoptotic"])
    additional_genes_of_interest = set(["akt","ras"])

    output_genes = proliferation_grow_tfs.union(apoptosis_tfs)


    rank_df = pd.DataFrame()


    for output_gene in output_genes:

        print ("processing", output_gene)

        df = pd.read_csv(os.path.join("results", "{}_expressions_with_ras.csv".format(output_gene)), 
            index_col=0)
        
        original_expressions = df.loc["original"]

        modification_p_values = {}
        
        for modification in df.index[1:]:

            modification_expressions = df.loc[modification]

            _, p_value = kruskal(original_expressions, 
                modification_expressions, 
                nan_policy="omit")

            # print (gene, "\t", p_value)
            modification_p_values.update({modification: p_value})
        
        sorted_p_values = sorted(set(modification_p_values.values()))

        rank_dict = {modification: sorted_p_values.index(p_value) + 1 
            for modification, p_value in modification_p_values.items()}
        
        rank_df = rank_df.append(pd.Series(rank_dict, name=output_gene))

    mean_over_all_outputs = rank_df.mean(axis=0).to_dict()
    sorted_modifications = sorted(mean_over_all_outputs, key= mean_over_all_outputs.get)
    for modification in sorted_modifications:
        print ("{:15s}\t{}".format(modification, mean_over_all_outputs[modification]))




if __name__ == "__main__":
    main()
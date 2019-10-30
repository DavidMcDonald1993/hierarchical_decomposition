import os

import numpy as np 
import pandas as pd 

from scipy.stats import kruskal, mannwhitneyu, ttest_rel

def main():


    proliferation_grow_tfs = set(["elk1", "creb", "ap1", "cmyc", 
        "p70s6_2", "hsp27"])
    apoptosis_tfs = set(["pro_apoptotic"])

    output_genes = proliferation_grow_tfs.union(apoptosis_tfs)


    p_value_df = pd.DataFrame()
    rank_df = pd.DataFrame()


    for output_gene in output_genes:

        print ("processing", output_gene)

        dataframe_filename = os.path.join("results", "{}_mutations.csv".format(output_gene))
        print ("reading", dataframe_filename)
        df = pd.read_csv(dataframe_filename, index_col=0)
        assert df.shape[1] == 10000
        
        original_expressions = df.loc["original"]

        modification_p_values = {}
        assert (len(modification_p_values)) == 0
        
        for modification in df.index[1:]:

            assert modification not in ["cancer", "original"]

            print ("looking at modification", modification)

            modification_expressions = df.loc[modification]

            # _, p_value = kruskal(original_expressions, 
            #     modification_expressions, 
            #     nan_policy="omit")

            assert np.nan not in original_expressions
            assert np.nan not in modification_expressions

            t_statistic, p_value = ttest_rel(original_expressions, 
                modification_expressions, 
                nan_policy="omit")

            if output_gene == "pro_apoptotic":
                # want positive t statistic
                if t_statistic < 0:
                    p_value = 1
            else:
                # expect negative t statistic
                if t_statistic > 0:
                    p_value = 1

            if np.isnan(p_value):
                p_value = 1

            # if p_value < 1e-15:
            #     p_value = 0


            assert not np.isnan(p_value)
            assert modification not in modification_p_values
            modification_p_values.update({modification: p_value})
        
        p_value_df = p_value_df.append(pd.Series(modification_p_values, name=output_gene))
        
        sorted_p_values = sorted(set(modification_p_values.values()))

        rank_dict = {modification: sorted_p_values.index(p_value) + 1 
            for modification, p_value in modification_p_values.items()}
        
        rank_df = rank_df.append(pd.Series(rank_dict, name=output_gene))

    p_value_df.to_csv("cancer_p_values.csv")
    rank_df.to_csv("cancer_rank_dataframe.csv")

    mean_over_all_outputs = rank_df.mean(axis=0).to_dict()
    sorted_modifications = sorted(mean_over_all_outputs, key= mean_over_all_outputs.get)
    for modification in sorted_modifications:
        print ("{:15s}\t{}".format(modification, mean_over_all_outputs[modification]))

if __name__ == "__main__":
    main()
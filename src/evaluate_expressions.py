import os

import glob

import numpy as np 
import pandas as pd 

from scipy.stats import kruskal, mannwhitneyu, ttest_rel

import argparse


def load_dfs(df_directory):

    print ("loading data")

    dfs = {}

    full_dfs = glob.glob(os.path.join(df_directory, 
        "*full.csv"))
    if len(full_dfs) == 0:

        print ("full dataframes not found, loading chunks")

        for f in glob.iglob(os.path.join(df_directory, "*chunk*.csv")):
            _file = f.split("/")[-1]
            output_gene = _file.split("_expressions")[0]
            df = pd.read_csv(f, index_col=0)
            if output_gene in dfs:
                dfs[output_gene] = dfs[output_gene].append(df)
            else:
                dfs[output_gene] = df
            print ("read", _file)

        for output_gene, df in dfs.items():
            df.to_csv(os.path.join(df_directory, 
                "{}_expressions_full.csv".format(output_gene)))

    else:
        print ("loading full dataframes")
        for full_df in full_dfs:
            _file = full_df.split("/")[-1]
            output_gene = _file.split("_expressions")[0]
            dfs[output_gene] = pd.read_csv(full_df, index_col=0)
    
    return dfs

def evaluate_egfr(dfs):

    p_value_df = pd.DataFrame()
    rank_df = pd.DataFrame()

    for output_gene, df in dfs.items():

        modification_p_values = {}

        cancer_expressions = df.loc["cancer"]

        mutations = df.loc[~(df.index=="cancer")]

        for mutation in mutations.index:
            mutation_expressions = df.loc[mutation]

            t_statistic, p_value = ttest_rel(cancer_expressions, 
                mutation_expressions, 
                nan_policy="omit")


            if output_gene == "pro_apoptotic":
                # want negative value
                if t_statistic > 0:
                    p_value = 1
            else:
                # expect positive
                if t_statistic < 0:
                    p_value = 1

            if np.isnan(p_value):
                p_value = 1

            # if p_value < 1e-15:
            #     p_value = 0

            assert not np.isnan(p_value)
            assert modification not in modification_p_values

            modification_p_values.update({modification: p_value})
        
        p_value_df = p_value_df.append(
            pd.Series(modification_p_values, 
            name=output_gene))
        
        sorted_p_values = sorted(set(modification_p_values.values()))

        rank_dict = {modification: sorted_p_values.index(p_value) + 1 
            for modification, p_value in modification_p_values.items()}
        
        rank_df = rank_df.append(pd.Series(rank_dict, 
            name=output_gene))

    # p_value_df.to_csv(os.path.join(expressions_dir, 
    #     "p_values.csv"))
    # rank_df.to_csv(os.path.join(expressions_dir, 
    #     "rank_dataframe.csv"))

    # mean_rank_filename = os.path.join(expressions_dir,
    #     "mean_ranks.csv")


    

def evaluate_gastric(dfs):

    anti_survival = ["Caspase8", "Caspase9", "FOXO"]
    pro_survival = ["RSK", "TCF", "cMYC"]

    for gene in anti_survival + pro_survival:
        assert gene in dfs 

    print ("computing mean activation for each mutation")
    activation_means = {output_gene: df.mean(1) 
        for output_gene, df in dfs.items()}

    print (activation_means["FOXO"])
    raise SystemExit

        
    growth_scores = 0
    for ps_output in pro_survival:
        growth_scores += activation_means[ps_output]
    for as_output in anti_survival:
        growth_scores -= activation_means[as_output]

    growth_scores = growth_scores.sort_values()

    return growth_scores
    

def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(
        description="Evaluate mutated networks")

    # parser.add_argument("--edgelist", 
    #     dest="edgelist", type=str, default=None,
    #     help="edgelist to load.")

    # parser.add_argument("--primes", 
    #     dest="primes", type=str, default=None,
    #     help="primes to load.")
    # parser.add_argument("--cancer_mutation", 
    #     dest="cancer_mutation", type=str, default=None, nargs="*",
    #     help="genes to set to overexpressed.")
    

    parser.add_argument("-d", "--df_directory", 
        dest="df_directory", 
        type=str, default=None,
        help="Directory to load results.")

    # parser.add_argument("--output_genes", dest="output_genes", 
    #     type=str, default=None, nargs="+",
    #     help="Directory to save results.")

    return parser.parse_args()


def main():

    args = parse_args()
    directory = args.df_directory

    dfs = load_dfs(directory)



    # for gastric dataset
    growth_scores = evaluate_gastric(dfs)

    growth_scores.to_csv(os.path.join(directory, 
        "growth_scores.csv"))



    # mutations = ["erbb11"]

    # expressions_dir = os.path.join("results", "_".join(mutations))
     
    # core_of_the_core = set(["pi3k", "pip3", "gab1"])

    # proliferation_grow_tfs = set(["elk1", "creb", "ap1", "cmyc", 
    #     "p70s6_2", "hsp27"])
    # apoptosis_tfs = set(["pro_apoptotic"])

    # output_genes = proliferation_grow_tfs.union(apoptosis_tfs)

    # p_value_df = pd.DataFrame()
    # rank_df = pd.DataFrame()

    # for output_gene in output_genes:

    #     print ("processing", output_gene)

    #     dataframe_filename = os.path.join(expressions_dir,
    #         "{}_expressions.csv".format(output_gene))
    #     print ("reading", dataframe_filename)
    #     df = pd.read_csv(dataframe_filename, index_col=0)
    #     assert df.shape[1] == 10000
        
    #     cancer_expressions = df.loc["cancer"]

    #     modification_p_values = {}
    #     assert (len(modification_p_values)) == 0
        
    #     for modification in df.index[2:]:

    #         assert modification not in ["cancer", "original"]

    #         print ("looking at modification", modification)

    #         modification_expressions = df.loc[modification]

    #         # _, p_value = kruskal(original_expressions, 
    #         #     modification_expressions, 
    #         #     nan_policy="omit")

    #         assert np.nan not in cancer_expressions
    #         assert np.nan not in modification_expressions

    #         t_statistic, p_value = ttest_rel(cancer_expressions, 
    #             modification_expressions, 
    #             nan_policy="omit")

    #         if output_gene == "pro_apoptotic":
    #             # want negative value
    #             if t_statistic > 0:
    #                 p_value = 1
    #         else:
    #             # expect positive
    #             if t_statistic < 0:
    #                 p_value = 1

    #         if np.isnan(p_value):
    #             p_value = 1

    #         # if p_value < 1e-15:
    #         #     p_value = 0

    #         assert not np.isnan(p_value)
    #         assert modification not in modification_p_values

    #         modification_p_values.update({modification: p_value})
        
    #     p_value_df = p_value_df.append(pd.Series(modification_p_values, name=output_gene))
        
    #     sorted_p_values = sorted(set(modification_p_values.values()))

    #     rank_dict = {modification: sorted_p_values.index(p_value) + 1 
    #         for modification, p_value in modification_p_values.items()}
        
    #     rank_df = rank_df.append(pd.Series(rank_dict, name=output_gene))

    # p_value_df.to_csv(os.path.join(expressions_dir, 
    #     "p_values.csv"))
    # rank_df.to_csv(os.path.join(expressions_dir, 
    #     "rank_dataframe.csv"))

    # mean_rank_filename = os.path.join(expressions_dir,
    #     "mean_ranks.csv")

    # mean_over_all_outputs = rank_df.mean(axis=0).to_dict()
    # sorted_modifications = sorted(mean_over_all_outputs, key= mean_over_all_outputs.get)
    # with open(mean_rank_filename, "w") as f:
    #     for modification in sorted_modifications:
    #         # if any([core_gene in modification
    #         #     for core_gene in core_of_the_core]):
    #         # if "_" not in modification:
    #             f.write("{},{}\n".format(modification,
    #                 mean_over_all_outputs[modification]))
    #             print ("{:15s}\t{}".format(modification, 
    #                 mean_over_all_outputs[modification]))




if __name__ == "__main__":
    main()
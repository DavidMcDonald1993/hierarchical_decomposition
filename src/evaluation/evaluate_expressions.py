import os

import glob

import numpy as np 
import pandas as pd 

from scipy.stats import kruskal, mannwhitneyu, ttest_rel

import argparse

def load_dfs(df_directory):

    print ("loading data from directory", df_directory)

    dfs = {}

    full_dfs = glob.glob(os.path.join(df_directory, 
        "*expressions_full.csv"))
    if len(full_dfs) == 0:

        print ("full dataframes not found, loading chunks")

        for f in glob.iglob(os.path.join(df_directory, "chunks",
            "*chunk*.csv")):
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
            print ("read", _file)

    return dfs

def evaluate_egfr(dfs, output_dir):

    print ("evaluating EGFR")

    p_value_df = pd.DataFrame()

    for output_gene, df in dfs.items():

        print ("processing output", output_gene)

        target_set_p_values = {}

        cancer_expressions = df.loc["cancer"]

        target_sets = set(df.index) - {"cancer"} #df.loc[~(df.index=="cancer")]

        print ("number of target sets:", len(target_sets))

        for target_set in target_sets:

            print ("processing target set", target_set,
                "for output", output_gene)

            target_set_expressions = df.loc[target_set]

            print ("computing t-statistic and p-value")
            t_statistic, p_value = ttest_rel(
                cancer_expressions, 
                target_set_expressions, 
                nan_policy="omit")


            if output_gene == "pro_apoptotic":
                # expect negative value
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
            assert target_set not in target_set_p_values

            target_set_p_values.update({target_set: p_value})
        
        p_value_df[output_gene] = pd.Series(target_set_p_values) 
        
      
    rank_df = p_value_df.rank(axis=0, ascending=True, method="first")
    mean_rank_df = rank_df.mean(axis=1) # mean over all outputs
    mean_rank_df = mean_rank_df.sort_values(ascending=True)

    p_value_df_filename = os.path.join(output_dir, 
        "p_values.csv")
    print ("writing p-values to", p_value_df_filename)
    p_value_df.to_csv(p_value_df_filename)

    rank_df_filename = os.path.join(output_dir, 
        "rank_dataframe.csv")
    print ("writing ranks to", rank_df_filename)
    rank_df.to_csv(rank_df_filename)

    mean_rank_filename = os.path.join(output_dir,
        "mean_ranks.csv")
    print ("writing mean ranks to", mean_rank_filename)
    mean_rank_df.to_csv(mean_rank_filename)

    



def evaluate_gastric(dfs, output_dir):

    print ("evaluating gastric")

    anti_survival = ["Caspase8", "Caspase9", "FOXO"]
    pro_survival = ["RSK", "TCF", "cMYC"]

    for gene in anti_survival + pro_survival:
        assert gene in dfs 

    print ("computing mean activation for each mutation")
    activation_means = {output_gene: df.mean(1) # mean across attractors
        for output_gene, df in dfs.items()}

    growth_scores = 0 # sum for all target_groups
    for ps_output in pro_survival:
        growth_scores += activation_means[ps_output]
    for as_output in anti_survival:
        growth_scores -= activation_means[as_output]

    growth_scores = growth_scores.sort_values()


    growth_scores.to_csv(os.path.join(output_dir, 
            "growth_scores.csv"))


def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(
        description="Evaluate mutated networks")

    parser.add_argument("-d", "--df_directory", 
        dest="df_directory", 
        type=str, default=None,
        help="Directory to load results.")

    return parser.parse_args()


def main():

    args = parse_args()
    directory = args.df_directory

    dfs = load_dfs(directory)

    if "gastric" in directory:
        # for gastric dataset
        evaluate_gastric(dfs, directory)

    else:
        evaluate_egfr(dfs, directory)

if __name__ == "__main__":
    main()
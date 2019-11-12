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
            print ("read", _file)

    # remove any mutations caused by output genes
    regex = "|".join(dfs.keys())
    for output_gene, df in dfs.items():
        dfs[output_gene] = df.loc[~df.index.str.contains(regex)]

    return dfs

def evaluate_egfr(dfs, output_dir):

    print ("evaluating EGFR")

    p_value_df = pd.DataFrame()
    rank_df = pd.DataFrame()

    for output_gene, df in dfs.items():

        mutation_p_values = {}

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
            assert mutation not in mutation_p_values

            mutation_p_values.update({mutation: p_value})
        
        p_value_df = p_value_df.append(
            pd.Series(mutation_p_values, 
            name=output_gene))
        
        sorted_p_values = sorted(set(mutation_p_values.values()))

        rank_dict = {mutation: sorted_p_values.index(p_value) + 1 
            for mutation, p_value in mutation_p_values.items()}
        
        rank_df = rank_df.append(pd.Series(rank_dict, 
            name=output_gene))

    p_value_df.to_csv(os.path.join(output_dir, 
        "p_values.csv"))
    rank_df.to_csv(os.path.join(output_dir, 
        "rank_dataframe.csv"))

    mean_rank_filename = os.path.join(output_dir,
        "mean_ranks.csv")


    mean_over_all_outputs = rank_df.mean(axis=0).to_dict()
    sorted_modifications = sorted(mean_over_all_outputs, key= mean_over_all_outputs.get)
    with open(mean_rank_filename, "w") as f:
        for modification in sorted_modifications:
                f.write("{},{}\n".format(modification,
                    mean_over_all_outputs[modification]))
                print ("{:15s}\t{}".format(modification, 
                    mean_over_all_outputs[modification]))



    

def evaluate_gastric(dfs):

    print ("evaluating gastric")

    anti_survival = ["Caspase8", "Caspase9", "FOXO"]
    pro_survival = ["RSK", "TCF", "cMYC"]

    for gene in anti_survival + pro_survival:
        assert gene in dfs 

    print ("computing mean activation for each mutation")
    activation_means = {output_gene: df.mean(1) 
        for output_gene, df in dfs.items()}

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
        growth_scores = evaluate_gastric(dfs)
        growth_scores.to_csv(os.path.join(directory, 
            "growth_scores.csv"))

    else:
        evaluate_egfr(dfs, directory)

if __name__ == "__main__":
    main()
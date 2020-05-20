import os

import glob

import numpy as np 
import pandas as pd 
import networkx as nx

from scipy.stats import ttest_ind, ttest_rel

import argparse

def load_dfs(df_directory):


    print ("loading data from directory", df_directory)

    dfs = {}

    full_dfs = glob.glob(os.path.join(df_directory, 
        "*expressions_full.csv"))
    if len(full_dfs) == 0:

        print ("full dataframes not found, loading chunks from directory",
            os.path.join(df_directory, "chunks"))

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

def evaluate_change_significance_ranks(dfs, 
    cell_growth_tfs, 
    cell_death_tfs,
    output_dir,
    one_tailed=True):

    print ("evaluating change siginificance")

    p_value_df = pd.DataFrame()

    for output_gene, df in dfs.items():

        print ("processing output", output_gene)

        target_set_p_values = {}

        cancer_expressions = df.loc["cancer"]

        target_sets = set(df.index) - {"cancer"} 

        print ("number of target sets:", len(target_sets))

        for target_set in target_sets:

            print ("processing target set", target_set,
                "for output", output_gene)

            target_set_expressions = df.loc[target_set]
            
            if cancer_expressions.equals(target_set_expressions):
                p_value = 1. # no change in expression
            else:
                print ("computing t-statistic and p-value")
                if output_gene in cell_death_tfs:
                    print ("looking for increase")
                    # expect increase in expression
                    # target set > cancer
                    t_statistic, p_value = ttest_ind(
                        target_set_expressions, 
                        cancer_expressions,
                        nan_policy="omit",
                        equal_var=False
                        )
                else:
                    assert output_gene in cell_growth_tfs, output_gene
                    print ("looking for decrease")

                    # expect decrease in expression
                    # cancer_expressions > target_set

                    t_statistic, p_value = ttest_ind(
                        cancer_expressions, 
                        target_set_expressions, 
                        nan_policy="omit",
                        equal_var=False
                        )

                assert not np.isnan(p_value), t_statistic

                if one_tailed:
                    p_value /= 2 # one tailed ttest
                    if t_statistic < 0:
                        p_value = 1 - p_value

            assert target_set not in target_set_p_values

            target_set_p_values.update(
                {target_set: p_value})
        
        p_value_df[output_gene] = pd.Series(target_set_p_values) 
    

    rank_df = p_value_df.rank(axis=0, 
        ascending=True, method="min")
    mean_rank_df = rank_df.mean(axis=1) # mean over all outputs
    mean_rank_df = mean_rank_df.sort_values(ascending=True)
    
    p_value_df_filename = os.path.join(output_dir, 
        "p_values_all_targets.csv")
    print ("writing p-values to", p_value_df_filename)
    p_value_df.to_csv(p_value_df_filename)

    rank_df_filename = os.path.join(output_dir, 
        "rank_dataframe_all_targets.csv")
    print ("writing ranks to", rank_df_filename)
    rank_df.to_csv(rank_df_filename)

    mean_rank_filename = os.path.join(output_dir,
        "mean_ranks_all_targets.csv")
    print ("writing mean ranks to", mean_rank_filename)
    mean_rank_df.to_csv(mean_rank_filename)

    # split up by number of genes
    splits = [s.split("+") for s in p_value_df.index]

    for n_genes in (1, 2, 3):

        idx = [len(s) == n_genes for s in splits]

        p_value_df_n_genes = p_value_df.loc[idx]

        rank_df_n_genes = p_value_df_n_genes.rank(
            axis=0, ascending=True, method="min")
        mean_rank_df_n_genes = rank_df_n_genes.mean(axis=1) # mean over all outputs
        mean_rank_df_n_genes = mean_rank_df_n_genes.\
            sort_values(ascending=True)

        p_value_df_filename = os.path.join(output_dir, 
            "p_values_{}_targets.csv".format(n_genes))
        print ("writing p-values to", p_value_df_filename)
        p_value_df_n_genes.to_csv(p_value_df_filename)

        rank_df_filename = os.path.join(output_dir, 
            "rank_dataframe_{}_targets.csv".format(n_genes))
        print ("writing ranks to", rank_df_filename)
        rank_df_n_genes.to_csv(rank_df_filename)

        mean_rank_filename = os.path.join(output_dir,
            "mean_ranks_{}_targets.csv".format(n_genes))
        print ("writing mean ranks to", mean_rank_filename)
        mean_rank_df_n_genes.to_csv(mean_rank_filename)


def evaluate_gastric(dfs, output_dir):

    print ("evaluating gastric")
    
    pro_survival = ["RSK", "TCF", "cMYC"]
    anti_survival = ["Caspase8", "Caspase9", "FOXO"]

    for gene in anti_survival + pro_survival:
        assert gene in dfs.keys()

    growth_scores = 0 # sum for all target sets
    for ps_output in pro_survival:
        growth_scores += dfs[ps_output]
    for as_output in anti_survival:
        growth_scores -= dfs[as_output]

    # mean across attractors
    print ("computing mean activation across attractors for each target set")
    # activation_means = {output_gene: growth_scores.mean(axis=1)
    #     for output_gene, df in dfs.items()}
    growth_scores = growth_scores.mean(axis=1) # mean across all attractors

    growth_scores = growth_scores.sort_values()
    growth_scores.to_csv(os.path.join(output_dir, 
            "growth_scores_all_targets.csv"))

    splits = [s.split("+") for s in growth_scores.index]

    for n_genes in (1, 2, 3):

        idx = [len(s) == n_genes for s in splits]
        growth_scores_n_genes = growth_scores.loc[idx]

        growth_scores_n_genes.to_csv(os.path.join(
            output_dir, 
            "growth_scores_{}_targets.csv".format(n_genes)))

def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(
        description="Evaluate targeted networks")

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
        cell_growth_tfs = {"RSK", "TCF", "cMYC"}
        cell_death_tfs = {"Caspase8", "Caspase9", "FOXO"}
        one_tailed = True
        # evaluate_gastric(dfs, directory)
    elif "egfr" in directory:
        cell_growth_tfs = {"elk1", "creb", "ap1",
            "cmyc", "p70s6_2", "hsp27"}
        cell_death_tfs = {"pro_apoptotic"}
        one_tailed = True

    else:
        assert "synthetic" in directory
        cell_growth_tfs = {"n13", "n14", "n15"}
        cell_death_tfs = set()
        one_tailed = False # do not consider direction


    for tf in cell_growth_tfs.union(cell_death_tfs):
        assert tf in dfs.keys(), tf 

    evaluate_change_significance_ranks(dfs,
        cell_growth_tfs, cell_death_tfs, directory, one_tailed=one_tailed)


if __name__ == "__main__":
    main()
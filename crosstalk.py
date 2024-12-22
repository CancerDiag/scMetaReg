import re
import random
import operator
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
from random import sample
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from collections import defaultdict

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


# Module enrichment
def perform_enrichment(gene_list1, gene_list2, pop_size):
    # Calculate the intersection of gene lists
    intersection = set(gene_list1) & set(gene_list2)
    len_genes1 = len(gene_list1)
    len_genes2 = len(gene_list2)
    len_intersection = len(intersection)

    # Perform hypergeometric test
    p_value = 1- hypergeom.cdf(len_intersection - 1, pop_size, len_genes1, len_genes2)

    return p_value

def get_pvalue_dct(given_genelst):
    total_genes =  # Number of features after scTransform normalisation
    pvalue_dct = {}
    for tf, target_lst in tf_target_dct.items():
        pvalue = perform_enrichment(given_genelst, target_lst, total_genes)
        pvalue_dct[tf] = pvalue
    return pvalue_dct

def get_sig_tfs(pvalue_dct): 
    sig_tfs_df = pd.DataFrame.from_dict(pvalue_dct, orient='index')
    sig_tfs_df.columns = ["p_value"]
    
    pvalue_lst = list(sig_tfs_df["p_value"].values)
    adjusted_p_values = multipletests(pvalue_lst, method='fdr_bh')[1]
    sig_tfs_df["adj_pvalue"] = list(adjusted_p_values)
    
    sig_tfs_df = sig_tfs_df.sort_values('adj_pvalue', ascending=True)    
    sig_tfs_df = sig_tfs_df[sig_tfs_df["adj_pvalue"] < 0.01]
    print(sig_tfs_df.shape)
    print(sig_tfs_df.index)
    return sig_tfs_df


def get_module_results(nodes):
    module_nodes = list(nodes)
    print("module size:",len(module_nodes))
    module_pvalue_dct = get_pvalue_dct(module_nodes)
    module_tfs_df = get_sig_tfs(module_pvalue_dct)
    return module_tfs_df




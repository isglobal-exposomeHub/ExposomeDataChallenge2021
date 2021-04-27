# Author Jean-Baptiste Guimbaud
# Meersens

import matplotlib.pyplot as plt
import numpy as np
import scipy
import pandas as pd
import scipy.cluster.hierarchy as sch

def plot_corr(corr,size=10):
    '''Plot a graphical correlation matrix for a dataframe'''
    fig, ax = plt.subplots(figsize=(size, size))
    fig.subplots_adjust(bottom=0.05, top=0.7)
    cax = ax.matshow(corr, cmap='RdYlGn')
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90, fontsize=6)
    plt.yticks(range(len(corr.columns)), corr.columns)
    
    # Add the colorbar legend
    cbar = fig.colorbar(cax, ticks=[-1, 0, 1], aspect=40, shrink=.8)
    plt.show()


def cluster_corr(df, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    """
    cluster_th = 4

    X = df.corr().values
    d = sch.distance.pdist(X)
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')

    columns = [df.columns.tolist()[i] for i in list(np.argsort(ind))]
    df = df.reindex(columns, axis=1)

    unique, counts = np.unique(ind, return_counts=True)
    counts = dict(zip(unique, counts))

    i = 0
    j = 0
    columns = []
    for cluster_l1 in set(sorted(ind)):
        j += counts[cluster_l1]
        sub = df[df.columns.values[i:j]]
        if counts[cluster_l1]>cluster_th:        
            X = sub.corr().values
            d = sch.distance.pdist(X)
            L = sch.linkage(d, method='complete')
            ind = sch.fcluster(L, 0.5*d.max(), 'distance')
            col = [sub.columns.tolist()[i] for i in list((np.argsort(ind)))]
            sub = sub.reindex(col, axis=1)
        cols = sub.columns.tolist()
        columns.extend(cols)
        i = j
    return df.reindex(columns, axis=1)
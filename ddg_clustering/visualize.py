from .clustering import AA_COLS, LOS, GOS
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np
from .clustering import ClusterType, Metric
import matplotlib.colors as mcolors
# DDG HEATMAP
def ddg_heatmap(gene, gene_data, chain='A', start=None, end=None):
    df = gene_data

    chain_data = df[df['chain'] == chain]
    if start != None:
        chain_data = chain_data[chain_data["residue_id"] >= start]
    if end != None:
        chain_data = chain_data[chain_data["residue_id"] <= end]   
    
    data = chain_data[AA_COLS]
    data_mean = data.mean(axis=1)
    data_mean.name = "avg"

    height = 4
    # width is proportional to the number of residues
    width = len(data) / 100
    figsize = (width, height)
    plt.figure(figsize=figsize)

    # Create the heatmap
    combined_data = pd.concat([data, data_mean], axis=1).T
    ax = sns.heatmap(combined_data, 
                    center=0,
                    vmin=GOS-0.5, vmax=LOS+0.5,
                    cmap='jet')

    ax.set_yticks(range(len(combined_data.index)))
    labels = combined_data.index
    new_labels = [label if i % 2 == 0 else f"{label}   " for i, label in enumerate(labels)]
    ax.set_yticklabels(new_labels, rotation=0)

    # make x-axis labels appear on top and rotate 180
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('Alternative Residue')
    ax.set_title(f'{gene} ∆∆G values from in silico mutagenesis')

    # Optional: Adjust font size for better readability
    plt.tick_params(axis='y', which='major', labelsize=8)

    plt.show()

## VISUALIZE THRESHOLD HITS
def visualize_threshold_hits(gene, clust_df, cluster_type=ClusterType.LOSS_OF_STABILITY, metric=Metric.SIGNIFICANCE,
                              pvalue=0.05, distribution=stats.gengamma):

    # combine all to a single array and 1 column
    all_values = np.array(clust_df[AA_COLS]).flatten()

    if cluster_type == ClusterType.LOSS_OF_STABILITY:
        all_values = all_values[all_values >= 0]
    else:
        all_values = all_values[all_values < 0]

    if metric == Metric.SIGNIFICANCE:
        # fit a gengamma distribution to the positive values
        fit_values = all_values if cluster_type == ClusterType.LOSS_OF_STABILITY else -all_values
        params = distribution.fit(fit_values)

        quantile = distribution.ppf(1-pvalue, *params)
        cutoff = 1-pvalue
        print(f"{round(cutoff*100)}% quantile: {quantile}")


    # replot the positive distribution with a vertical red line at the 95% quantile
    sns.histplot(all_values)
    if metric == Metric.SIGNIFICANCE:
        quantile = quantile if cluster_type == ClusterType.LOSS_OF_STABILITY else -quantile
        plt.axvline(quantile, color='r', linestyle='--', label=f'{round(cutoff*100)}% quantile: {quantile:.2f}')
    else:
        cutoff = LOS if cluster_type == ClusterType.LOSS_OF_STABILITY else GOS
        plt.axvline(cutoff, color='r', linestyle='--', label=f'cutoff: {cutoff:.2f}')
    plt.legend()
    if metric == Metric.SIGNIFICANCE:
        plt.title(f"{gene} {str(cluster_type)} ∆∆G distribution with significance cutoff at p={pvalue}")
    else:
        plt.title(f"{gene} {str(cluster_type)} ∆∆G distribution with cutoff at {cutoff}")
    plt.xlabel("∆∆G")
    if cluster_type == ClusterType.LOSS_OF_STABILITY:
        plt.xlim(-0.5, 10)
    else:
        plt.xlim(-5, 0.5)
    plt.show()


## VISUALIZE CLUSTERS ON SEQUENCE
def visualize_clusters_on_sequence(gene, gene_data, cluster_type=ClusterType.LOSS_OF_STABILITY):
        top_clusters = gene_data[gene_data['cluster'].notnull()]

        number_of_clusters = len(top_clusters['cluster'].unique())
        colors = sns.color_palette('hls', number_of_clusters - 1)
        # set nans to be number_of_clusters
        gene_data['cluster'] = gene_data['cluster'].fillna(number_of_clusters)
        custom_palette = colors + [mcolors.to_rgba('grey')] * (1)
        discrete_cmap = mcolors.ListedColormap(custom_palette)

        height = 4
        # width is proportional to the number of residues
        width = len(gene_data) / 100
        figsize = (width, height)
        plt.figure(figsize=figsize)

        sns.heatmap(gene_data['cluster'].values.reshape(1, -1), cmap=discrete_cmap, yticklabels=False)

        plt.title(f'{str(cluster_type)} clusters for {gene}')
        plt.ylabel(str(cluster_type)) 

        plt.show()



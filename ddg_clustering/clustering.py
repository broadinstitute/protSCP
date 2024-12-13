from enum import Enum
import numpy as np
from scipy import stats
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import pandas as pd

# suppress the pandas warning for assigning a slice
pd.options.mode.chained_assignment = None

class Metric(Enum):
    SIGNIFICANCE = 1
    THRESHOLD = 2

class ClusterType(Enum):
    LOSS_OF_STABILITY = 1
    GAIN_OF_STABILITY = 2

    # Adding a tostring function
    def __str__(self):
        return self.name.replace("_", " ").title()

AA_COLS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

LOS = 1.38
GOS = -0.98


def plot_clusters(gene, gene_data, metric, los_threshold, figsize=(20,5), plot_counts=False):
        top_clusters = gene_data[gene_data['cluster'].notnull()]
        if plot_counts:
            # print residue type value counts
            top_5 = top_clusters['aa'].value_counts().nlargest(3)
            # Replace other values with 'other'
            top_clusters['aa'] = top_clusters['aa'].apply(lambda x: x if x in top_5.index else 'other')

            # quick barplot of residue type counts
            top_clusters['aa'].value_counts().plot(kind='bar')
            # title it with gene, threshold, metric, and residue type counts
            plt.title(f'{gene} threshold = {los_threshold} {str(metric)} residue type counts')
            # annotate the bars with the number
            for i, v in enumerate(top_clusters['aa'].value_counts()):
                plt.text(i - 0.1, v + 0.5, str(v))

            plt.show()

            # quick barplot of ss3 counts
            top_clusters['ss3'].value_counts().plot(kind='bar')
            # title it with gene, threshold, metric, and ss3 counts
            plt.title(f'{gene} threshold = {los_threshold} {str(metric)} ss3 counts')
            # annotate the bars with the number
            for i, v in enumerate(top_clusters['ss3'].value_counts()):
                plt.text(i - 0.1, v + 0.5, str(v))
            plt.show()

        number_of_clusters = len(top_clusters['cluster'].unique())
        colors = sns.color_palette('hls', number_of_clusters - 1)
        # set nans to be number_of_clusters
        gene_data['cluster'] = gene_data['cluster'].fillna(number_of_clusters)
        custom_palette = colors + [mcolors.to_rgba('grey')] * (1)
        discrete_cmap = mcolors.ListedColormap(custom_palette)

        plt.figure(figsize=figsize)
        sns.heatmap(gene_data['cluster'].values.reshape(1, -1), cmap=discrete_cmap, yticklabels=False)

        plt.title(f'{str(metric)} for {gene} threshold = {los_threshold}')
        plt.ylabel(str(metric)) 

        plt.show()


## EXTERNAL: PRIORITIZING SITES
def prioritize_sites(gene_data, pvalue=0.05, distribution=None, metric=None, type=ClusterType.LOSS_OF_STABILITY, cutoff=15):
    LoS = (type == ClusterType.LOSS_OF_STABILITY)

    gene_data = add_los_gos(gene_data)


    if distribution is None:
        distribution = stats.gengamma
    if metric is None:
        metric = Metric.SIGNIFICANCE

    ddg_cols = gene_data[AA_COLS]
    all_values = np.array(ddg_cols)

    # Step 1: Separate the positive and negative numbers
    positive_values = all_values[all_values > 0]
    negative_values = all_values[all_values < 0]

    # flatten the values to 1d
    params_positive = distribution.fit(positive_values, floc=0)  # Fit with location fixed at 0
    params_negative = distribution.fit(-negative_values, floc=0)  # Fit with location fixed at 0

    if metric == Metric.SIGNIFICANCE:
        if LoS:
            significance_level = 1 - pvalue
            cutoff = distribution.ppf(significance_level, *params_positive)
            over_cutoff = gene_data[(gene_data[AA_COLS] > cutoff).any(axis=1)]
        else:
            significance_level = 1 - pvalue
            cutoff = -1 * distribution.ppf(significance_level, *params_negative)
            over_cutoff = gene_data[(gene_data[AA_COLS] < cutoff).any(axis=1)]
    else:
        if LoS:
            over_cutoff = gene_data[gene_data["los"] > cutoff]
        else:
            over_cutoff = gene_data[gene_data["gos"] < cutoff]
    
    return over_cutoff


## EXTERNAL: CLUSTERING PRIORITIZED SITES
def cluster_prioritized_sites(prioritized_sites, gene_data, distance_threshold=6, min_cluster_size=5):
    xyz = prioritized_sites[['xca', 'yca', 'zca']].values

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=distance_threshold, linkage='single')
    clustering.fit(xyz)

    prioritized_sites['cluster'] = clustering.labels_

    new_gene_data = gene_data.copy()
    new_gene_data.loc[prioritized_sites.index, 'cluster'] = clustering.labels_

    cluster_labels = clustering.labels_
    unique_clusters, counts = np.unique(cluster_labels, return_counts=True)
    sorted_indices = np.argsort(-counts)

    rank_mapping = {unique_clusters[sorted_indices[i]]: i + 1 for i in range(len(unique_clusters))}
    new_gene_data['cluster'] = new_gene_data['cluster'].apply(lambda x: rank_mapping.get(x, len(unique_clusters)))

    # get clusters with at least 5 residues
    # step 1: get the value counts and filter for counts >= min_cluster_size. find the biggest cluster number from those
    valid_clusters = new_gene_data['cluster'].value_counts()[new_gene_data['cluster'].value_counts() >= min_cluster_size].index
    # filter out cluster with index = unique_clusters.size
    valid_clusters = valid_clusters[valid_clusters < unique_clusters.size]

    smallest_valid_cluster = valid_clusters.max()

    # in gene_data, set all rows with a cluster > smallest_valid_cluster to nan
    new_gene_data.loc[new_gene_data['cluster'] > smallest_valid_cluster, 'cluster'] = np.nan

    cluster_sizes = {rank_mapping[cluster]: count for cluster, count in zip(unique_clusters, counts)}
    
    return new_gene_data, cluster_sizes


## EXTERNAL: SAVING CLUSTERS
def save_clusters(output_filename, cluster_output, keep_aa=False):
    case_clusters = cluster_output[0]

    if not keep_aa:
        # remove AA cols and xca, yca, zca, and plddt
        case_clusters = case_clusters.drop(AA_COLS + ['xca', 'yca', 'zca', 'plddt'], axis=1)

    case_clusters.to_csv(output_filename, index=False)

# Formatting DF for hit prioritization
def add_los_gos(df):
    new_df = df.copy()

    # Counting the number of columns where the value is greater than LOS
    new_df['los'] = new_df[AA_COLS].apply(lambda row: sum(value > LOS for value in row), axis=1)
    # Counting the number of columns where the value is less than GOS
    new_df['gos'] = new_df[AA_COLS].apply(lambda row: sum(value < GOS for value in row), axis=1)

    return new_df

#################### INTERNAL/UNUSED FOR PIPELINE: CLUSTERING WITH CONTROL ####################

def find_gene_data(filename=None, residue_wise_df=None):     
    if filename is not None:
        # Load the gene data
        df = pd.read_csv(filename)
    elif residue_wise_df is not None:
        df = residue_wise_df

    # If the df has all the cols in AA_COLS, then get these and add an los and gos column for the number of columns where
    # the value is > LOS or < GOS
    if all(col in df.columns for col in AA_COLS):
        # Counting the number of columns where the value is greater than LOS
        df['los'] = df[AA_COLS].apply(lambda row: sum(value > LOS for value in row), axis=1)
        # Counting the number of columns where the value is less than GOS
        df['gos'] = df[AA_COLS].apply(lambda row: sum(value < GOS for value in row), axis=1)

    num_chains = df['chain'].nunique()
    if num_chains > 1:
        df.index = df.apply(lambda row: f"{row['chain']}:{row['aa']}{row['residue_id']}", axis=1)
    else:
        df.index = df['aa'] + df['residue_id'].astype(str) # Fix the index to be residue_index + type

    return df


def perform_clustering_wrapped(gene_data, los_threshold=10, metric=Metric.SIGNIFICANCE, LoS=True):
    ddg_cols = gene_data[AA_COLS]
    all_values = np.array(ddg_cols)

    # Step 1: Separate the positive and negative numbers
    positive_values = all_values[all_values > 0]
    negative_values = all_values[all_values < 0]

    # flatten the values to 1d
    params_positive = stats.gengamma.fit(positive_values, floc=0)  # Fit with location fixed at 0
    params_negative = stats.gengamma.fit(-negative_values, floc=0)  # Fit with location fixed at 0

    if metric == Metric.SUM or metric == Metric.SUM_CAPPED:
        pass
    elif metric == Metric.SIGNIFICANCE:
        if LoS:
            significance_level = 0.95
            cutoff = stats.gengamma.ppf(significance_level, *params_positive)
            over_cutoff = gene_data[(gene_data[AA_COLS] > cutoff).any(axis=1)]
        else:
            significance_level = 0.95
            cutoff = -1 * stats.gengamma.ppf(significance_level, *params_negative)
            over_cutoff = gene_data[(gene_data[AA_COLS] < cutoff).any(axis=1)]
        
    elif metric == Metric.MEAN_STD:
        # get the mean and std dev for each residue
        pass
    else:
        over_cutoff = gene_data[gene_data["los"] > los_threshold]
    xyz = over_cutoff[['xca', 'yca', 'zca']].values

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=6, linkage='single')
    clustering.fit(xyz)

    gene_data.loc[over_cutoff.index, 'cluster'] = clustering.labels_

    cluster_labels = clustering.labels_
    unique_clusters, counts = np.unique(cluster_labels, return_counts=True)
    sorted_indices = np.argsort(-counts)

    rank_mapping = {unique_clusters[sorted_indices[i]]: i + 1 for i in range(len(unique_clusters))}
    gene_data['cluster'] = gene_data['cluster'].apply(lambda x: rank_mapping.get(x, len(unique_clusters)))

    # get clusters with at least 5 residues
    # step 1: get the value counts and filter for counts >= 5. find the biggest cluster number from those
    valid_clusters = gene_data['cluster'].value_counts()[gene_data['cluster'].value_counts() >= 5].index
    # filter out cluster with index = unique_clusters.size
    valid_clusters = valid_clusters[valid_clusters < unique_clusters.size]

    smallest_valid_cluster = valid_clusters.max()

    # in gene_data, set all rows with a cluster > smallest_valid_cluster to nan
    gene_data.loc[gene_data['cluster'] > smallest_valid_cluster, 'cluster'] = np.nan

    cluster_sizes = {rank_mapping[cluster]: count for cluster, count in zip(unique_clusters, counts)}
    
    return gene_data, cluster_sizes


def perform_clustering(gene="BRCA1", filename=None, residue_wise_df=None, type=ClusterType.LOSS_OF_STABILITY, cutoff=None):
    LoS = (type == ClusterType.LOSS_OF_STABILITY)

    if cutoff is None:
        cutoff = 5 # default cutoff
    
    return cluster_case_and_control(gene=gene, filename=filename, 
                                    residue_wise_df=residue_wise_df, plot_seq=True, LoS=LoS,
                                    metric=Metric.SIGNIFICANCE, threshold=cutoff, return_clusters=True)


def cluster_case_and_control(gene="TAOK1", metric=Metric.SIGNIFICANCE, LoS=True, threshold=10, repeats=2, plot=False, return_clusters=False, rotate_background=False, plot_seq=False,
                       filename=None, residue_wise_df=None):
    # 1. load the gene data
    # 2. split positive and negative values
    # 3. background: get 20xN values samples from the positive distribution (with 0s for negatives in the distribution)
    # 4. Get clusters for 20xN background
    # 5. Get clusters for positive values
    gene_data = find_gene_data(filename=filename, residue_wise_df=residue_wise_df)
    ddg_cols = gene_data[AA_COLS]

    # get the distribution of all values in ddg cols
    all_values = ddg_cols.values.flatten()

    if plot:
        # plot the distribution of all values
        sns.histplot(all_values, bins=100)
        plt.show()


    # construct the background
    background_clusters = []
    background_cluster_sizes = []
    if rotate_background:
        # rotate all AA cols by 100 each way
        # Rotate rows downward by 100
        rotate_down = ddg_cols.shift(100, axis=0)  # Shift down by 100 rows
        rotate_down.iloc[:100, :] = ddg_cols.iloc[-100:, :].values  # Fill the top 100 rows with the last 100 rows from original

        # Rotate rows upward by 100
        rotate_up = ddg_cols.shift(-100, axis=0)  # Shift up by 100 rows
        rotate_up.iloc[-100:, :] = ddg_cols.iloc[:100, :].values  # Fill the bottom 100 rows with the first 100 rows from original


        background_df_left = pd.DataFrame(rotate_up, columns=AA_COLS)
        background_df_right = pd.DataFrame(rotate_down, columns=AA_COLS)

        background_df_left.index = gene_data.index
        background_df_right.index = gene_data.index

        background_df_left = pd.concat([background_df_left, gene_data.drop(AA_COLS, axis=1)], axis=1)
        background_df_right = pd.concat([background_df_right, gene_data.drop(AA_COLS, axis=1)], axis=1)

        # get the clusters for the background
        clusters_left, cluster_sizes_left = perform_clustering_wrapped(background_df_left, metric=metric, los_threshold=threshold, LoS=LoS)
        clusters_right, cluster_sizes_right = perform_clustering_wrapped(background_df_right, metric=metric, los_threshold=threshold, LoS=LoS)
        background_clusters.append(clusters_left)
        background_clusters.append(clusters_right)
        background_cluster_sizes.append(cluster_sizes_left)
        background_cluster_sizes.append(cluster_sizes_right)

        if plot:
            plot_clusters(gene, clusters_left, metric + f" control (left)", threshold, plot_counts=True)
            plot_clusters(gene, clusters_right, metric + f" control (right)", threshold, plot_counts=True)
    else:
        # sample 20xN values from the distribution
        # use a hardcoded random seed for reproducibility
        np.random.seed(42)
        for i in range(repeats):

            # Example data: Replace with your actual data
            all_values = np.array([ddg_cols])  # Replace with your actual data array

            # Step 1: Separate the positive and negative numbers
            positive_values = all_values[all_values > 0]
            negative_values = all_values[all_values < 0]

            # Step 2: Fit Gamma distributions to each subset
            pos_params = stats.gengamma.fit(positive_values, floc=0)  # Fit with location fixed at 0
            neg_params = stats.gengamma.fit(-negative_values, floc=0)  # Fit with location fixed at 0

            # Step 3: Calculate proportions of positive and negative numbers
            num_positive = len(positive_values)
            num_negative = len(negative_values)
            total = num_positive + num_negative

            prob_positive = num_positive / total
            prob_negative = num_negative / total

            # Step 4: Generate random samples
            sample_size = ddg_cols.shape  # Desired shape of the output

            # Generate an array indicating which distribution to sample from
            sample_distribution = np.random.choice(
                ['positive', 'negative'], size=sample_size, p=[prob_positive, prob_negative])

            # Initialize an array for the sampled values
            background = np.zeros(sample_size)

            # Sample from the appropriate distribution based on the sample_distribution array
            for index, dist in np.ndenumerate(sample_distribution):
                if dist == 'positive':
                    background[index] = stats.gengamma.rvs(*pos_params)
                else:
                    background[index] = -1 * stats.gengamma.rvs(*neg_params)

            background_df = pd.DataFrame(background, columns=AA_COLS)
            # set background_df index to gene_data index for concating
            background_df.index = gene_data.index

            # add back all other columns
            background_df = pd.concat([background_df, gene_data.drop(AA_COLS, axis=1)], axis=1)

            # get the clusters for the background
            clusters, cluster_sizes = perform_clustering_wrapped(background_df, metric=metric, los_threshold=threshold, LoS=LoS)
            background_clusters.append(clusters)
            background_cluster_sizes.append(cluster_sizes)

            if plot or plot_seq:
                plot_clusters(gene, clusters, str(metric) + f" control ({i})", threshold, plot_counts=plot)

    # get the clusters for the positive values
    positive_clusters, positive_cluster_sizes = perform_clustering_wrapped(gene_data, metric=metric, los_threshold=threshold, LoS=LoS)
    if plot or plot_seq:
        plot_clusters(gene, positive_clusters, metric, threshold, plot_counts=plot)

    #print(background_cluster_sizes, positive_cluster_sizes)
    if return_clusters:
        return (background_clusters, positive_clusters)
    return (background_cluster_sizes, positive_cluster_sizes)

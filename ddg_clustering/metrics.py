from __future__ import annotations

from typing import List, Optional, Dict

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score as _silhouette_score
from scipy import stats

from .clustering import (
    AA_COLS,
    Metric,
    ClusterType,
    prioritize_sites,
    cluster_prioritized_sites,
)


######################### SILHOUETTE #########################
def compute_silhouette_score(cluster_df: pd.DataFrame, min_cluster_size: Optional[int] = None) -> Optional[float]:
    """Compute silhouette score on 3D coordinates for clustered residues.

    Requires columns: "xca", "yca", "zca", and "cluster".
    Returns None if fewer than 2 clusters are present or if computation fails.
    """
    labeled = cluster_df.dropna(subset=["cluster"]).copy()
    if min_cluster_size is not None and min_cluster_size > 1:
        counts = labeled["cluster"].value_counts()
        large_clusters = counts[counts >= min_cluster_size].index
        labeled = labeled[labeled["cluster"].isin(large_clusters)]
    num_clusters = labeled["cluster"].nunique()
    if num_clusters < 2:
        return None
    try:
        return float(_silhouette_score(labeled[["xca", "yca", "zca"]], labeled["cluster"]))
    except Exception:
        return None


######################### GAP STATISTIC #########################
def _within_cluster_sum_of_squares(clusters_df: pd.DataFrame) -> float:
    target_cols = ["xca", "yca", "zca", "cluster"]
    clustered_points = clusters_df.dropna(subset=["cluster"])[target_cols]
    per_cluster_sums: List[float] = []
    for cluster_id in clustered_points["cluster"].unique():
        cluster_points = clustered_points[clustered_points["cluster"] == cluster_id]
        cx, cy, cz = cluster_points["xca"].mean(), cluster_points["yca"].mean(), cluster_points["zca"].mean()
        per_cluster_sums.append(
            float(np.sum((cluster_points["xca"] - cx) ** 2 +
                         (cluster_points["yca"] - cy) ** 2 +
                         (cluster_points["zca"] - cz) ** 2))
        )
    if not per_cluster_sums:
        return float("nan")
    return float(np.mean(per_cluster_sums))


def compute_gap_statistic(cluster_df: pd.DataFrame, randomized_cluster_dfs: List[pd.DataFrame], adjusted: bool = True) -> Optional[float]:
    """Compute Gap statistic adjusted by number of clusters (k^(2/3) scaling).
    Returns None on failure.
    """
    try:
        labeled = cluster_df.dropna(subset=["cluster"])  # exclude unclustered
        k = labeled["cluster"].nunique()
        if k == 0:
            return None
        scale = k ** (2.0 / 3.0)
        if not adjusted:
            scale = 1
        within = scale * _within_cluster_sum_of_squares(labeled)

        randomized_scaled: List[float] = []
        for df in randomized_cluster_dfs:
            r_labeled = df.dropna(subset=["cluster"])  # exclude unclustered
            rk = r_labeled["cluster"].nunique()
            rscale = rk ** (2.0 / 3.0) if rk > 0 else 0.0
            if not adjusted:
                rscale = 1
            r_within = _within_cluster_sum_of_squares(r_labeled)
            if np.isfinite(r_within) and r_within > 0 and rscale > 0:
                randomized_scaled.append(rscale * r_within)

        if not randomized_scaled or not np.isfinite(within) or within <= 0:
            return None
        return float(np.log(np.mean(randomized_scaled)) - np.log(within))
    except Exception:
        return None


######################### RANDOMIZED CLUSTERS #########################
def generate_randomized_clusters(
    gene_data: pd.DataFrame,
    num_random: int = 10,
    seed: Optional[int] = 42,
    LoS: bool = True,
    metric: Metric = Metric.SIGNIFICANCE,
    threshold: int = 10,
    distance_threshold: int = 6,
    min_cluster_size: int = 5,
    pvalue: float = 0.05,
) -> List[pd.DataFrame]:
    """Generate randomized cluster DataFrames for Gap computation.

    - Fits generalized gamma to positive and negative ∆∆G distributions over AA columns
    - Samples a background matrix with same shape as AA columns
    - Restores other columns (including coordinates)
    - Clusters each background and returns a list of cluster DataFrames

    Requires gene_data to include AA columns and 3D coords ("xca","yca","zca").
    """
    rng = np.random.default_rng(seed)

    ddg_cols_df = gene_data[AA_COLS]
    ddg_values = ddg_cols_df.to_numpy()

    positive_values = ddg_values[ddg_values > 0]
    negative_values = ddg_values[ddg_values < 0]

    # If no positive or negative values, bail early
    if positive_values.size == 0 or negative_values.size == 0:
        return []

    # Fit distributions with loc fixed at 0 as in pipeline
    pos_params = stats.gengamma.fit(positive_values, floc=0)
    neg_params = stats.gengamma.fit(-negative_values, floc=0)

    num_positive = positive_values.size
    num_negative = negative_values.size
    total = num_positive + num_negative
    prob_positive = num_positive / total
    prob_negative = num_negative / total

    randomized_clusters: List[pd.DataFrame] = []
    other_cols_df = gene_data.drop(columns=AA_COLS)

    for _ in range(num_random):
        # Decide per-entry whether to sample from pos or neg
        sample_distribution = rng.choice(
            ["positive", "negative"], size=ddg_cols_df.shape, p=[prob_positive, prob_negative]
        )

        background = np.zeros(ddg_cols_df.shape, dtype=float)
        # Vectorized sampling by mask
        pos_mask = sample_distribution == "positive"
        neg_mask = ~pos_mask
        pos_count = int(pos_mask.sum())
        neg_count = int(neg_mask.sum())
        if pos_count > 0:
            background[pos_mask] = stats.gengamma.rvs(*pos_params, size=pos_count, random_state=rng)
        if neg_count > 0:
            background[neg_mask] = -1.0 * stats.gengamma.rvs(*neg_params, size=neg_count, random_state=rng)

        background_df = pd.DataFrame(background, columns=AA_COLS, index=gene_data.index)
        background_df = pd.concat([background_df, other_cols_df], axis=1)

        # Prioritize then cluster using the public API
        try:
            cluster_type = ClusterType.LOSS_OF_STABILITY if LoS else ClusterType.GAIN_OF_STABILITY
            prioritized = prioritize_sites(
                background_df,
                pvalue=pvalue,
                distribution=None,
                metric=metric,
                type=cluster_type,
                cutoff=threshold,
            )

            # If no prioritized residues, return a dataframe with NaN clusters
            if prioritized is None or prioritized.empty:
                empty_df = background_df.copy()
                empty_df["cluster"] = np.nan
                randomized_clusters.append(empty_df)
                continue

            clustered_df, _ = cluster_prioritized_sites(
                prioritized,
                background_df,
                distance_threshold=distance_threshold,
                min_cluster_size=min_cluster_size,
            )
            randomized_clusters.append(clustered_df)
        except Exception:
            # On failure, fall back to returning background with NaN clusters to keep list length
            empty_df = background_df.copy()
            empty_df["cluster"] = np.nan
            randomized_clusters.append(empty_df)

    return randomized_clusters


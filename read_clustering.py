from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import Birch

from enum import Enum, auto
import numpy as np


class ClusterType(Enum):
    KMEANS = auto()
    BIRCH = auto()
    AGGLOMERATIVE = auto()
    DENSITY = auto()


def clustering(mod_table, types=None, clusters=2, eps_iter=None, min_samples=None,
               min_divisor=5, eps_divisor=4, eps_step=1, eps_range=10, neg_weight=3):
    print("Clustering")

    if not types:
        types = [ClusterType.KMEANS, ClusterType.BIRCH, ClusterType.AGGLOMERATIVE, ClusterType.DENSITY]
    if not min_samples:
        min_samples = mod_table.shape[0] // min_divisor
    if not eps_iter:
        eps_iter = _get_eps_iter(mod_table, eps_divisor, eps_step, eps_range)

    clustered = _get_clustered(_get_transformed(mod_table), types, clusters, eps_iter, min_samples, neg_weight)

    return _split_ids(clustered, mod_table.index)


def _get_eps_iter(mod_table, divisor, step, eps_range):
    middle = mod_table.shape[1] // divisor
    low = middle - eps_range
    if low < 1:
        low = 1
    return range(low, middle + eps_range, step)


def _split_ids(clustered, mod_index):
    ids_clustered = {}

    for clust_type in clustered:
        ids_clustered[clust_type] = {}
        cluster_list = clustered[clust_type]
        for cluster, r_name in zip(cluster_list, mod_index):
            if cluster not in ids_clustered[clust_type].keys():
                ids_clustered[clust_type][cluster] = []
            ids_clustered[clust_type][cluster].append(r_name)

    return ids_clustered


def _get_transformed(mod_table):
    transform_pipeline = make_pipeline(
        StandardScaler(),
        KNNImputer(n_neighbors=len(mod_table) // 4),
    )

    return transform_pipeline.fit_transform(mod_table)


def _get_clustered(transformed_table, types, clusters, eps_iter, min_samples, neg_weight):
    result = {}

    if ClusterType.KMEANS in types:
        result[ClusterType.KMEANS] = KMeans(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.BIRCH in types:
        result[ClusterType.BIRCH] = Birch(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.AGGLOMERATIVE in types:
        result[ClusterType.AGGLOMERATIVE] = AgglomerativeClustering(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.DENSITY in types:
        result[ClusterType.DENSITY] = __density_cluster(transformed_table, eps_iter, min_samples, neg_weight)

    return result


def __density_cluster(transformed_table, eps_iter, min_samples, neg_weight):
    results = []

    for eps in eps_iter:
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        dbscan.fit(transformed_table)
        results.append(dbscan.labels_)

    return __choose_best_density(results, neg_weight)


def __choose_best_density(results, neg_weight):
    best = None
    score = -np.inf

    for possibility in results:
        temp_score = __count_density_score(possibility, neg_weight)
        if temp_score >= score:
            score = temp_score
            best = possibility

    return best


def __count_density_score(possibility, neg_weight):
    unique, counts = np.unique(possibility, return_counts=True)
    dict_counts = dict(zip(unique, counts))
    negative = dict_counts.get(-1, 0) * neg_weight

    return sum([-abs(dict_counts[x] - len(possibility) / max(len(dict_counts), 2))
                for x in dict_counts if x >= 0]) - negative

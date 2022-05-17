import numpy as np

from enum import Enum, auto

from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.impute import KNNImputer
from sklearn.metrics import silhouette_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class ClusterType(Enum):
    """
    Enum of available types of clustering
    """

    KMEANS = auto()
    MINIK = auto()
    AGGLOMERATIVE = auto()
    DENSITY = auto()


def clustering(mod_table, types, imputer_divisor, clusters, min_clusters, max_clusters, hac_linkage,
               min_samples, min_divisor, eps_divisor, eps_step, eps_range, neg_weight):
    """
    Function for clustering epigenetic data.

    Parameters:
    :param mod_table: pandas table of methylation on specific positions and reads
    :param types: string with characters A,B,D,K defining which types of clustering will be used
    :param imputer_divisor: divisor used for choosing amount of neighbours when imputing
    :param clusters: Number of clusters to be used in clustering (does not apply to DBSCAN)
    :param min_clusters: Minimum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param max_clusters: Maximum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param hac_linkage: Linkage to be used in HAC clustering
    :param min_samples: minimum samples in DBSCAN cluster
    :param min_divisor: divisor for automatic calculation of minimum samples (only applies to DBSCAN)
    :param eps_divisor: divisor for automatic calculation of epsilon value (only applies to DBSCAN)
    :param eps_step: distance between tested epsilon values (only applies to DBSCAN)
    :param eps_range: broadness of tested epsilon values (only applies to DBSCAN)
    :param neg_weight: weight of outliers when choosing best DBSCAN clustering


    :return: dictionary of types of clustering and dictionary of clusters and split ids
    """

    print("Clustering")

    types = types.upper()
    if not min_samples:
        min_samples = mod_table.shape[0] // min_divisor
    eps_iter = _get_eps_iter(mod_table, eps_divisor, eps_step, eps_range)

    clustered = _get_clustered(_get_transformed(mod_table, imputer_divisor),
                               types, clusters, min_clusters, max_clusters, hac_linkage,
                               eps_iter, min_samples, neg_weight)

    return _split_ids(clustered, mod_table.index)


def _get_eps_iter(mod_table, divisor, step, eps_range):
    """
    Function for creating iterator over chosen eps values

    Parameters:
    :param mod_table:  pandas table of methylation on specific positions and reads
    :param divisor: divisor for automatic calculation of epsilon value
    :param step: distance between tested epsilon values
    :param eps_range: broadness of tested epsilon values

    :return: iterator based on parameters
    """

    middle = mod_table.shape[1] // divisor
    low = middle - eps_range
    if low < 1:
        low = 1

    return range(low, middle + eps_range, step)


def _split_ids(clustered, mod_index):
    """
    Function for splitting ids based on a list with information about adherence to a cluster

    Parameters:
    :param clustered: dictionary of each used clust type and list with information about adherence to a cluster
    :param mod_index: list of read ids

    :return: dictionary of types of clustering and dictionary of clusters and split ids
    """

    ids_clustered = {}

    for clust_type in clustered:
        ids_clustered[clust_type] = {}
        cluster_list = clustered[clust_type]
        for cluster, r_name in zip(cluster_list, mod_index):
            if cluster not in ids_clustered[clust_type].keys():
                ids_clustered[clust_type][cluster] = []
            ids_clustered[clust_type][cluster].append(r_name)

    return ids_clustered


def _get_transformed(mod_table, imputer_divisor):
    """
    Function for normalizing dataset and imputing missing values (using KNNImputer)

    Parameters:
    :param mod_table:  pandas table of methylation on specific positions and reads
    :param imputer_divisor: divisor used for choosing amount of neighbours when imputing

    :return: numpy matrix of scaled and imputed values
    """

    transform_pipeline = make_pipeline(
        StandardScaler(),
        KNNImputer(n_neighbors=len(mod_table) // imputer_divisor),
    )

    return transform_pipeline.fit_transform(mod_table)


def _get_clustered(transformed_table, types, clusters, min_clusters, max_clusters,
                   hac_linkage, eps_iter, min_samples, neg_weight):
    """
    Function for clustering pre-processed epigenetic data.

    Parameters:
    :param transformed_table: numpy matrix with normalized and non-NaN data
    :param types: string with characters A,B,D,K defining which types of clustering will be used
    :param clusters: Number of clusters to be used in clustering (does not apply to DBSCAN)
    :param min_clusters: Minimum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param max_clusters: Maximum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param hac_linkage: Linkage to be used in HAC clustering
    :param eps_iter: iterable over eps values which should be used in DBSCAN clustering
    :param min_samples: minimum samples in DBSCAN cluster
    :param neg_weight: weight of outliers when choosing best DBSCAN clustering

    :return: dictionary of each used clust type and list with information about adherence to a cluster
    """

    result = {}

    if "K" in types:
        clusters = _find_cluster_count(clusters, min_clusters, max_clusters, transformed_table,
                                       lambda x:  KMeans(n_clusters=x))
        result[ClusterType.KMEANS] = KMeans(n_clusters=clusters).fit_predict(transformed_table)
        print(f"Clustered with kmeans, number of clusters found was: {clusters}, silhouette score was: "
              f"{silhouette_score(transformed_table, result[ClusterType.KMEANS])}")

    if "M" in types:
        clusters = _find_cluster_count(clusters, min_clusters, max_clusters, transformed_table,
                                       lambda x: MiniBatchKMeans(n_clusters=x))
        result[ClusterType.MINIK] = MiniBatchKMeans(n_clusters=clusters).fit_predict(transformed_table)
        print(f"Clustered with Mini batch k means, number of clusters found was: {clusters}, silhouette score was: "
              f"{silhouette_score(transformed_table, result[ClusterType.MINIK])}")

    if "A" in types:
        clusters = _find_cluster_count(clusters, min_clusters, max_clusters, transformed_table,
                                       lambda x: AgglomerativeClustering(n_clusters=x, compute_full_tree=True,
                                                                         linkage=hac_linkage))
        result[ClusterType.AGGLOMERATIVE] = AgglomerativeClustering(n_clusters=clusters, compute_full_tree=True,
                                                                    linkage=hac_linkage).fit_predict(transformed_table)
        print(f"Clustered with HAC, number of clusters found was: {clusters}, silhouette score was: "
              f"{silhouette_score(transformed_table, result[ClusterType.AGGLOMERATIVE])}")

    if "D" in types:
        result[ClusterType.DENSITY] = _density_cluster(transformed_table, eps_iter, min_samples, neg_weight)
        if 1 not in result[ClusterType.DENSITY]:
            print("No clusters were found with DBSCAN")
        else:
            outliers = np.count_nonzero(result[ClusterType.DENSITY] == -1)
            clusters = len(set(result[ClusterType.DENSITY]))
            if outliers > 0:
                clusters -= 1
            print(f"Clustered with DBSCAN, found {clusters} clusters,"
                  f"and {outliers} outlier(s),"
                  f" silhouette score was: {silhouette_score(transformed_table, result[ClusterType.DENSITY])}")

    return result


def _density_cluster(transformed_table, eps_iter, min_samples, neg_weight):
    """
    Function for clustering data with DBSCAN multiple times(only best result returned)

    Parameters:
    :param transformed_table: numpy matrix with normalized and non-NaN data
    :param eps_iter: iterable over eps values which should be used
    :param min_samples: minimum samples
    :param neg_weight: weight of outliers when choosing best DBSCAN clustering

    :return: list with
    """

    results = []

    for eps in eps_iter:
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        dbscan.fit(transformed_table)
        results.append(dbscan.labels_)

    return _choose_best_density(results, neg_weight, transformed_table)


def _choose_best_density(results, neg_weight, transformed_table):
    """
    Function for deciding the best density clustering

    Parameters:
    :param results: list of list to be chosen from
    :param neg_weight: weight of outliers
    :param transformed_table: numpy matrix with normalized and non-NaN data

    :return: chosen clustering (list of values indicative of belonging to specific cluster)
    """

    best = None
    score = -np.inf

    for possibility in results:
        temp_score = _count_density_score(possibility, neg_weight, transformed_table)
        if temp_score >= score:
            score = temp_score
            best = possibility

    return best


def _count_density_score(possibility, neg_weight, transformed_table):
    """
    Function for calculating score of clustering

    ----------
    -the less outliers the better
    ----------

    Parameters:
    :param possibility: list of values indicative of belonging to specific cluster
    :param neg_weight: weight of outliers
    :param transformed_table: numpy matrix with normalized and non-NaN data

    :return: silhouette score with outliers counted in
    """

    unique, counts = np.unique(possibility, return_counts=True)
    dict_counts = dict(zip(unique, counts))

    if len(dict_counts) <= 1:
        return -np.inf

    negative = dict_counts.get(-1, 0)/len(possibility) * neg_weight
    return silhouette_score(transformed_table, possibility) - negative


def _find_cluster_count(clusters, min_clusters, max_clusters, transformed_table, clustering_f):
    """
    Finds the best number of clusters for a given clustering function, uses silhouette score as metric

    :param clusters: Number of clusters to be used in clustering (does not apply to DBSCAN)
    :param min_clusters: Minimum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param max_clusters: Maximum number of clusters to be used in clustering (does not apply to DBSCAN)
    :param transformed_table: numpy matrix with normalized and non-NaN data
    :param clustering_f: function to be used for clustering, that has one parameter number of clusters

    :return: count of clusters based on the best silhouette score
    """

    score = -np.inf
    if not clusters:
        for i in range(min_clusters, max_clusters + 1):
            s_score = silhouette_score(transformed_table, clustering_f(i).fit_predict(transformed_table))
            if s_score > score:
                score = s_score
                clusters = i
    return clusters

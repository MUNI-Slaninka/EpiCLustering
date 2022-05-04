from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import Birch

from enum import Enum, auto


class ClusterType(Enum):
    KMEANS = auto()
    BIRCH = auto()
    AGGLOMERATIVE = auto()
    DENSITY = auto()


def clustering(mod_table, types=None, clusters=2, eps_list=None, min_samples=None):
    print("Clustering")

    if not types:
        types = [ClusterType.KMEANS, ClusterType.BIRCH, ClusterType.AGGLOMERATIVE, ClusterType.DENSITY]
    if not min_samples:
        min_samples = mod_table.shape[0] // 5
    if not eps_list:
        eps_list = _get_eps_list(mod_table)

    clustered = _get_clustered(_get_transformed(mod_table), types, clusters, eps_list, min_samples)

    return _split_ids(clustered, mod_table)


def _get_eps_list(mod_table):
    return 0
    # todo


def _split_ids(clustered, mod_table):
    return {}
    # todo


def _get_transformed(mod_table):
    transform_pipeline = make_pipeline(
        StandardScaler(),
        KNNImputer(n_neighbors=len(mod_table) // 4),
    )

    return transform_pipeline.fit_transform(mod_table)


def _get_clustered(transformed_table, types, clusters, eps_list, min_samples):
    result = {}

    if ClusterType.KMEANS in types:
        result[ClusterType.KMEANS] = KMeans(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.BIRCH in types:
        result[ClusterType.BIRCH] = Birch(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.AGGLOMERATIVE in types:
        result[ClusterType.AGGLOMERATIVE] = AgglomerativeClustering(n_clusters=clusters).fit_predict(transformed_table)

    if ClusterType.DENSITY in types:
        result[ClusterType.DENSITY] = __density_cluster(transformed_table, eps_list, min_samples)

    return result


def __density_cluster(transformed_table, eps_list, min_samples):
    result = []

    for eps in eps_list:
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        dbscan.fit(transformed_table)
        result.append(dbscan.labels_)

    return __choose_best_density(result)


def __choose_best_density(result):
    return []
    # todo

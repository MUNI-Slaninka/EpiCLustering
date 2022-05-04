from read_clustering import clustering
from read_processing import get_processed


if __name__ == '__main__':
    clustering(get_processed("../mod_bases_mapped.bam"))

import os
import pysam


def output(dict_clusters, file_path, output_dir, name, start, end, cg_mod_table, basic_mod_table, bam, positions):
    """
    Function for creating all necessary output files

    Parameters:
    :param dict_clusters: dictionary consisting of clusters, reads and type of clustering
    :param file_path: path to input bam file
    :param output_dir: baseline for path
    :param name: name of region from reads to be used
    :param start: start position to be used
    :param end: end position to be used
    :param cg_mod_table: pandas table of CG positions and their modifications for each read
    :param basic_mod_table: pandas table of C positions and their modifications for each read
    :param bam: True if bam files should be outputted
    :param positions: True if csv with pandas frames for positions should be outputted
    :return:
    """

    if positions:
        os.mkdir(output_dir)
        cg_mod_table.to_csv(f"{output_dir}/cg_positions.csv")
        basic_mod_table.to_csv(f"{output_dir}/c_positions.csv")
    if bam:
        _output_bam(dict_clusters, file_path, output_dir, name, start, end)


def _output_bam(dict_clusters, file_path, output_dir, name, start, end):
    """
    Function for saving bam files of clustered reads by their clustering

    Parameters:
    :param dict_clusters: dictionary consisting of clusters, reads and type of clustering
    :param file_path: path to input bam file
    :param output_dir: baseline for path
    :param name: name of region from reads to be used
    :param start: start position to be used
    :param end: end position to be used

    :return: None
    """

    with pysam.AlignmentFile(file_path, "rb") as bam_file:
        for dir_1 in dict_clusters:
            for dir_2 in dict_clusters[dir_1]:
                os.makedirs(f"{output_dir}/{dir_1}/{dir_2}")
        for key in dict_clusters:
            for read in bam_file.fetch(name, start, end):
                with pysam.AlignmentFile(__get_read_path(read.qname, dict_clusters, key, output_dir),
                                         "wb", template=bam_file) as out_file:
                    out_file.write(read)


def __get_read_path(read_id, dict_clusters, key, output_dir):
    """
    Function for choosing the right path to store the bam file of a read

    Parameters:
    :param read_id: id of a read
    :param dict_clusters: dictionary consisting of clusters, reads and type of clustering
    :param key: type of clustering
    :param output_dir: baseline for path

    :return: path to store the bam file
    """

    for k2, ids in dict_clusters[key].items():
        if read_id in ids:
            return f"{output_dir}/{key}/{k2}/{read_id}.bam"

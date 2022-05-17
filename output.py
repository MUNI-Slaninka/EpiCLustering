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

    dict_reads={}
    with pysam.AlignmentFile(file_path, "rb") as bam_file:
        for key in dict_clusters:
            for read in bam_file.fetch(name, start, end):
                path = _get_read_path(read.qname, dict_clusters, key, output_dir)
                dict_reads[path] = dict_reads.get(path, []) + [read]
        for dir_name in dict_clusters:
            os.makedirs(f"{output_dir}/{dir_name}")
        for path in dict_reads:
            out_file = pysam.AlignmentFile(path, "wb", template=bam_file)
            for read in dict_reads[path]:
                out_file.write(read)
            out_file.close()


def _get_read_path(read_id, dict_clusters, key, output_dir):
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
            return f"{output_dir}/{key}/{k2}.bam"

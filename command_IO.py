import os
import pysam






def output_bam(dict_clusters, bam_file, output_dir="../result_bam", name='chrX_v0.7', start=72280000, end=72282000):
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    for dir_name in dict_clusters:
        os.makedirs(f"{output_dir}/{dir_name}")
        for file_name in dict_clusters[dir_name]:
            with open(f"{output_dir}/{dir_name}/{file_name}.bam", mode='w'):
                pass
    for read in bam_file.fetch(name, start, end):
        print(__get_read_path(read.qname, dict_clusters, output_dir))
        with pysam.AlignmentFile(__get_read_path(read.qname, dict_clusters, output_dir), "wb", template=bam_file)\
                as out_file:
            out_file.write(read)


def __get_read_path(read_id, dict_clusters, output_dir):
    for k1, d in dict_clusters.items():
        for k2, ids in d.items():
            if read_id in ids:
                return f"{output_dir}/{k1}/{k2}.bam"

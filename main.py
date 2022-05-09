from read_clustering import clustering
from read_processing import get_processed
from command_IO import output_bam

import os
import sys
import argparse


def parse_arguments():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument()
    args = my_parser.parse_args()

    input_path = args.Path


if not os.path.isdir(input_path):
    print('The path specified does not exist')
    sys.exit()
if __name__ == '__main__':
    clust_dict = clustering(get_processed("../mod_bases_mapped.bam"))
    print(clust_dict)
    output_bam(clust_dict, "../mod_bases_mapped.bam")

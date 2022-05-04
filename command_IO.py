import os
import sys
import argparse


def parse_arguments():
    # todo
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument()
    args = my_parser.parse_args()

    input_path = args.Path

    if not os.path.isdir(input_path):
        print('The path specified does not exist')
        sys.exit()


def output_bam():
    pass
    #todo

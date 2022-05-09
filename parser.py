import os
import sys
import argparse


def parse_arguments():
    """
    Parses command line arguments

    Formats help

    :return: parsed arguments
    """

    arg_p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Tool for clustering reads with the use of epigenetic data.",
                                    epilog="Bam files should be mapped and have Mm and M1 tags filled."
                                           "Cluster types should be specified by string containing letters "
                                           "A, B, D and K, thus choosing clustering by "
                                           "Agglomeration, Birch, Kmeans and DBScan"
                                    )

    arg_p.add_argument("-i", "--input", required=True,
                       help="Specify path to input BAM file")

    arg_p.add_argument("-o", "--out-dir", required=True,
                       help="Specify path to directory in which output will be stored")

    arg_p.add_argument("-n", "--ref-name", required=True,
                       help="Specify path to input BAM file")

    arg_p.add_argument("-s", "--start", required=True, type=int,
                       help="Specify path to directory in which output will be stored")

    arg_p.add_argument("-e", "--end", required=True, type=int,
                       help="Specify path to input BAM file")

    arg_p.add_argument("-t", "--clust-types", required=False, default="ABDK",
                       help="Specify types of clustering algorithms")

    arg_p.add_argument("-d", "--imputer-divisor", required=False, type=int, default=4,
                       help="Specify types of clustering algorithms")

    arg_p.add_argument("-c", "--clusters", required=False, type=int, default=2,
                       help="Specify number of expected clusters (for non-Density based clustering)")

    arg_p.add_argument("-m", "--min-samples", required=False, type=int,
                       help="Specify number of minimal samples in cluster (for Density based clustering)")

    arg_p.add_argument("--min-divisor", required=False, type=int, default=5,
                       help="Specify divisor for automatic minimum samples calculation (for Density based clustering)")

    arg_p.add_argument("--eps-divisor", required=False, type=int, default=4,
                       help="Specify divisor for automatic epsilon value calculation (for Density based clustering)")

    arg_p.add_argument("--eps-step", required=False, type=int, default=1,
                       help="Specify step size for automatic epsilon value calculation (for Density based clustering)")

    arg_p.add_argument("--eps-range", required=False, type=int, default=10,
                       help="Specify broadness of epsilon range for automatic "
                            "epsilon value calculation (for Density based clustering)")

    arg_p.add_argument("--neg-weight", required=False, type=int, default=3,
                       help="Specify weight of outliers for choice of the best Density based clustering")

    arg_p.add_argument("--out-bam", required=False, type=bool, default=True,
                       help="Specify if clustering results will be output as bam file")

    arg_p.add_argument("--out-positions", required=False, type=bool, default=False,
                       help="Specify if positions used for clustering will be output as csv file")

    args = vars(arg_p.parse_args())

    if not os.path.isfile(args['input']):
        print('The path specified does not exist')
        sys.exit()

    if os.path.isdir(args['out_dir']):
        print('Output dir already exists')
        sys.exit()

    return args



from clustering import clustering
from output import output
from input_processing import get_processed
from parser import parse_arguments


def main():
    """
    main loop

    :return: None
    """

    args = parse_arguments()

    cg_mod_table, basic_mod_table = get_processed(args["input"], args["ref_name"], args["start"], args["end"])

    clust_dict = clustering(cg_mod_table, args["clust_types"], args["imputer_divisor"], args["clusters"],
                            args["min_clusters"], args["max_clusters"], args["hac_linkage"],
                            args["min_samples"], args["min_divisor"], args["eps_divisor"],
                            args["eps_step"], args["eps_range"], args["neg_weight"])

    output(clust_dict, args["input"], args["out_dir"], args["ref_name"], args["start"], args["end"],
           cg_mod_table, basic_mod_table, args["out_bam"], args["out_positions"])


if __name__ == '__main__':
    main()

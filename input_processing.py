import pandas as pd
import pysam

from itertools import chain
from warnings import simplefilter


def get_processed(bam_file, name, start, end):
    """
    Function for retrieving  methylation per CG position from bam file

    Parameters:
    :param bam_file: path to bam file
    :param name: name of region from reads to be used
    :param start: start position to be used
    :param end: end position to be used

    :return: Tuple(pandas table of CG positions and their modifications for each read,
                   pandas table of C positions and their modifications for each read)
    """

    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    dict_mod_pos, dict_cit_pos = _get_pos_cit(bam_file, name, start, end)
    mod_table = _create_mod_table(dict_mod_pos, dict_cit_pos)

    return _get_cg_pairs(mod_table), mod_table


def _get_pos_cit(bam_file, name, start, end):
    """
    Function for retrieving modified and cytosine positions reads from bam file

    Parameters:
    :param bam_file: path to bam_file
    :param name: name of region from reads to be used
    :param start: start position to be used
    :param end: end position to be used

    :return: tuple:(dictionary of modified positions per read, dictionary of cytosine positions per read)
    """

    print("Processing reads")

    dict_mod_pos = {}
    dict_cit_pos = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        for read in bam_file.fetch(name, start, end):
            dict_mod_pos[read.qname] = _process_read(read, start, end, dict_cit_pos)

    return dict_mod_pos, dict_cit_pos


def _process_read(read, start, end, dict_cit_pos):
    """
    Function for retrieving modified and cytosine positions from single read

    Parameters:
    :param read: single read
    :param start: start position to be used
    :param end: end position to be used
    :param dict_cit_pos: dictionary of cytosine positions per read !input-output parameter

    :return:
    """

    mod_bases = []
    ref_pos_list = read.get_reference_positions(full_length=True)
    mod_bases_list = list(read.modified_bases.values())

    if mod_bases_list:
        mod_bases = mod_bases_list[0]
        mod_bases = [(ref_pos_list[x[0]], x[1]) for x in mod_bases if
                     ref_pos_list[x[0]] and start < ref_pos_list[x[0]] < end]

    _fill_other_cit(read, dict_cit_pos)

    return mod_bases


def _fill_other_cit(read, dict_cit_pos):
    """
    Function for filling dictionary of cytosine positions with non modified cytosines

    Parameters:
    :param read: single read
    :param dict_cit_pos: dictionary of cytosine positions per read !input-output parameter

    :return: None
    """

    cit_pos = []
    sequence = read.query_sequence
    ref_pos_list = read.get_reference_positions(full_length=True)

    if sequence:
        for i in range(len(sequence)):
            if (sequence[i] == "C" and read.is_forward) or (sequence[i] == "G" and read.is_reverse):
                cit_pos.append(ref_pos_list[i])
        dict_cit_pos[read.qname] = cit_pos


def _create_mod_table(dict_mod_pos, dict_cit_pos):
    """
    Function for creating a pandas dataframe with methylation for each read and cytosine position

    Parameters:
    :param dict_mod_pos: dictionary of modified positions per read
    :param dict_cit_pos: dictionary of cytosine positions per read

    :return: pandas dataframe with methylation for each read and cytosine position
    """

    print("Creating table of modified bases")

    columns = sorted(set(map(lambda x: x[0], chain(*dict_mod_pos.values()))))
    mod_table = pd.DataFrame(columns=columns, index=dict_mod_pos.keys())

    for row_index in mod_table.index:
        print(".", end="")
        for column_index in mod_table.columns:
            if column_index in dict_cit_pos[row_index]:
                mod_table.at[row_index, column_index] = 1

    for read_id in dict_mod_pos:
        for position, prob in dict_mod_pos[read_id]:
            mod_table.at[read_id, position] = prob

    print()

    return mod_table


def _get_cg_pairs(mod_table):
    """
    Function for filtering and combining pandas data frame based on CG pairs

    Parameters:
    :param mod_table: pandas data frame with per position cytosine modifications

    :return: pandas data frame with modifications per CG pair
    """

    print("Filtering for CG pairs")

    pair_cg = [x for x in mod_table.columns if x+1 in mod_table.columns]
    pair_cg_table = pd.DataFrame()

    for pos in pair_cg:
        combined = mod_table[pos].combine_first(mod_table[pos + 1])
        if combined.count():
            pair_cg_table[pos] = combined

    pair_cg_table = pair_cg_table.copy()

    return pair_cg_table

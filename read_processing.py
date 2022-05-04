import pysam
import random
import pandas as pd

from itertools import chain
from warnings import simplefilter


def get_processed(bam_file, name='chrX_v0.7', start=72280000, end=72282000):
    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    dict_mod_pos, dict_cit_pos = _get_pos_cit(bam_file, name, start, end)
    mod_table = _create_mod_table(dict_mod_pos, dict_cit_pos)

    return _get_cg_pairs(mod_table)


def _get_pos_cit(bam_file, name, start, end):
    print("Processing reads")

    dict_mod_pos = {}
    dict_cit_pos = {}
    bam_file = pysam.AlignmentFile(bam_file, "rb")

    for read in bam_file.fetch(name, start, end):
        dict_mod_pos[read.qname] = __process_read(read, start, end, dict_cit_pos)

    bam_file.close()
    return dict_mod_pos, dict_cit_pos


def __process_read(read, start, end, dict_cit_pos):
    mod_bases = []
    ref_pos_list = read.get_reference_positions(full_length=True)
    mod_bases_list = list(read.modified_bases.values())

    if mod_bases_list:
        mod_bases = mod_bases_list[0]
        mod_bases = [(ref_pos_list[x[0]], x[1]) for x in mod_bases if
                     ref_pos_list[x[0]] and start < ref_pos_list[x[0]] < end]

    __fill_other_cit(read, dict_cit_pos)

    return mod_bases


def __fill_other_cit(read, dict_cit_pos):
    cit_pos = []
    sequence = read.query_sequence
    ref_pos_list = read.get_reference_positions(full_length=True)

    if sequence:
        for i in range(len(sequence)):
            if (sequence[i] == "C" and read.is_forward) or (sequence[i] == "G" and read.is_reverse):
                cit_pos.append(ref_pos_list[i])
        dict_cit_pos[read.qname] = cit_pos


def _create_mod_table(dict_mod_pos, dict_cit_pos):
    print("Creating table of modified bases")

    columns = sorted(set(map(lambda x: x[0], chain(*dict_mod_pos.values()))))
    mod_table = pd.DataFrame(columns=columns, index=dict_mod_pos.keys())

    for row_index in mod_table.index:
        print(".", end="")
        for column_index in mod_table.columns:
            if column_index in dict_cit_pos[row_index]:
                mod_table.at[row_index, column_index] = random.randint(1, 16)

    for read_id in dict_mod_pos:
        for position, prob in dict_mod_pos[read_id]:
            mod_table.at[read_id, position] = prob

    print()
    return mod_table


def _get_cg_pairs(mod_table):
    print("Filtering for CG pairs")

    pair_cg = [x for x in mod_table.columns if x+1 in mod_table.columns]
    pair_cg_table = pd.DataFrame()

    for pos in pair_cg:
        combined = mod_table[pos].combine_first(mod_table[pos + 1])
        if combined.count():
            pair_cg_table[pos] = combined

    pair_cg_table = pair_cg_table.copy()
    return pair_cg_table

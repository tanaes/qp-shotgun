# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# This file contains functions related to generating default parameters for
# Shogun
# ------------------------------------------------------------------------------

import os
from os.path import join, isdir

ALIGNERS = ["utree", "burst", "bowtie2"]


def get_dbs(db_folder):
    dbs = {}
    # Loop through the databases and create a dict of them
    for folder in os.listdir(db_folder):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs[folder] = folder_path

    return(dbs)


def get_dbs_list(db_folder):
    dbs = []
    # Loop through the databases and create a list string
    for folder in sorted(os.listdir(db_folder)):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            dbs.append(folder_path)
    dbs_formatted = (', '.join('"' + item + '"' for item in dbs))

    return(dbs_formatted)


def generate_shogun_dflt_params():
    dflt_param_set = {}
    db_parent_path = os.environ["QC_SHOGUN_DB_DP"]
    # Get a the databases available and the database name
    dbs = get_dbs(db_parent_path)
    # Create dict with command options per database
    for db in dbs:
        for aligner in ALIGNERS:
            dflt_param_set[db+'_'+aligner] = {'Database': dbs[db],
                                              'Aligner tool': aligner,
                                              'Number of threads': 1}

    return(dflt_param_set)


def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

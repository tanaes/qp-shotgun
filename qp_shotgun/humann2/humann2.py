# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import basename, join

from future.utils import viewitems
import pandas as pd


def generate_humann2_analysis_commands(forward_seqs, reverse_seqs, map_file,
                                       out_dir, parameters):
    """Generates the humann2 commands

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    list of str
        The humann2 commands

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    The forward and reverse files are going to have different filenames but
    the sample name is going to be the same so the results are merged when
    joining the outputs
    """
    # making sure the forward and reverse reads are in the same order
    forward_seqs.sort()
    if reverse_seqs:
        if len(forward_seqs) != len(reverse_seqs):
            raise ValueError('Your reverse and forward files are of different '
                             'length. Forward: %s. Reverse: %s.' %
                             (', '.join(forward_seqs),
                              ', '.join(reverse_seqs)))
        reverse_seqs.sort()

    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # we match sample name and forward filename
    samples = []
    for i, fname in enumerate(forward_seqs):
        f = basename(fname)
        # removing extentions: fastq or fastq.gz
        if 'fastq' in f.lower().rsplit('.', 2):
            f = f[:f.lower().rindex('.fastq')]
        # this try/except block is simply to retrieve all possible errors
        # and display them in the next if block
        try:
            samples.append((fname, f, sn_by_rp[f]))
            if reverse_seqs:
                fr = basename(reverse_seqs[i])
                if 'fastq' in fr.lower().rsplit('.', 2):
                    fr = fr[:fr.lower().rindex('.fastq')]
                samples.append((reverse_seqs[i], fr, sn_by_rp[f]))
            del sn_by_rp[f]
        except KeyError:
                pass

    if sn_by_rp:
        raise ValueError(
            'Some run_prefix values do not match your sample names: %s'
            % ', '.join(sn_by_rp.keys()))

    cmds = []
    params = ' '.join(["%s %s" % (k, v) for k, v in viewitems(parameters)])
    for ffn, fn, s in samples:
        cmds.append('humann2 --input %s --output %s --output-basename %s '
                    '--output-format biom %s' % (ffn, join(out_dir, fn), s,
                                                 params))

    return cmds


def get_sample_names_by_run_prefix(mapping_file):
    """Generates a dictionary of run_prefix and sample names

    Parameters
    ----------
    mapping_file : str
        The mapping file

    Returns
    -------
    dict
        Dict mapping run_prefix to sample id

    Raises
    ------
    ValueError
        If there is more than 1 sample per run_prefix
    """
    qiime_map = pd.read_csv(mapping_file, delimiter='\t', dtype=str,
                            encoding='utf-8')
    qiime_map.set_index('#SampleID', inplace=True)

    samples = {}
    errors = []
    for prefix, df in qiime_map.groupby('run_prefix'):
        len_df = len(df)
        if len_df != 1:
            errors.append('%s has %d samples (%s)' % (prefix, len_df,
                                                      ', '.join(df.index)))
        else:
            samples[prefix] = df.index.values[0]

    if errors:
        raise ValueError("You have run_prefix values with multiple "
                         "samples: %s" % ' -- '.join(errors))

    return samples

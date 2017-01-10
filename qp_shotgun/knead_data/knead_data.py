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
                            encoding='utf-8', keep_default_na=False,
                            na_values=[])
    qiime_map.set_index('#SampleID', inplace=True)

    samples = {}
    errors = []
    for prefix, df in qiime_map.groupby('run_prefix'):
        len_df = len(df)
        if len_df != 1:
            errors.append('%s has %d samples (%s)' % (prefix, len_df,
                                                      ', '.join(df.index)))
        else:
            samples[prefix] = df.index[0]

    if errors:
        raise ValueError("You have run_prefix values with multiple "
                         "samples: %s" % ' -- '.join(errors))

    return samples

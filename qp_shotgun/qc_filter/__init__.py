# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .qc_filter import qc_filter

__all__ = ['qc_filter']

# Define the qc_filter command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'Bowtie2 database to filter': ['choice:["Human"]','Human'],
    'Number of threads to be used': ['integer', '4']
    }
outputs = {'Filtered files': 'per_sample_FASTQ'}
dflt_param_set = {
    'Human Filtering': {
        'Bowtie2 database to filter': 'Human',
        'Number of threads to be used': 4
        }
}

qc_filter_cmd = QiitaCommand(
    'QC_Filter', "Sequence QC - Filtering", qc_filter,
    req_params, opt_params, outputs, dflt_param_set)

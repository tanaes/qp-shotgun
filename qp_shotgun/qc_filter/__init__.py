# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand
from .qc_filter import (qc_filter, generate_qc_filter_dflt_params,
                        get_dbs_list)
from os.path import join
from os import environ
__all__ = ['qc_filter']

# Define the qc_filter command
default_db = join(environ["QC_FILTER_DB_DP"], 'phix', 'phix')
default_db_list = get_dbs_list(environ["QC_FILTER_DB_DP"])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'Bowtie2 database to filter': ["choice: [%s]" % default_db_list,
                                   default_db],
    'Number of threads to be used': ['integer', '4']
    }
outputs = {'Filtered files': 'per_sample_FASTQ'}
dflt_param_set = generate_qc_filter_dflt_params()

qc_filter_cmd = QiitaCommand(
    'QC_Filter', "Sequence QC - Filtering", qc_filter,
    req_params, opt_params, outputs, dflt_param_set)

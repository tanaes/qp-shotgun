# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand
from .shogun import shogun
from .utils import (generate_shogun_dflt_params, get_dbs_list)
from os.path import join
from os import environ


__all__ = ['Shogun']

# Define the shogun command
default_db = join(environ["QC_SHOGUN_DB_DP"], 'ref82')
default_db_list = get_dbs_list(environ["QC_SHOGUN_DB_DP"])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ["choice: [%s]" % default_db_list,
                 default_db],
    # aligner
    'Aligner tool': ['choice:["utree", "burst", "bowtie2"]', 'bowtie2'],
    # threads
    'Number of threads': ['integer', '5'],
    }
outputs = {'Shogun Alignment Profile': 'BIOM'}
dflt_param_set = generate_shogun_dflt_params()

shogun_cmd = QiitaCommand(
    'Shogun', "Functional and Taxonomic Predictions", shogun,
    req_params, opt_params, outputs, dflt_param_set)

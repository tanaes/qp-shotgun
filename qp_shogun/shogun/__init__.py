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
default_db = join(environ["QC_SHOGUN_DB_DP"], 'shogun')
default_db_list = get_dbs_list(environ["QC_SHOGUN_DB_DP"])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ["choice: [%s]" % default_db_list,
                 default_db],
    # aligner
    'Aligner tool': ['choice:["utree", "burst", "bowtie2"]', 'bowtie2'],
    # threads
    'Number of threads': ['integer', '1'],
    }
outputs = {
    'Functional Predictions - genus, KEGG Modules Coverage': 'BIOM',
    'Functional Predictions - genus, KEGG Modules': 'BIOM',
    'Functional Predictions - genus, KEGG Pathways Coverage': 'BIOM',
    'Functional Predictions - genus, KEGG Pathways': 'BIOM',
    'Functional Predictions - genus, KEGG': 'BIOM',
    'Functional Predictions - genus, Normalized': 'BIOM',
    'Functional Predictions - species, KEGG Modules Coverage': 'BIOM',
    'Functional Predictions - species, KEGG Modules': 'BIOM',
    'Functional Predictions - species, KEGG Pathways Coverage': 'BIOM',
    'Functional Predictions - species, KEGG Pathways': 'BIOM',
    'Functional Predictions - species, KEGG': 'BIOM',
    'Functional Predictions - species, Normalized': 'BIOM',
    'Functional Predictions - strain, KEGG Modules Coverage': 'BIOM',
    'Functional Predictions - strain, KEGG Modules': 'BIOM',
    'Functional Predictions - strain, KEGG Pathways Coverage': 'BIOM',
    'Functional Predictions - strain, KEGG Pathways': 'BIOM',
    'Functional Predictions - strain, KEGG': 'BIOM',
    'Functional Predictions - strain, Normalized': 'BIOM',
    'Taxonomic Predictions - genus': 'BIOM',
    'Taxonomic Predictions - species': 'BIOM',
    'Taxonomic Predictions - strain': 'BIOM'}
dflt_param_set = generate_shogun_dflt_params()

shogun_cmd = QiitaCommand(
    'Shogun', "Functional and Taxonomic Predictions", shogun,
    req_params, opt_params, outputs, dflt_param_set)

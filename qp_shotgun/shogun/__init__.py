# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .shogun import shogun

__all__ = ['qc_shogun']

# Define the qc_trim command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ['string', ''],
    # aligner
    'Aligner tool': ['choice:["utree", "burst", "bowtie2"]', 'bowtie2'],
    # taxonomic levels
    'Taxonomy Level': ['choice:["kingdom", "phylum", "class", "order", '
                       '"family", "genus", "species", "strain", "all"]',
                       'all'],
    # threads
    'Number of threads': ['integer', '1'],
    }
outputs = {'Functional Predictions': 'BIOM', 'Taxonomic Predictions': 'BIOM'}
dflt_param_set = {
    'Shogun': {
        'Database': '',
        'Aligner tool': 'bowtie2',
        'Taxonomy Level': 'all',
        'Number of threads': 1,
        }
}

shogun = QiitaCommand(
    'Shogun', "Functional and Taxonomic Predictions", shogun,
    req_params, opt_params, outputs, dflt_param_set)

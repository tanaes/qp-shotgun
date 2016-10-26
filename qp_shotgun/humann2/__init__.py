# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .humann2 import humann2

__all__ = ['humann2']

# Define the HUMAnN2 command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # there are other parameters not included that will be ignored in this
    # configuration as we assume that the correct values were set in the env
    # by the admin installing the tools:
    # metaphlan
    # bowtie2
    # usearch
    # rapsearch
    # diamond
    # id-mapping
    # pathways-database
    # o-log
    # input-format
    # search-mode
    'nucleotide-database': ['choice:["default"]', 'default'],
    'protein-database': ['choice:["default"]', 'default'],
    'bypass-prescreen': ['boolean', 'False'],
    'bypass-nucleotide-index': ['boolean', 'False'],
    'bypass-translated-search': ['boolean', 'False'],
    'bypass-nucleotide-search': ['boolean', 'False'],
    'annotation-gene-index': ['integer', '8'],
    'evalue': ['float', '1.0'],
    'metaphlan-options': ['string', '-t rel_ab'],
    'log-level': ['choice:["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]',
                  'DEBUG'],
    'remove-temp-output': ['boolean', 'False'],
    'threads': ['integer', '1'],
    'prescreen-threshold': ['float', '0.01'],
    'identity-threshold': ['float', '50.0'],
    'translated-subject-coverage-threshold': ['float', '50.0'],
    'translated-query-coverage-threshold': ['float', '90.0'],
    'translated-alignment': ['choice:["usearch", "rapsearch", "diamond"]',
                             'diamond'],
    'xipe': ['choice:["on", "off"]', 'off'],
    'minpath': ['choice:["on", "off"]', 'on'],
    'pick-frames': ['choice:["on", "off"]', 'off'],
    'gap-fill': ['choice:["on", "off"]', 'off'],
    'output-format': ['choice:["tsv", "biom"]', 'biom'],
    'output-max-decimals': ['integer', '10'],
    'remove-stratified-output': ['boolean', 'False'],
    'pathways': ['choice:["metacyc", "unipathway"]', 'metacyc'],
    'memory-use': ['choice:["minimum", "maximum"]', 'minimum'],
    'remove-column-description-output': ['boolean', 'True']}
outputs = {
    'Gene family table': 'BIOM',
    'Path coverage table': 'BIOM',
    'Path abundance table': 'BIOM',
    'Gene family CMP table': 'BIOM',
    'Path coverage RELAB table': 'BIOM',
    'Path abundance RELAB table': 'BIOM',
    'Gene family CMP table - stratified': 'BIOM',
    'Path coverage RELAB table - stratified': 'BIOM',
    'Path abundance RELAB table - stratified': 'BIOM',
    'Gene family CMP table - unstratified': 'BIOM',
    'Path coverage RELAB table - unstratified': 'BIOM',
    'Path abundance RELAB table - unstratified': 'BIOM'}
dflt_param_set = {
    'Defaults': {
        'nucleotide-database': 'default', 'protein-database': 'default',
        'bypass-prescreen': False, 'bypass-nucleotide-index': False,
        'bypass-translated-search': False, 'bypass-nucleotide-search': False,
        'annotation-gene-index': 8, 'evalue': 1.0, 
        'metaphlan-options': '-t rel_ab', 'log-level': 'DEBUG',
        'remove-temp-output': False, 'threads': 1,
        'prescreen-threshold': 0.01, 'identity-threshold': 50.0,
        'translated-subject-coverage-threshold': 50.0,
        'translated-query-coverage-threshold': 90.0,
        'translated-alignment': 'diamond',
        'xipe': 'off', 'minpath': 'on', 'pick-frames': 'off',
        'gap-fill': 'off', 'output-format': 'biom',
        'output-max-decimals': 10, 'remove-stratified-output': False,
        'pathways': 'metacyc', 'memory-use': 'minimum'}
}
humann2_cmd = QiitaCommand(
    "HUMAnN2 0.9.1", "Community profiling", humann2, req_params, opt_params,
    outputs, dflt_param_set)

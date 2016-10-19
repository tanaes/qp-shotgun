# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .humann2 import humann2

__all__ = ['humann2']


# Initialize the plugin
plugin = QiitaPlugin(
    'HUMAnN2', '0.9.1', 'HUMAnN2 is the next generation of HUMAnN (HMP '
    'Unified Metabolic Analysis Network)')

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
    'nucleotide-database': ['choice:["chocophlan"]', 'chocophlan'],
    'protein-database': ['choice:["uniref"]', 'uniref'],
    'bypass-prescreen': ['boolean', 'False'],
    'bypass-nucleotide-index': ['boolean', 'False'],
    'bypass-translated-search': ['boolean', 'False'],
    'bypass-nucleotide-search': ['boolean', 'False'],
    'annotation-gene-index': ['integer', '8'],
    'evalue': ['float', '1.0'],
    'search-mode': ['choice:["", "uniref50", "uniref90"]', ''],
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
    'input-format': ['choice:["", "fastq", "fastq.gz", "fasta", "fasta.gz", '
                     '"sam", "bam", "blastm8", "genetable", "biom"]', ''],
    'pathways': ['choice:["metacyc", "unipathway"]', 'metacyc'],
    'memory-use': ['choice:["minimum", "maximum"]', 'minimum'],
    'remove-column-description-output': ['boolean', 'True']}
outputs = {'BIOM': 'HUMAnN2 output'}
dflt_param_set = {
    'Defaults': {
        'nucleotide-database': 'chocophlan', 'protein-database': 'uniref',
        'bypass-prescreen': False, 'bypass-nucleotide-index': False,
        'bypass-translated-search': False, 'bypass-nucleotide-search': False,
        'annotation-gene-index': 8, 'evalue': 1.0, 'search-mode': '',
        'metaphlan-options': '-t rel_ab', 'log-level': 'DEBUG',
        'remove-temp-output': False, 'threads': 1,
        'prescreen-threshold': 0.01, 'identity-threshold': 50.0,
        'translated-subject-coverage-threshold': 50.0,
        'translated-query-coverage-threshold': 90.0,
        'translated-alignment': 'diamond',
        'xipe': 'off', 'minpath': 'on', 'pick-frames': 'off',
        'gap-fill': 'off', 'output-format': 'biom',
        'output-max-decimals': 10, 'remove-stratified-output': False,
        'input-format': '', 'pathways': 'metacyc', 'memory-use': 'minimum'}
}
humann2_cmd = QiitaCommand(
    "HUMAnN2", "Community profiling", humann2, req_params, opt_params,
    outputs, dflt_param_set)
plugin.register_command(humann2_cmd)

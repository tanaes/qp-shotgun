# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

#from .fastqc import fastqc

__all__ = ['fastqc']

# Initialize the plugin
plugin = QiitaPlugin(
    'FastQC', '0.11.5', 'FastQC is a quality control tool for high throughput'
    ' sequence data.')

# Define the FastQC command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # there are other parameters not included that will be ignored in this
    # configuration as we assume that the correct values were set in the env
    # by the admin installing the tools:
    # --input # input FASTQ file (add a second argument for paired input files)
    # --outdir # directory to write output files
    'extract': ['boolean', 'False'],  # uncompress zipped output file
    'noextract': ['boolean', 'True'],  # compress zipped output file
    'threads': ['integer', '1'],  # threads to run
    'kmers': ['integer', '7'],  # kmer content module kmer length
    }
outputs = {'html_summary': 'html_summary'}
dflt_param_set = {
    'Defaults': {
        'extract': False, 'noextract': True, 'threads': 1, 'kmers': 7
        }
}
fastqc_cmd = QiitaCommand(
    "FastQC", "Sequence QC", fastqc, req_params, opt_params,
    outputs, dflt_param_set)
plugin.register_command(fastqc_cmd)

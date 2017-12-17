# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .qc_trim import qc_trim

__all__ = ['qc_trim']

# Define the kneaddata command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'adapter': ['string', 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],  # 3' adapter
    'A': ['string', 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT'],  # 3' adapter for rev
    'quality-cutoff': ['integer', '15'],  # 3' quality cutoff
    'minimum-length': ['integer', '80'],  # min length after trimming
    'pair-filter': ['choice:["any", "both"]', 'any'],  # drop pairs whose mates
    # are filtered out
    'max-n': ['integer', '80'],  # maximum Ns to drop sequence
    'trim-n': ['boolean', 'True'],  # trim Ns on end of read
    'nextseq-trim': ['boolean', 'False'],  # NextSeq-specific quality trimming
    }
outputs = {'Adapter trimmed files': 'per_sample_FASTQ'}
dflt_param_set = {
    'KAPA HyperPlus with iTru': {
        'adapter': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        'A': 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT',
        'quality-cutoff': 15, 'minimum-length': 80, 'pair-filter': 'any',
        'max-n': 80, 'trim-n': True, 'nextseq-trim': False
        }
}

qc_trim_cmd = QiitaCommand(
    'Atropos v1.1.15', "Sequence QC - adapter trimming", qc_trim,
    req_params, opt_params, outputs, dflt_param_set)

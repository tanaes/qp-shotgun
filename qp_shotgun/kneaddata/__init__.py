# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaCommand

from .kneaddata import kneaddata

__all__ = ['kneaddata']

# Define the kneaddata command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # there are other parameters not included that will be ignored in this
    # configuration as we assume that the correct values were set in the env
    # by the admin installing the tools:
    # trimmomatic
    # bowtie2
    # --input # input FASTQ file (add a second argument for paired input files)
    # --output # directory to write output files
    # --output-prefix # prefix for output files [ DEFAULT : $SAMPLE_kneaddata ]
    # --log # filepath for log [ DEFAULT : $OUTPUT_DIR/$SAMPLE_kneaddata.log ]
    # --bowtie2 # path to bowtie executable
    # --bmtagger # path to bmtagger exectuable
    # --trf # path to TRF executable
    # --store-temp-output # store temp output files
    # --trimmomatic # the location of the trimmomatic (automatically
    #               # installed by pip install kneaddata)
    # --log-level
    'reference-db': ['choice:["default"]', 'default'],  # ref db
    'bypass-trim': ['boolean', 'False'],  # bypass the trim step
    'threads': ['integer', '1'],  # threads to run
    'processes': ['integer', '1'],  # processes to run
    'quality-scores': ['choice:["phred33", "phred64"]', 'phred33'],  # quality
    'run-bmtagger': ['boolean', 'False'],  # run BMTagger instead of Bowtie2
    'run-trf': ['boolean', 'False'],  # run TRF repeat finder tool
    'run-fastqc-start': ['boolean', 'True'],  # run FastQC on original data
    'run-fastqc-end': ['boolean', 'False'],  # run FastQC on filtered data
    # Trimmomatic options
    'max-memory': ['string', '500m'],  # max memory in mb [ DEFAULT : 500m ]
    # leaving as empty string to simply. note that in tests is not empty.
    'trimmomatic-options': ['string', ''],

    # Bowtie2 options
    # 'bowtie2-options': ['string', '--very-sensitive']
    # BMTagger options
    # # TRF options
    # 'match': ['integer', '2'], # matching weight
    # 'mismatch': ['integer', '7'], # mismatching penalty
    # 'delta': ['integer', '7'], # indel penalty
    # 'pm': ['integer', '80'], # match probability
    # 'pi': ['integer', '10'], # indel probability
    # 'minscore': ['integer', '50'], # mimimum alignment score to report
    # 'maxperiod': ['integer', '500m'] # maximum period size to report
    # FastQC options
    }
outputs = {
    'KneadData clean paired': 'per_sample_FASTQ',
    'KneadData clean unmatched R1': 'per_sample_FASTQ',
    'KneadData clean unmatched R2': 'per_sample_FASTQ',
    'KneadData clean R1': 'per_sample_FASTQ',
    }
dflt_param_set = {
    'Defaults': {
        'reference-db': 'default', 'bypass-trim': False, 'threads': 1,
        'processes': 1, 'quality-scores': 'phred33', 'run-bmtagger': False,
        'run-trf': False, 'run-fastqc-start': True, 'run-fastqc-end': False,
        'max-memory': '500m', 'trimmomatic-options': ''
        }
}

kneaddata_cmd = QiitaCommand(
    'KneadData 0.5.1', "Sequence QC - Human sequence removal", kneaddata,
    req_params, opt_params, outputs, dflt_param_set)

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


# Initialize the plugin
plugin = QiitaPlugin(
    'KneadData', '0.5.1', 'KneadData is a tool designed to perform quality '
    'control on metagenomic and metatranscriptomic sequencing data, '
    'especially data from microbiome experiments.')

# Define the HUMAnN2 command
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # there are other parameters not included that will be ignored in this
    # configuration as we assume that the correct values were set in the env
    # by the admin installing the tools:
    # trimmomatic
    # bowtie2

    # --input # input FASTQ file (add a second argument instance to run with paired input files)
    # --output # directory to write output files
    # --output-prefix # prefix for all output files [ DEFAULT : $SAMPLE_kneaddata ]
    # --log # filepath for log [ DEFAULT : $OUTPUT_DIR/$SAMPLE_kneaddata.log ]
    # --trimmomatic # path to trimmomatic executable
    # --bowtie2 # path to bowtie executable
    # --bmtagger # path to bmtagger exectuable
    # --trf # path to TRF executable
    'reference-db': ['choice:["human_genome"]', 'human_genome'], # ref db
    'bypass-trim': ['bool', 'False'], # bypass the trim step
    'threads': ['integer', '1'], # threads to run 
    'processes': ['integer', '1'], # processes to run
    'quality-scores': ['choice:["phred33","phred64"]', 'phred33'], # quality mapping
    'run-bmtagger': ['bool', 'False'], # run BMTagger instead of Bowtie2
    'run-trf': ['bool', 'False'], # run TRF repeat finder tool
    'run-fastqc-start': ['bool', 'True'], # run FastQC on original data
    'run-fastqc-end': ['bool', 'True'], # run FastQC on filtered data
    'store-temp-output': ['bool', 'False'], # store temp output files
    'log-level': ['choice:["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]',
                  'DEBUG'],

    # Trimmomatic options
    'max-memory': ['integer', '500'], # max memory in mb [ DEFAULT : 500 ]
    'trimmomatic-options': ['string', 'ILLUMINACLIP:$trimmomatic/adapters/'
                            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 '
                            'SLIDINGWINDOW:4:15 MINLEN:36'],

    # Bowtie2 options
    'bowtie2-options': ['string', '--very-sensitive']

    # BMTagger options
    
    # TRF options
    'match': ['integer', '2'], # matching weight
    'mismatch': ['integer', '7'], # mismatching penalty
    'delta': ['integer', '7'], # indel penalty
    'pm': ['integer', '80'], # match probability
    'pi': ['integer', '10'], # indel probability
    'minscore': ['integer', '50'], # mimimum alignment score to report
    'maxperiod': ['integer', '500'] # maximum period size to report

    # FastQC options
    }
outputs = {'per_sample_FASTQ': 'per_sample_FASTQ'}
dflt_param_set = {
    'Defaults': {
        'reference-db': 'human_genome', 'bypass-trim': False, 'threads': 1,
        'processes': 1, 'quality-scores': 'phred33', 'run-bmtagger': False,
        'run-trf': False, 'run-fastqc-start': True, 'run-fastqc-end': True,
        'store-temp-output': False, 'log-level': 'DEBUG', 'max-memory': 500,
        'trimmomatic-options': 'ILLUMINACLIP:$trimmomatic/adapters/'
                                'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 '
                                'SLIDINGWINDOW:4:15 MINLEN:36',
        'bowtie2-options': '--very-sensitive', 'match': 2, 'mismatch': 7,
        'delta': 7, 'pm': 80, 'pi': 10, 'minscore': 50, 'maxperiod': '500'}
}
kneaddata_cmd = QiitaCommand(
    "KneadData", "Sequence QC", kneaddata, req_params, opt_params,
    outputs, dflt_param_set)
plugin.register_command(kneaddata_cmd)
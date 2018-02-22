# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import os
from os.path import join
from tempfile import TemporaryDirectory
from qp_shotgun.utils import (
    _format_qc_params, make_read_pairs_per_sample,
    _run_commands, _per_sample_ainfo)
    
SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Taxonomic Level': 'levels', 'Number of threads': 'threads'}

def _read_fastq_seqs(filepath):
    # This reading method is adopted from qiime2 demux
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual.strip())

def generate_fna_file(fwd_seqs, rev_seqs, temp_path, map_file):
    # Combines reverse and forward seqs per sample
    # Returns filepaths of new combined files
    both_fp = []
    samples = make_read_pairs_per_sample(fwd_seqs, rev_seqs, map_file)
    for run_prefix, sample, f_fp, r_fp in samples:
        counter = 0
        output_fp = join(temp_path, sample, '_both.fna')
        output = open(output_fp, "a")
        # Loop through forward file
        for seq_header, seq, qual_header, qual in _read_fastq_seqs(f_fp):
            output.write("%s_%d\n" % (sample, counter))
            output.write("%s\n" % seq)
            count+=1
        # Loop through reverse file
        for seq_header, seq, qual_header, qual in _read_fastq_seqs(r_fp):
            output.write("%s_%d\n" % (sample, counter))
            output.write("%s\n" % seq)
            count+=1
        output.close()
        both_fp.append(output_fp)

    return both_fp

def shogun(qclient, job_id, parameters, out_dir):
    """Run Shogun with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    # Step 1 get the rest of the information need to run Atropos
    qclient.update_job_step(job_id, "Step 1 of 6: Collecting information")
    artifact_id = parameters['input']
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    qiime_map = prep_info['qiime-map']

    # Step 2 converting to fna
    print(fps, qiime_map)
    qclient.update_job_step(
        job_id, "Step 2 of 6: Converting to FNA for Shogun")
    temp_path = os.environ['QC_SHOGUN_TEMP_DP']
    with TemporaryDirectory(dir=temp_path, prefix='shogun_') as temp_dir:
        rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
        both_fps = generate_fna_file(
            fps['raw_forward_seqs'], rs, temp_dir, qiime_map)

    # Step 3 align
    qclient.update_job_step(job_id, "Step 3 of 6: Aligning FNA with Shogun")
    # shogun align \
    # --aligner utree \
    # --threads 2 \
    # --database ~/git_sw/oecophylla/test_data/test_dbs/shogun \
    # --input both.fna \
    # --output out

    # Step 4 taxonomic profile
    qclient.update_job_step(
        job_id, "Step 4 of 6: Taxonomic profile with Shogun")
    # shogun assign_taxonomy \
    # --aligner utree \
    # --database ~/git_sw/oecophylla/test_data/test_dbs/shogun \
    # --input out/alignment.utree.tsv \
    # --output out/profile.tsv

    # Step 5 functional profile
    qclient.update_job_step(
        job_id, "Step 5 of 6: Functional profile with Shogun")
    # shogun functional \
    # --database {params.db} \
    # --input alignment.utree.tsv \
    # --output {temp_dir} \
    # --level $level \

    # Step 6 functional profile
    qclient.update_job_step(job_id, "Step 6 of 6: Converting results to BIOM")
    # biom convert -i otu_table.taxonomy.txt -o otu_table.from_txt.biom
    # --table-type="OTU table" --process-obs-metadata taxonomy --to-hdf5

    ainfo = {}

    return True, ainfo, ""

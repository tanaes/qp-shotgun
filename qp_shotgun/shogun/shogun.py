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
from .utils import readfq
from qp_shotgun.utils import (make_read_pairs_per_sample)
import gzip
import tarfile

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Taxonomic Level': 'levels', 'Number of threads': 'threads'}


def generate_fna_file(fwd_seqs, rev_seqs, temp_path, map_file):
    # Combines reverse and forward seqs per sample
    # Returns filepaths of new combined files
    samples = make_read_pairs_per_sample(fwd_seqs, rev_seqs, map_file)
    output_fp = join(temp_path, '_combined.fna')
    output = open(output_fp, "a")
    count = 0
    for run_prefix, sample, f_fp, r_fp in samples:
        with gzip.open(f_fp, 'rt') as fp:
            # Loop through forward file
            for header, seq, qual in readfq(fp):
                output.write("%s_%d\n" % (sample, count))
                output.write("%s\n" % seq)
                count += 1
        with gzip.open(r_fp, 'rt') as fp:
            # Loop through reverse file
            for header, seq, qual in readfq(fp):
                output.write("%s_%d\n" % (sample, count))
                output.write("%s\n" % seq)
                count += 1
    output.close()

    return output_fp


def _format_params(parameters, func_params):
    params = {}
    # Loop through all of the commands alphabetically
    for param in func_params:
        # Find the value using long parameter names
        parameter = func_params[param]
        value = parameters[param]
        params[parameter] = value

    return params


def generate_shogun_align_commands(input_fp, out_dir, temp_dir, parameters):
    cmds = []
    cmds.append(
        'shogun align --aligner {aligner} --threads {threads} '
        '--database {database} --input {input} --output {output}'.format(
            aligner=parameters['aligner'],
            threads=parameters['threads'],
            database=parameters['database'],
            input=input_fp,
            output=out_dir))

    return cmds


def generate_shogun_assign_taxonomy_commands(input_fp, out_dir,
                                             temp_dir, parameters):
    cmds = []
    aln2ext = {'utree': 'tsv', 'burst': 'b6', 'bowtie2': 'sam'}
    ext = aln2ext[parameters['aligner']]
    cmds.append(
        'shogun assign_taxonomy --aligner {aligner} '
        '--database {database} --input {input} --output {output}'.format(
            aligner=parameters['aligner'],
            database=parameters['database'],
            input=join(out_dir, 'alignment.%s.%s' % (parameters['aligner'],
                                                     ext)),
            output=join(out_dir, 'profile.tsv')))

    return cmds


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
    qclient.update_job_step(
        job_id, "Step 2 of 6: Converting to FNA for Shogun")
    temp_path = os.environ['QC_SHOGUN_TEMP_DP']

    with TemporaryDirectory(dir=temp_path, prefix='shogun_') as temp_dir:
        rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
        comb_fp = generate_fna_file(
            fps['raw_forward_seqs'], rs, temp_dir, qiime_map)
        # Combining files
        parameters = _format_params(parameters, SHOGUN_PARAMS)
        with tarfile.open(parameters['database'], 'r:bz2') as db:
            db.extractall(temp_dir)
            parameters['database'] = join(temp_dir, db.getnames()[0])
            # Step 3 align
            qclient.update_job_step(
                job_id, "Step 3 of 6: Aligning FNA with Shogun")
            generate_shogun_align_commands(
                comb_fp, out_dir, temp_dir, parameters)
            # Step 4 taxonomic profile
            qclient.update_job_step(
                job_id, "Step 4 of 6: Taxonomic profile with Shogun")
            generate_shogun_assign_taxonomy_commands(
                comb_fp, out_dir, temp_dir, parameters)
            # Step 5 functional profile
            qclient.update_job_step(
                job_id, "Step 5 of 6: Functional profile with Shogun")
            # shogun functional \
            # --database {params.db} \
            # --input alignment.utree.tsv \
            # --output {temp_dir} \
            # --level $level \

            # Step 6 functional profile
            qclient.update_job_step(
                job_id, "Step 6 of 6: Converting results to BIOM")
            # biom convert -i otu_table.taxonomy.txt -o otu_table.from_txt.biom
            # --table-type="OTU table" --process-obs-metadata taxonomy
            # --to-hdf5

    ainfo = {}

    return True, ainfo, ""

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os.path import join
from tempfile import TemporaryDirectory
from .utils import readfq
from qp_shotgun.utils import (make_read_pairs_per_sample)
import gzip

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Taxonomic Level': 'levels', 'Number of threads': 'threads'}


def generate_fna_file(temp_path, samples):
    # Combines reverse and forward seqs per sample
    # Returns filepaths of new combined files
    output_fp = join(temp_path, 'combined.fna')
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


def generate_shogun_align_commands(input_fp, temp_dir, parameters):
    cmds = []
    cmds.append(
        'shogun align --aligner {aligner} --threads {threads} '
        '--database {database} --input {input} --output {output}'.format(
            aligner=parameters['aligner'],
            threads=parameters['threads'],
            database=parameters['database'],
            input=input_fp,
            output=temp_dir))

    return cmds


def generate_shogun_assign_taxonomy_commands(temp_dir, parameters):
    cmds = []
    aln2ext = {'utree': 'tsv', 'burst': 'b6', 'bowtie2': 'sam'}
    ext = aln2ext[parameters['aligner']]
    output_fp = join(temp_dir, 'profile.tsv')
    cmds.append(
        'shogun assign_taxonomy '
        '--aligner {aligner} '
        '--database {database} '
        '--input {input} --output {output}'.format(
            aligner=parameters['aligner'],
            database=parameters['database'],
            input=join(temp_dir, 'alignment.%s.%s' % (parameters['aligner'],
                                                      ext)),
            output=output_fp))

    return cmds, output_fp


def generate_shogun_functional_commands(profile_dir, temp_dir,
                                        parameters, sel_level):
    cmds = []
    cmds.append(
        'shogun functional '
        '--database {database} '
        '--input {input} '
        '--output {output} '
        '--level {level}'.format(
            database=parameters['database'],
            input=profile_dir,
            output=join(temp_dir, 'profile.functional.%s.tsv' % sel_level),
            level=sel_level))

    return cmds


def generate_biom_conversion_commands(input_fp, output_dir, level):
    cmds = []
    cmds.append(
        'biom convert -i {input} '
        '-o {output} '
        '--table-type="OTU table" '
        '--process-obs-metadata taxonomy --to-hdf5'.format(
            input=input_fp,
            output=join(output_dir, 'otu_table.%s.biom' % level)))

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

    with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:
        rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
        samples = make_read_pairs_per_sample(
            fps['raw_forward_seqs'], rs, qiime_map)

        comb_fp = generate_fna_file(temp_dir, samples)
        # Combining files
        parameters = _format_params(parameters, SHOGUN_PARAMS)
        # Step 3 align
        qclient.update_job_step(
            job_id, "Step 3 of 6: Aligning FNA with Shogun")
        generate_shogun_align_commands(
            comb_fp, temp_dir, parameters)
        # Step 4 taxonomic profile
        qclient.update_job_step(
            job_id, "Step 4 of 6: Taxonomic profile with Shogun")
        cmd, profile_fp = generate_shogun_assign_taxonomy_commands(
            temp_dir, parameters)
        # Step 5 functional profile
        qclient.update_job_step(
            job_id, "Step 5 of 6: Functional profile with Shogun")
        levels = ['genus', 'species', 'level']
        for level in levels:
            generate_shogun_functional_commands(
                profile_fp, temp_dir, parameters, level)
        # Step 6 functional profile
        qclient.update_job_step(
            job_id, "Step 6 of 6: Converting results to BIOM")
        for level in levels:
            input_fp = ('profile.functional.%s.tsv' % level)
            generate_biom_conversion_commands(input_fp, out_dir, level)

    ainfo = {}

    return True, ainfo, ""

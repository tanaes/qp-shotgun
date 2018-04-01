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
from qp_shotgun.utils import (make_read_pairs_per_sample, _run_commands)
import gzip
from qiita_client import ArtifactInfo

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
    output = join(temp_dir, 'profile.functional.%s.tsv' % sel_level)
    cmds.append(
        'shogun functional '
        '--database {database} '
        '--input {input} '
        '--output {output} '
        '--level {level}'.format(
            database=parameters['database'],
            input=profile_dir,
            output=output,
            level=sel_level))

    return cmds


def generate_shogun_redist_commands(profile_dir, temp_dir,
                                    parameters, sel_level):
    cmds = []
    output = join(temp_dir, 'profile.redist.%s.tsv' % sel_level)
    cmds.append(
        'shogun redistribute '
        '--database {database} '
        '--level {level} '
        '--input {input} '
        '--output {output}'.format(
            database=parameters['database'],
            input=profile_dir,
            output=output,
            level=sel_level))

    return cmds, output


def generate_biom_conversion_commands(input_fp, output_dir, level, version):
    cmds = []
    output = join(output_dir, 'otu_table.%s.%s.biom' % (level, version))
    cmds.append(
        'biom convert -i {input} '
        '-o {output} '
        '--table-type="OTU table" '
        '--process-obs-metadata taxonomy --to-hdf5'.format(
            input=input_fp,
            output=output))

    return cmds, output


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
    qclient.update_job_step(job_id, "Step 1 of 7: Collecting information")
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
        job_id, "Step 2 of 7: Converting to FNA for Shogun")

    with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:
        rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
        samples = make_read_pairs_per_sample(
            fps['raw_forward_seqs'], rs, qiime_map)

        # Combining files
        comb_fp = generate_fna_file(temp_dir, samples)

        # Formatting parameters
        parameters = _format_params(parameters, SHOGUN_PARAMS)

        # Step 3 align
        msg = "Step 3 of 7: Aligning FNA with Shogun"
        align_cmd = generate_shogun_align_commands(
            comb_fp, temp_dir, parameters)
        success, msg = _run_commands(
            qclient, job_id, align_cmd, msg, 'Shogun Align')
        if not success:
            return False, None, msg

        # Step 4 taxonomic profile
        msg = "Step 4 of 7: Taxonomic profile with Shogun"
        assign_cmd, profile_fp = generate_shogun_assign_taxonomy_commands(
            temp_dir, parameters)
        success, msg = _run_commands(
            qclient, job_id, assign_cmd, msg, 'Shogun taxonomy assignment')
        if not success:
            return False, None, msg

        # Step 5 redistribute profile
        msg = "Step 5 of 7: Redistributed profile with Shogun"
        levels = ['genus', 'species', 'level']
        redist_fps = []
        for level in levels:
            redist_cmd, output = generate_shogun_redist_commands(
                profile_fp, temp_dir, parameters, level)
            redist_fps.append(output)
            success, msg = _run_commands(
                qclient, job_id, redist_cmd, msg, 'Shogun redistribute')
            if not success:
                return False, None, msg

        # Step 6 functional profile
        msg = "Step 6 of 7: Functional profile with Shogun"
        levels = ['genus', 'species', 'level']
        func_fps = []
        for level in levels:
            func_cmd, output = generate_shogun_functional_commands(
                profile_fp, temp_dir, parameters, level)
            func_fps.append(output)
            success, msg = _run_commands(
                qclient, job_id, func_cmd, msg, 'Shogun functional')
            if not success:
                return False, None, msg

        # Step 6 functional profile
        qclient.update_job_step(
            job_id, "Step 7 of 7: Converting results to BIOM")
        func_biom_outputs = []
        redist_biom_outputs = []
        for func_fp, redist_fp in zip(func_fps, redist_fps):
            # Coverting funcitonal files to biom
            biom_cmd, output = generate_biom_conversion_commands(
                func_fp, out_dir, level, 'functional')
            success, msg = _run_commands(
                qclient, job_id, biom_cmd, msg, ' Functional Biom conversion')
            if not success:
                return False, None, msg
            else:
                func_biom_outputs.append(output)
            # Converting redistributed files to biom
            biom_cmd, output = generate_biom_conversion_commands(
                redist_fp, out_dir, level, 'redist')
            redist_biom_outputs.append(output)
            success, msg = _run_commands(
                qclient, job_id, biom_cmd, msg, 'Redistribute Biom conversion')
            if not success:
                return False, None, msg
            else:
                redist_biom_outputs.append(output)
    func_files_type_name = 'Functional Predictions'
    redist_files_type_name = 'Taxonomic Predictions'
    ainfo = [ArtifactInfo(func_files_type_name, 'BIOM', func_biom_outputs),
             ArtifactInfo(redist_files_type_name, 'BIOM', redist_biom_outputs)]

    return True, ainfo, ""

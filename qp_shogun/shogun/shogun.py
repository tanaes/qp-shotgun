# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os.path import join
from tempfile import TemporaryDirectory
from .utils import readfq, import_shogun_biom, shogun_db_functional_parser
from qp_shogun.utils import (make_read_pairs_per_sample, _run_commands)
import gzip
from qiita_client import ArtifactInfo
from biom import Table, util

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Number of threads': 'threads'}


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
                output.write(">%s_%d\n" % (sample, count))
                output.write("%s\n" % seq)
                count += 1
        with gzip.open(r_fp, 'rt') as fp:
            # Loop through reverse file
            for header, seq, qual in readfq(fp):
                output.write(">%s_%d\n" % (sample, count))
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
    output = join(temp_dir, 'functional')
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

    return cmds, output


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


def run_shogun_to_biom(in_fp, biom_in, out_dir, level, version):
    if version == 'redist':
        output_fp = join(out_dir, 'otu_table.%s.%s.biom'
                        % (level, version))
    else:
        output_fp = join(out_dir, 'otu_table.%s.%s.%s.biom'
                        % (level, version, biom_in[0]))
    tb = import_shogun_biom(in_fp, biom_in[1],
                            biom_in[2], biom_in[3])
    with util.biom_open(output_fp, 'w') as f:
        tb.to_hdf5(f, "shogun")

    return output_fp


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
        sys_msg = "Step 3 of 7: Aligning FNA with Shogun (%d/{0})"
        align_cmd = generate_shogun_align_commands(
            comb_fp, temp_dir, parameters)
        success, msg = _run_commands(
            qclient, job_id, align_cmd, sys_msg, 'Shogun Align')

        if not success:
            return False, None, msg

        # Step 4 taxonomic profile
        sys_msg = "Step 4 of 7: Taxonomic profile with Shogun (%d/{0})"
        assign_cmd, profile_fp = generate_shogun_assign_taxonomy_commands(
            temp_dir, parameters)
        success, msg = _run_commands(
            qclient, job_id, assign_cmd, sys_msg, 'Shogun taxonomy assignment')
        if not success:
            return False, None, msg

        # Step 5 redistribute profile
        sys_msg = "Step 5 of 7: Redistributed profile with Shogun (%d/{0})"
        levels = ['genus', 'species', 'strain']
        redist_fps = []
        for level in levels:
            redist_cmd, output = generate_shogun_redist_commands(
                profile_fp, temp_dir, parameters, level)
            redist_fps.append(output)
            success, msg = _run_commands(
                qclient, job_id, redist_cmd, sys_msg, 'Shogun redistribute')
            if not success:
                return False, None, msg

        # Step 6 functional profile
        sys_msg = "Step 6 of 7: Functional profile with Shogun (%d/{0})"
        levels = ['species']
        func_fp = ''
        for level in levels:
            func_cmd, output = generate_shogun_functional_commands(
                profile_fp, temp_dir, parameters, level)
            func_fp = output
            success, msg = _run_commands(
                qclient, job_id, func_cmd, sys_msg, 'Shogun functional')
            if not success:
                return False, None, msg
        # Step 6 functional profile
        sys_msg = "Step 7 of 7: Converting results to BIOM (%d/{0})"
        func_biom_outputs = []
        redist_biom_outputs = []
        # Converting redistributed files to biom
        redist_levels = ['genus', 'species', 'strain']
        for redist_fp, level in zip(redist_fps, redist_levels):
            biom_in = ["redist", None, '', True]
            output = run_shogun_to_biom(
                redist_fp, biom_in, out_dir, level, 'redist')
            redist_biom_outputs.append(output)
        # Coverting funcitonal files to biom
        func_db_fp = shogun_db_functional_parser(parameters['database'])
        for level in levels:

            func_to_biom_fps = [
                ["kegg.modules.coverage", func_db_fp['module'], 'module', False],
                ["kegg.modules", func_db_fp['module'], 'module', False],
                ["kegg.pathways.coverage", func_db_fp['pathway'],
                'pathway', False],
                ["kegg.pathways", func_db_fp['pathway'], 'pathway', False],
                ["kegg", func_db_fp['enzyme'], 'enzyme', True],
                ["normalized", func_db_fp['enzyme'], 'pathway', True]]

            for biom_in in func_to_biom_fps:
                biom_in_fp = join(func_fp, "profile.%s.%s.txt"
                                  % (level, biom_in[0]))
                output = run_shogun_to_biom(biom_in_fp, biom_in, out_dir,
                                            level, 'func')
                func_biom_outputs.append(output)

    func_files_type_name = 'Functional Predictions'
    redist_files_type_name = 'Taxonomic Predictions'
    ainfo = [ArtifactInfo(func_files_type_name, 'BIOM', func_biom_outputs),
             ArtifactInfo(redist_files_type_name, 'BIOM', redist_biom_outputs)]

    return True, ainfo, ""

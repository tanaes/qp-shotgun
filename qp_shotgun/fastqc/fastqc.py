# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import basename, join, splitext
from qp_shotgun.kneaddata.kneaddata import make_read_pairs_per_sample

from qiita_client import ArtifactInfo
from qiita_client.util import system_call


def format_fastqc_params(parameters):

    params = []

    for param in sorted(parameters):
        value = parameters[param]

        if str(value) == 'True':
            params.append('--%s' % param)
            continue
        elif str(value) == 'False':
            continue
        elif value:
            params.append('--%s %s' % (param, value))
            continue

    param_string = ' '.join(params)

    return(param_string)


def generate_fastqc_commands(forward_seqs, reverse_seqs, map_file, out_dir,
                             parameters):
    """Generates the FastQC commands

    Parameters
    ----------
    seqs : list of str
        The list of seqs filepaths
    map_file : str
        The path to the mapping file
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    list of str
        The FastQC commands

    Raises
    ------

    Notes
    -----
    """

    samples = make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file)

    cmds = []
    prefixes = []

    param_string = format_fastqc_params(parameters)
    fps = []
    for run_prefix, sample, f_fp, r_fp in samples:
        prefixes.append(run_prefix)
        if r_fp is None:
            cmds.append('mkdir -p %s; fastqc --outdir "%s" %s %s' %
                        (join(out_dir, run_prefix), join(out_dir, run_prefix),
                         param_string, f_fp))
            fps.append((f_fp, None))
        else:
            cmds.append('mkdir -p %s; fastqc --outdir "%s" %s %s %s' %
                        (join(out_dir, run_prefix), join(out_dir, run_prefix),
                         param_string, f_fp, r_fp))
            fps.append((f_fp, r_fp))

    return cmds, samples


def _guess_fastqc_filename(fp):
    f_p = basename(fp)

    exts = ['.fastq', '.fq', '.gz', '.gzip']
    while splitext(f_p)[1] in exts:
        f_p = splitext(f_p)[0]

    return "%s_fastqc.html" % f_p, "%s_fastqc.zip" % f_p


def _per_sample_ainfo(out_dir, samples):
    ainfo = []
    for run_prefix, sample, f_fp, r_fp in samples:
        sam_out_dir = join(out_dir, run_prefix)

        ainfo += [
            ArtifactInfo('FastQC html summary', 'html_summary',
                         [(join(sam_out_dir, _guess_fastqc_filename(f_fp)[0]),
                          'html_summary')]),
            ArtifactInfo('FastQC data summary', 'zip_file',
                         [(join(sam_out_dir, _guess_fastqc_filename(f_fp)[1]),
                          'zip_file')])]
        if r_fp:
            ainfo += [
                ArtifactInfo('FastQC html summary', 'html_summary',
                             [(join(sam_out_dir,
                               _guess_fastqc_filename(r_fp)[0]),
                               'html_summary')]),
                ArtifactInfo('FastQC data summary', 'zip_file',
                             [(join(sam_out_dir,
                               _guess_fastqc_filename(r_fp)[1]),
                               'zip_file')])]
    print(ainfo)
    return ainfo


def _run_commands(qclient, job_id, commands, msg):
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % i)
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = ("Error running FastQC:\nStd out: %s\nStd err: %s"
                         "\n\nCommand run was:\n%s"
                         % (std_out, std_err, cmd))
            return False, error_msg

    return True, ""


def fastqc(qclient, job_id, parameters, out_dir):
    """Run FastQC with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        Yhe path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    # Step 1 get the rest of the information need to run FastQC
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['input']

    # removing input from parameters so it's not part of the final command
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    qiime_map = prep_info['qiime-map']

    # Step 2 generating command FastQC
    qclient.update_job_step(job_id, "Step 2 of 3: Generating"
                                    " FastQC command")
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    commands, samples = generate_fastqc_commands(fps['raw_forward_seqs'],
                                                 rs, qiime_map, out_dir,
                                                 parameters)

    # Step 3 execute FastQC
    msg = "Step 3 of 3: Executing FastQC job (%d/{0})".format(len(commands))
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Step 4 generating artifacts
    ainfo = _per_sample_ainfo(out_dir, samples)

    return True, ainfo, ""

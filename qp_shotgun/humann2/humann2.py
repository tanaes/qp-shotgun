# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir
from os.path import basename, join, exists

from future.utils import viewitems
from functools import partial

from qiita_client import ArtifactInfo
from qiita_client.util import system_call, get_sample_names_by_run_prefix


def generate_humann2_analysis_commands(forward_seqs, reverse_seqs, map_file,
                                       out_dir, parameters, assure_od=True):
    """Generates the HUMAnN2 commands

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name
    assure_od : boolean, optional
        Assure that the output dirs for HUMAnN2 exist. Useful for testing.

    Returns
    -------
    list of str
        The HUMAnN2 commands

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    The forward reads filename has to be a perfect match with the value stored
    in run_prefix. This is not a requirement for the reverse reads. However,
    note that we assume that the filenames of the forward and reverse reads
    are pretty similar so if sorted they will have matching orders.
    """
    # making sure the forward and reverse reads are in the same order
    forward_seqs.sort()
    if reverse_seqs:
        if len(forward_seqs) != len(reverse_seqs):
            raise ValueError('Your reverse and forward files are of different '
                             'length. Forward: %s. Reverse: %s.' %
                             (', '.join(forward_seqs),
                              ', '.join(reverse_seqs)))
        reverse_seqs.sort()

    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # we match sample name and forward filename
    samples = []
    for i, fname in enumerate(forward_seqs):
        f = basename(fname)
        # removing extentions: fastq or fastq.gz
        if 'fastq' in f.lower().rsplit('.', 2):
            f = f[:f.lower().rindex('.fastq')]
        # this try/except block is simply to retrieve all possible errors
        # and display them in the next if block
        try:
            samples.append((fname, f, sn_by_rp[f]))
            if reverse_seqs:
                fr = basename(reverse_seqs[i])
                if 'fastq' in fr.lower().rsplit('.', 2):
                    fr = fr[:fr.lower().rindex('.fastq')]
                samples.append((reverse_seqs[i], fr, sn_by_rp[f]))
            del sn_by_rp[f]
        except KeyError:
            pass

    if sn_by_rp:
        raise ValueError(
            'Some run_prefix values do not match your sample names: %s'
            % ', '.join(sn_by_rp.keys()))

    cmds = []
    params = ['--%s "%s"' % (k, v) for k, v in viewitems(parameters) if v]
    for ffn, fn, s in samples:
        od = join(out_dir, fn)
        # just making sure the output directory exists
        if assure_od and not exists(od):
            mkdir(od)
        cmds.append('humann2 --input "%s" --output "%s" --output-basename '
                    '"%s" --output-format biom %s' % (ffn, od, s,
                                                      ' '.join(params)))

    return cmds


def _run_commands(qclient, job_id, commands, msg):
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % i)
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            print '-->', cmd
            error_msg = ("Error running HUMANn2:\nStd out: %s\nStd err: %s"
                         % (std_out, std_err))
            return False, error_msg

    return True, ""


def humann2(qclient, job_id, parameters, out_dir):
    """Run humann2 with the given parameters

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run HUMAnN2
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job
    """
    # Step 1 get the rest of the information need to run humann2
    qclient.update_job_step(job_id, "Step 1 of 5: Collecting information")
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

    # Step 2 generating command humann2
    qclient.update_job_step(job_id, "Step 2 of 5: Generating HUMANn2 command")
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    commands = generate_humann2_analysis_commands(fps['raw_forward_seqs'], rs,
                                                  qiime_map, out_dir,
                                                  parameters)

    # Step 3 execute humann2
    msg = "Step 3 of 5: Executing HUMANn2 job (%d/{0})".format(len(commands))
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Step 4 merge tables
    commands = []
    commands.append(('humann2_join_tables -i {0} -o {0}/genefamilies.biom '
                     '--file_name genefamilies').format(out_dir))
    commands.append(('humann2_join_tables -i {0} -o {0}/pathcoverage.biom '
                     '--file_name pathcoverage').format(out_dir))
    commands.append(('humann2_join_tables -i {0} -o {0}/pathabundance.biom '
                    '--file_name pathabundance_relab').format(out_dir))
    msg = "Step 4 of 5: Merging resulting tables job (%d/3)"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Step 5 generating re-normalized tables
    commands = []
    commands.append(('humann2_renorm_table -i {0}/genefamilies.biom -u cpm '
                     '-o {0}/genefamilies_cpm.biom').format(out_dir))
    commands.append(('humann2_renorm_table -i {0}/pathcoverage.biom -u relab '
                     '-o {0}/pathcoverage_relab.biom').format(out_dir))
    commands.append(('humann2_renorm_table -i {0}/pathabundance.biom -u relab '
                     '-o {0}/pathabundance_relab.biom').format(out_dir))
    msg = "Step 5 of 5: Re-normalizing tables (%d/3)"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Generating artifact
    pb = partial(join, out_dir)
    filepaths = [(pb('genefamilies.biom'), 'biom'),
                 (pb('pathcoverage.biom'), 'biom'),
                 (pb('pathabundance.biom'), 'biom'),
                 (pb('genefamilies_cpm.biom'), 'biom'),
                 (pb('pathcoverage_relab.biom'), 'biom'),
                 (pb('pathabundance_relab.biom'), 'biom')]
    ainfo = [ArtifactInfo('OTU table', 'BIOM', filepaths)]

    return True, ainfo, ""

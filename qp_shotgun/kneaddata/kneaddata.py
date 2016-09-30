# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import basename, join

from future.utils import viewitems

from qiita_client.util import get_sample_names_by_run_prefix


def make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file):
    """Recovers read pairing information

    Parameters
    ----------
    forward_seqs : list of str
        The list of forward seqs filepaths
    reverse_seqs : list of str
        The list of reverse seqs filepaths
    map_file : str
        The path to the mapping file

    Returns
    -------
    samples: list of tup
        list of 3-tuples with run prefix, sample name, fwd read fp, rev read fp

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    """

    # sort forward seqs
    forward_seqs.sort()

    # check that rev seqs are same len
    if reverse_seqs:
        if len(forward_seqs) != len(reverse_seqs):
            raise ValueError('Your reverse and forward files are of different '
                             'length. Forward: %s. Reverse: %s.' %
                             (', '.join(forward_seqs),
                              ', '.join(reverse_seqs)))
        reverse_seqs.sort()

    # get run prefixes
    # These are prefixes that should match uniquely to forward reads
    # sn_by_rp is dict of samples keyed by run prefixes
    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # make pairings
    samples = []
    used_prefixes = []
    for i, f_fp in enumerate(forward_seqs):
        # f_fp is the fwd read filepath
        f_fn = basename(f_fp)

        # iterate over run prefixes and make sure only one matches
        run_prefix = None
        for rp in sn_by_rp:
            if f_fn.startswith(rp) and run_prefix is None:
                run_prefix = rp
            elif f_fn.startswith(rp) and run_prefix is not None:
                raise ValueError('Multiple run prefixes match this fwd read: '
                                 '%s' % f_fn)

        # make sure that we got one matching run prefix:
        if run_prefix is None:
            raise ValueError('No run prefix matching this fwd read: %s'
                             % f_fn)

        if run_prefix in used_prefixes:
            raise ValueError('This run prefix matches multiple fwd reads: '
                             '%s' % run_prefix)

        # if we have reverse reads, make sure the matching pair also
        # matches the run prefix:
        if (reverse_seqs and not
                basename(reverse_seqs[i]).startswith(run_prefix)):
            raise ValueError('Reverse read does not match this run prefix. '
                             'Run prefix: %s\nForward read: %s\n'
                             'Reverse read: %s\n' %
                             (run_prefix, f_fn, basename(reverse_seqs[i])))

        used_prefixes.append(run_prefix)
        # create the tuple for this read set
        if reverse_seqs:
            samples.append((run_prefix, sn_by_rp[run_prefix], f_fp,
                            reverse_seqs[i]))
        else:
            samples.append((run_prefix, sn_by_rp[run_prefix], f_fp, None))

    return(samples)


def format_kneaddata_params(parameters):

    params = []

    for param in sorted(parameters):
        value = parameters[param]

        if value == True
            params.append('--%s' % param)
            continue
        elif value == False
            continue
        elif value
            params.append('--%s %s' % (param, value))
            continue

    param_string = ' '.join(params)

    return(param_string)


def generate_kneaddata_commands(forward_seqs, reverse_seqs, map_file,
                                out_dir, parameters):
    """Generates the KneadData commands

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

    Returns
    -------
    list of str
        The KneadData commands

    Raises
    ------
    ValueError
        If the rev is not an empty list and the same length than fwd seqs
        The prefixes of the run_prefix don't match the file names

    Notes
    -----
    """
    # making sure the forward and reverse reads are in the same order

    # we match filenames, samples, and run prefixes
    samples = make_read_pairs_per_sample(forward_seqs, reverse_seqs, map_file)

    cmds = []

    param_string = format_kneaddata_params(parameters)
    for run_prefix, sample, f_fp, r_fp in samples:
        if r_fp is None:
            cmds.append('kneaddata --input "%s" --output "%s" '
                        '--output-prefix "%s" %s' % 
                        (f_fp, join(out_dir, run_prefix),
                         run_prefix, param_string))
        else:
            cmds.append('kneaddata --input "%s" --input "%s" --output "%s" '
                        '--output-prefix "%s" %s' % 
                        (f_fp, join(out_dir, run_prefix),
                         run_prefix, param_string))

    return cmds


def kneaddata(qclient, job_id, parameters, out_dir):
    """Run kneaddata with the given parameters

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
    # Step 1 get the rest of the information need to run kneaddata
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['input_data']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    qiime_map = prep_info['qiime-map']

    # Step 2 generating command kneaddata
    qclient.update_job_step(job_id, "Step 2 of 3: "
                            "Generating kneaddata command")
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    commands = generate_kneaddata_commands(fps['raw_forward_seqs'], rs,
                                           qiime_map, out_dir,
                                           parameters)
    # Step 3 execute kneaddata: TODO
    qclient.update_job_step(job_id, "Step 3 of 3: Executing kneaddata")

    commands_len = len(commands)

    # artifacts_info = []

    return True, artifacts_info, ""

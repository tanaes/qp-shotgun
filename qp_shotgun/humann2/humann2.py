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

def make_read_sets_per_sample(files, map_file):
    """Recovers read set information from kneaddata output

    Parameters
    ----------
    files : list of str
        The list of sequence filepaths from the kneaddata artifact
    map_file : str
        The path to the mapping file

    Returns
    -------
    read_sets: list of tup
        list of 7-tuples with run prefix, sample name, fwd paired read fp,
        rev paired read fp, fwd unpaired read fp, rev unpaired read fp, and
        single fwd read fp.

    Raises
    ------
    ValueError
        If there are files matching the kneaddata paired file naming convention
        (_paired_1.fastq.gz, _unmatched_1.fastq.gz) but which don't have all 4
        outputs.
    ValueError
        If there are files matching the kneaddata paired file naming convention
        (_paired_1.fastq.gz, _unmatched_1.fastq.gz) in addition to fastq.gz
        files that do not match naming convention (the latter are interpreted
        as single read files).
    ValueError
        If there are no *.fastq.gz files in the artifact

    Notes
    -----
    """

    # sort through the filenames and bin into sequence type lists
    fwd_paired = []
    fwd_unpaired = []
    rev_paired = []
    rev_unpaired = []
    single = []

    for fp in files:
        if str.endswith('_paired_1.fastq.gz'):
            fwd_paired.append(fp)
        elif str.endswith('_paired_2.fastq.gz'):
            rev_paired.append(fp)
        elif str.endswith('_unmatched_1.fastq.gz'):
            fwd_unpaired.append(fp)
        elif str.endswith('_unmatched_2.fastq.gz'):
            rev_unpaired.append(fp)
        elif str.endswith('.fastq.gz'):
            single.append(fp)

    # check that seq lists are same len
    if not (len(fwd_paired) == len(fwd_unpaired) == 
            len(rev_paired) == len(rev_unpaired)):
        raise ValueError('There are not equal numbers of forward paired, '
                         'forward unpaired, reverse paired, and reverse '
                         'unpaired sequences.')

    # check that there aren't both paired and single sequences
    if len(single) > 0 and len(fwd_paired) > 0:
        raise ValueError('There are both paired-end and single-end sequences.')

    # fill out unused seq file types with None and check that there exist files
    if len(fwd_paired) > 0:
        single = [None] * len(fwd_paired)
    elif len(single) > 0:
        fwd_paired = [None] * len(single)
        fwd_unpaired = [None] * len(single)
        rev_paired = [None] * len(single)
        rev_unpaired = [None] * len(single)
    else:
        raise ValueError('There are no *.fastq.gz files in the artifact')
    
    # make the 5-tuple of sequence filepaths
    seq_files = izip_longest(fwd_paired.sort(), rev_paired.sort(),
                             fwd_unpaired.sort(), rev_unpaired.sort(),
                             single.sort())

    # get run prefixes
    # These are prefixes that should match uniquely to forward reads
    # sn_by_rp is dict of samples keyed by run prefixes
    sn_by_rp = get_sample_names_by_run_prefix(map_file)

    # make sets
    read_sets = []
    used_prefixes = set()

    for i, f_p, r_p, f_u, r_u, s in enumerate(seq_files):
        # pick file basename
        if f_p is None:
            fn = basename(s)
        else:
            fn = basename(f_p)

        # iterate over run prefixes and make sure only one matches
        run_prefix = None
        for rp in sn_by_rp:
            if fn.startswith(rp) and run_prefix is None:
                run_prefix = rp
            elif fn.startswith(rp) and run_prefix is not None:
                raise ValueError('Multiple run prefixes match this read file: '
                                 '%s' % fn)

        # make sure that we got one matching run prefix:
        if run_prefix is None:
            raise ValueError('No run prefix matching this read file: %s'
                             % fn)

        if run_prefix in used_prefixes:
            raise ValueError('This run prefix matches multiple read files: '
                             '%s' % run_prefix)

        # if paired, check that all files match run prefix
        if s is None:
            if not (r_p.startsiwth(run_prefix) and
                    f_u.startswith(run_prefix) and
                    r_u.startswith(run_prefix)):
                raise ValueError('Not all read files match run prefix.'
                                 '\nRun prefix: %s\nForward paired: %s\n'
                                 'Reverse paired: %s\nForward unpaired: %s\n'
                                 'Reverse unpaired: %s\n' %
                                 (run_prefix, f_p, r_p, f_u, r_u))

        read_sets.append((run_prefix, sn_by_rp[run_prefix], f_p, r_p,
                          f_u, r_u, s ))

        used_prefixes.add(run_prefix)

    return(read_sets)

def make_single_fastq_gz(read_sets, out_dir, params):


def generate_humann2_analysis_commands(forward_seqs, reverse_seqs, map_file,
                                       out_dir, parameters):
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
    params = []
    for k, v in viewitems(parameters):
        if v is False or v in ['False', 'default', '']:
            continue
        if v is True or v == 'True':
            params.append('--%s' % k)
        else:
            params.append('--%s "%s"' % (k, v))
    for ffn, fn, s in samples:
        od = join(out_dir, fn)
        # just making sure the output directory exists
        if not exists(od):
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
    qclient.update_job_step(job_id, "Step 1 of 6: Collecting information")
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
    qclient.update_job_step(job_id, "Step 2 of 6: Generating HUMANn2 command")
    rs = fps['raw_reverse_seqs'] if 'raw_reverse_seqs' in fps else []
    commands = generate_humann2_analysis_commands(fps['raw_forward_seqs'], rs,
                                                  qiime_map, out_dir,
                                                  parameters)

    # Step 3 execute humann2
    msg = "Step 3 of 6: Executing HUMANn2 job (%d/{0})".format(len(commands))
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Step 4 merge tables
    commands = []
    commands.append(('humann2_join_tables -i {0} -o {0}/genefamilies.biom '
                     '--file_name genefamilies --search-subdirectories '
                     '--verbose').format(out_dir))
    commands.append(('humann2_join_tables -i {0} -o {0}/pathcoverage.biom '
                     '--file_name pathcoverage --search-subdirectories '
                     '--verbose').format(out_dir))
    commands.append(('humann2_join_tables -i {0} -o {0}/pathabundance.biom '
                     '--file_name pathabundance --search-subdirectories '
                     '--verbose').format(out_dir))
    msg = "Step 4 of 6: Merging resulting tables job (%d/3)"
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
    msg = "Step 5 of 6: Re-normalizing tables (%d/3)"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Step 6 stratifiying re-normalized tables
    commands = []
    pb = partial(join, out_dir)
    cmd = "humann2_split_stratified_table --input %s --output %s"
    commands.append(cmd % (pb(out_dir, 'genefamilies_cpm.biom'), out_dir))
    commands.append(cmd % (pb(out_dir, 'pathcoverage_relab.biom'), out_dir))
    commands.append(cmd % (pb(out_dir, 'pathabundance_relab.biom'), out_dir))
    msg = "Step 6 of 6: Stratifiying re-normalizing tables (%d/3)"
    success, msg = _run_commands(qclient, job_id, commands, msg)
    if not success:
        return False, None, msg

    # Generating 6 artifacts, separation is important for analysis
    ainfo = [
        ArtifactInfo('Gene family table', 'BIOM',
                     [(pb('genefamilies.biom'), 'biom')]),
        ArtifactInfo('Path coverage table', 'BIOM',
                     [(pb('pathcoverage.biom'), 'biom')]),
        ArtifactInfo('Path abundance table', 'BIOM',
                     [(pb('pathabundance.biom'), 'biom')]),
        ArtifactInfo('Gene family CMP table', 'BIOM',
                     [(pb('genefamilies_cpm.biom'), 'biom')]),
        ArtifactInfo('Path coverage RELAB table', 'BIOM',
                     [(pb('pathcoverage_relab.biom'), 'biom')]),
        ArtifactInfo('Path abundance RELAB table', 'BIOM',
                     [(pb('pathabundance_relab.biom'), 'biom')]),
        ArtifactInfo('Gene family CMP table - stratified', 'BIOM',
                     [(pb('genefamilies_cpm_stratified.biom'), 'biom')]),
        ArtifactInfo('Path coverage RELAB table - stratified', 'BIOM',
                     [(pb('pathcoverage_relab_stratified.biom'), 'biom')]),
        ArtifactInfo('Path abundance RELAB table - stratified', 'BIOM',
                     [(pb('pathabundance_relab_stratified.biom'), 'biom')]),
        ArtifactInfo('Gene family CMP table - unstratified', 'BIOM',
                     [(pb('genefamilies_cpm_unstratified.biom'), 'biom')]),
        ArtifactInfo('Path coverage RELAB table - unstratified', 'BIOM',
                     [(pb('pathcoverage_relab_unstratified.biom'), 'biom')]),
        ArtifactInfo('Path abundance RELAB table - unstratified', 'BIOM',
                     [(pb('pathabundance_relab_unstratified.biom'), 'biom')])]

    return True, ainfo, ""

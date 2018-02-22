# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Taxonomic Level': 'levels', 'Number of threads': 'threads'}


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
    # for sample in sample_S22282 sample_S22205
    # do
    #
    #     seqtk mergepe ${sample}.R1.*.gz ${sample}.R2.*.gz | \
    #     seqtk seq -A > ${sample}.fna
    #
    # done
    # sample = sys.argv[1]
    # with open('%s.out.fna' % sample, 'w') as o:
    #     with open('%s.fna' % sample) as f:
    #         seq = 1
    #         for line in f:
    #             if line.startswith('>'):
    #                 lineout = '>%s_%s\n' * (sample, seq)
    #                 seq += 1
    #                 o.write(lineout)
    #             else:
    #                 o.write(line)

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

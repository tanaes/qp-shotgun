# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from qiita_client.testing import PluginTestCase
import os
from os import remove
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import TemporaryDirectory
from qp_shotgun import plugin
from tempfile import mkdtemp
from json import dumps
from functools import partial
from qp_shotgun.shogun.utils import (
    get_dbs, get_dbs_list, generate_shogun_dflt_params)
from qp_shotgun.shogun.shogun import (
    generate_shogun_align_commands, _format_params,
    generate_shogun_assign_taxonomy_commands, generate_fna_file,
    generate_shogun_functional_commands, generate_shogun_redist_commands,
    generate_biom_conversion_commands, shogun)

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Number of threads': 'threads'}


class ShogunTests(PluginTestCase):

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self.db_path = os.environ["QC_SHOGUN_DB_DP"]
        self.params = {
            'Database': join(self.db_path, 'shogun'),
            'Aligner tool': 'bowtie2',
            'Number of threads': 1
        }
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_get_dbs(self):
        db_path = self.db_path
        obs = get_dbs(db_path)
        exp = {'shogun': join(db_path, 'shogun')}

        self.assertEqual(obs, exp)

    def test_get_dbs_list(self):
        db_path = self.db_path
        obs = get_dbs_list(db_path)
        exp = join(join('"'+db_path, 'shogun')+'"')

        self.assertEqual(obs, exp)

    def test_generate_shogun_dflt_params(self):
        obs = generate_shogun_dflt_params()
        exp = {
            'shogun_bowtie2': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'bowtie2',
                'Number of threads': 1},
            'shogun_utree': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'utree',
                'Number of threads': 1},
            'shogun_burst': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'burst',
                'Number of threads': 1}}

        self.assertEqual(obs, exp)

    def test_generate_fna_file(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as fp:
            sample = [
                ('s1', 'SKB8.640193', 'support_files/kd_test_1_R1.fastq.gz',
                 'support_files/kd_test_1_R2.fastq.gz')
                ]
            exp = join(fp, 'combined.fna')
            obs = generate_fna_file(fp, sample)

        self.assertEqual(obs, exp)

    def test_format_shogun_params(self):
        obs = _format_params(self.params, SHOGUN_PARAMS)
        exp = {
            'database': join(self.db_path, 'shogun'),
            'aligner': 'bowtie2',
            'threads': 1
        }

        self.assertEqual(obs, exp)

    def test_generate_shogun_align_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:

            exp_cmd = [
                ('shogun align --aligner bowtie2 --threads 1 '
                 '--database %sshogun --input %s/combined.fna '
                 '--output %s') %
                (self.db_path, temp_dir, temp_dir)
                ]

            params = _format_params(self.params, SHOGUN_PARAMS)
            obs_cmd = generate_shogun_align_commands(
                join(temp_dir, 'combined.fna'), temp_dir, params)

        self.assertEqual(obs_cmd, exp_cmd)

    def test_generate_shogun_assign_taxonomy_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:

            exp_cmd = [
                ('shogun assign_taxonomy --aligner bowtie2 '
                 '--database %sshogun --input %s/alignment.bowtie2.sam '
                 '--output %s/profile.tsv') %
                (self.db_path, temp_dir, temp_dir)
                ]
            exp_output_fp = join(temp_dir, 'profile.tsv')
            params = _format_params(self.params, SHOGUN_PARAMS)
            obs_cmd, obs_output_fp = generate_shogun_assign_taxonomy_commands(
                temp_dir, params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_output_fp, exp_output_fp)

    def test_generate_shogun_functional_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:

            exp_cmd = [
                ('shogun functional '
                 '--database %sshogun --input %s '
                 '--output %s --level species') %
                (self.db_path, join(temp_dir, 'profile.tsv'),
                 join(temp_dir, 'functional'))
                ]
            profile_dir = join(temp_dir, 'profile.tsv')
            params = _format_params(self.params, SHOGUN_PARAMS)
            obs_cmd, output = generate_shogun_functional_commands(
                profile_dir, temp_dir, params, 'species')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_generate_shogun_redist_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:

            exp_cmd = [
                ('shogun redistribute '
                 '--database %sshogun --level species --input %s '
                 '--output %s') %
                (self.db_path, join(temp_dir, 'profile.tsv'),
                 join(temp_dir, 'profile.redist.species.tsv'))
                ]
            profile_dir = join(temp_dir, 'profile.tsv')
            params = _format_params(self.params, SHOGUN_PARAMS)
            obs_cmd, output = generate_shogun_redist_commands(
                profile_dir, temp_dir, params, 'species')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_generate_biom_conversion_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:
            exp_cmd = [
                ('biom convert -i %s '
                 '-o %s '
                 '--table-type="OTU table" '
                 '--to-hdf5') %
                (join(temp_dir, 'profile.tsv'),
                 join(out_dir, 'otu_table.species.redist.biom'))
                ]
            profile_fp = join(temp_dir, 'profile.tsv')
            obs_cmd, output = generate_biom_conversion_commands(
                profile_fp, out_dir, 'species', 'redist')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_shogun_bt2(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
        fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')

        copyfile('support_files/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile('support_files/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
        copyfile('support_files/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        copyfile('support_files/S22282_S102_L001_R2_001.fastq.gz', fp2_2)

        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test Shogun artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shotgun', '0.0.1', 'Shogun']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 2 artifacts in total
        self.assertEqual(2, len(ainfo))

        obs_func_fps = []
        obs_redist_fps = []
        ainfo_func = ainfo[0]
        ainfo_redist = ainfo[1]
        self.assertEqual('BIOM', ainfo_func.artifact_type)
        self.assertEqual('BIOM', ainfo_redist.artifact_type)
        obs_func_fps = ainfo_func.files
        obs_redist_fps = ainfo_redist.files

        od = partial(join, out_dir)
        func_prefix = "species"
        exp_func_fps = [
            od("otu_table.%s.kegg.modules.coverage.biom" % func_prefix),
            od("otu_table.%s.kegg.modules.biom" % func_prefix),
            od("otu_table.%s.kegg.pathways.coverage.biom" % func_prefix),
            od("otu_table.%s.kegg.pathways.biom" % func_prefix),
            od("otu_table.%s.kegg.biom" % func_prefix),
            od("otu_table.%s.normalized.biom" % func_prefix)]

        exp_redist_fps = [
            od('otu_table.genus.redist.biom'),
            od('otu_table.species.redist.biom'),
            od('otu_table.strain.redist.biom')]

        self.assertEqual(obs_func_fps, exp_func_fps)
        self.assertEqual(obs_redist_fps, exp_redist_fps)

#    def test_shogun_burst(self):
#        # generating filepaths
#        in_dir = mkdtemp()
#        self._clean_up_files.append(in_dir)

#        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
#        fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
#        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
#        fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')

#        copyfile('support_files/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
#        copyfile('support_files/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
#        copyfile('support_files/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
#        copyfile('support_files/S22282_S102_L001_R2_001.fastq.gz', fp2_2)

#        # inserting new prep template
#        prep_info_dict = {
#            'SKB8.640193': {'run_prefix': 'S22205_S104'},
#            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
#        data = {'prep_info': dumps(prep_info_dict),
#                # magic #1 = testing study
#                'study': 1,
#                'data_type': 'Metagenomic'}
#        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

#        # inserting artifacts
#        data = {
#            'filepaths': dumps([
#                (fp1_2, 'raw_reverse_seqs'),
#                (fp2_1, 'raw_forward_seqs'),
#                (fp2_2, 'raw_reverse_seqs')]),
#            'type': "per_sample_FASTQ",
#            'name': "Test Shogun artifact",
#            'prep': pid}
#        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

#        self.params['input'] = aid
#        self.params['Aligner tool'] = 'burst'
#        data = {'user': 'demo@microbio.me',
#                'command': dumps(['qp-shotgun', '0.0.1', 'Shogun']),
#                'status': 'running',
#                'parameters': dumps(self.params)}
#        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

#        out_dir = mkdtemp()
#        self._clean_up_files.append(out_dir)

#        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

#        self.assertEqual("", msg)
#        self.assertTrue(success)

#        # we are expecting 2 artifacts in total
#        self.assertEqual(2, len(ainfo))

#        obs_func_fps = []
#        obs_redist_fps = []
#        ainfo_func = ainfo[0]
#        ainfo_redist = ainfo[1]
#        self.assertEqual('BIOM', ainfo_func.artifact_type)
#        self.assertEqual('BIOM', ainfo_redist.artifact_type)
#        obs_func_fps = ainfo_func.files
#        obs_redist_fps = ainfo_redist.files

#        od = partial(join, out_dir)
#        func_prefix = "species"
#        exp_func_fps = [
#            od("otu_table.%s.kegg.modules.coverage.biom" % func_prefix),
#            od("otu_table.%s.kegg.modules.biom" % func_prefix),
#            od("otu_table.%s.kegg.pathways.coverage.biom" % func_prefix),
#            od("otu_table.%s.kegg.pathways.biom" % func_prefix),
#            od("otu_table.%s.kegg.biom" % func_prefix),
#            od("otu_table.%s.normalized.biom" % func_prefix)]

#        exp_redist_fps = [
#            od('otu_table.genus.redist.biom'),
#            od('otu_table.species.redist.biom'),
#            od('otu_table.strain.redist.biom')]

#        self.assertEqual(obs_func_fps, exp_func_fps)
#        self.assertEqual(obs_redist_fps, exp_redist_fps)

    def test_shogun_utree(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
        fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')

        copyfile('support_files/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile('support_files/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
        copyfile('support_files/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        copyfile('support_files/S22282_S102_L001_R2_001.fastq.gz', fp2_2)

        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test Shogun artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        self.params['Aligner tool'] = 'utree'
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shotgun', '0.0.1', 'Shogun']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 2 artifacts in total
        self.assertEqual(2, len(ainfo))

        obs_func_fps = []
        obs_redist_fps = []
        ainfo_func = ainfo[0]
        ainfo_redist = ainfo[1]
        self.assertEqual('BIOM', ainfo_func.artifact_type)
        self.assertEqual('BIOM', ainfo_redist.artifact_type)
        obs_func_fps = ainfo_func.files
        obs_redist_fps = ainfo_redist.files

        od = partial(join, out_dir)
        func_prefix = "species"
        exp_func_fps = [
            od("otu_table.%s.kegg.modules.coverage.biom" % func_prefix),
            od("otu_table.%s.kegg.modules.biom" % func_prefix),
            od("otu_table.%s.kegg.pathways.coverage.biom" % func_prefix),
            od("otu_table.%s.kegg.pathways.biom" % func_prefix),
            od("otu_table.%s.kegg.biom" % func_prefix),
            od("otu_table.%s.normalized.biom" % func_prefix)]

        exp_redist_fps = [
            od('otu_table.genus.redist.biom'),
            od('otu_table.species.redist.biom'),
            od('otu_table.strain.redist.biom')]

        self.assertEqual(obs_func_fps, exp_func_fps)
        self.assertEqual(obs_redist_fps, exp_redist_fps)


if __name__ == '__main__':
    main()

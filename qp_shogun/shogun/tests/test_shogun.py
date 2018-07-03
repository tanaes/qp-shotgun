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
from qp_shogun import plugin
from tempfile import mkdtemp
from json import dumps
from biom import Table
import numpy as np
from io import StringIO
from qp_shogun.shogun.utils import (
    get_dbs, get_dbs_list, generate_shogun_dflt_params,
    import_shogun_biom, shogun_db_functional_parser, shogun_parse_module_table,
    shogun_parse_enzyme_table, shogun_parse_pathway_table)
from qp_shogun.shogun.shogun import (
    generate_shogun_align_commands, _format_params,
    generate_shogun_assign_taxonomy_commands, generate_fna_file,
    generate_shogun_functional_commands, generate_shogun_redist_commands,
    shogun)

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
            'Number of threads': 5
        }
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)
        self.enzymes = ('K00001\t'
                        '"1. Oxidoreductases"\t'
                        '"1.1  Acting on the CH-OH group of donors"\t'
                        '"1.1.1  With NAD+ or NADP+ as acceptor"\t'
                        '"1.1.1.1  alcohol dehydrogenase"\n'
                        'K00002\t'
                        '"1. Oxidoreductases"\t'
                        '"1.1  Acting on the CH-OH group of donors"\t'
                        '"1.1.1  With NAD+ or NADP+ as acceptor"\t'
                        '"1.1.1.2  alcohol dehydrogenase (NADP+)"\n'
                        'K00003\t'
                        '"1. Oxidoreductases"\t'
                        '"1.1  Acting on the CH-OH group of donors"\t'
                        '"1.1.1  With NAD+ or NADP+ as acceptor"\t'
                        '"1.1.1.3  homoserine dehydrogenase"')

        self.enz_md = {
            'K00001': {'taxonomy': ['1. Oxidoreductases',
                                    '1.1  Acting on the CH-OH group of donors',
                                    '1.1.1  With NAD+ or NADP+ as acceptor',
                                    '1.1.1.1  alcohol dehydrogenase']},
            'K00002': {'taxonomy': ['1. Oxidoreductases',
                                    '1.1  Acting on the CH-OH group of donors',
                                    '1.1.1  With NAD+ or NADP+ as acceptor',
                                    '1.1.1.2  alcohol dehydrogenase (NADP+)']},
            'K00003': {'taxonomy': ['1. Oxidoreductases',
                                    '1.1  Acting on the CH-OH group of donors',
                                    '1.1.1  With NAD+ or NADP+ as acceptor',
                                    '1.1.1.3  homoserine dehydrogenase']}}

        self.modules = (
            'K00003\t'
            '"Pathway module"\t'
            '"Nucleotide and amino acid metabolism"\t'
            '"Cysteine and methionine metabolism"\t'
            '"M00017  Methionine biosynthesis,'
            ' apartate => homoserine => methionine [PATH:map00270]"\n'
            'K00003\t'
            '"Pathway module"\t'
            '"Nucleotide and amino acid metabolism"\t'
            '"Serine and threonine metabolism"\t'
            '"M00018  Threonine biosynthesis, '
            'apartate => homoserine => threonine [PATH:map00260]"\n'
            'K00133\t'
            '"Pathway module"\t'
            '"Nucleotide and amino acid metabolism"\t'
            '"Cysteine and methionine metabolism"\t'
            '"M00017  Methionine biosynthesis,'
            ' apartate => homoserine => methionine [PATH:map00270]"')

        self.mod_md = {
            'M00017': {'taxonomy': ['Pathway module',
                                    'Nucleotide and amino acid metabolism',
                                    'Cysteine and methionine metabolism',
                                    'Methionine biosynthesis,' +
                                    ' apartate => homoserine => ' +
                                    'methionine [PATH:map00270]']},
            'M00018': {'taxonomy': ['Pathway module',
                                    'Nucleotide and amino acid metabolism',
                                    'Serine and threonine metabolism',
                                    'Threonine biosynthesis,' +
                                    ' apartate => homoserine => ' +
                                    'threonine [PATH:map00260]']}}

        self.pathways = ('K00271\t'
                         '"Enzymes"\t'
                         '"1. Oxidoreductases"\t'
                         '"1.4  Acting on the CH-NH2 group of donors"\t'
                         '"1.4.1  With NAD+ or NADP+ as acceptor"\t'
                         '"1.4.1.23  valine dehydrogenase (NAD+)"\n'
                         'K00272\t'
                         '"Enzymes"\t'
                         '"1. Oxidoreductases"\t'
                         '"1.4  Acting on the CH-NH2 group of donors"\t'
                         '"1.4.3  With oxygen as acceptor"\t'
                         '"1.4.3.1  D-aspartate oxidase"\n'
                         'K00273\t'
                         '"Enzymes"\t'
                         '"1. Oxidoreductases"\t'
                         '"1.4  Acting on the CH-NH2 group of donors"\t'
                         '"1.4.3  With oxygen as acceptor"\t'
                         '"1.4.3.3  D-amino-acid oxidase"')

        self.path_md = {
            '1.4.1  With NAD+ or NADP+ as acceptor': {
                'taxonomy': ['Enzymes',
                             '1. Oxidoreductases',
                             '1.4  Acting on the CH-NH2 group of donors']},
            '1.4.3  With oxygen as acceptor': {
                'taxonomy': ['Enzymes',
                             '1. Oxidoreductases',
                             '1.4  Acting on the CH-NH2 group of donors']}}

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
                'Number of threads': 5},
            'shogun_utree': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'utree',
                'Number of threads': 5},
            'shogun_burst': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'burst',
                'Number of threads': 5}}

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

    def test_shogun_db_functional_parser(self):
        db_path = self.params['Database']
        func_prefix = 'function/ko'
        exp = {
            'enzyme': join(db_path, '%s-enzyme-annotations.txt' % func_prefix),
            'module': join(db_path, '%s-module-annotations.txt' % func_prefix),
            'pathway': join(db_path, '%s-pathway-annotations.txt'
                            % func_prefix)}
        obs = shogun_db_functional_parser(db_path)

        self.assertEqual(obs, exp)

    def test_shogun_parse_enzyme_table(self):
        out_table = shogun_parse_enzyme_table(StringIO(self.enzymes))

        self.assertDictEqual(self.enz_md, out_table)

    def test_shogun_parse_module_table(self):
        out_table = shogun_parse_module_table(StringIO(self.modules))

        self.assertDictEqual(self.mod_md, out_table)

    def test_shogun_parse_pathway_table(self):
        out_table = shogun_parse_pathway_table(StringIO(self.pathways))

        self.assertDictEqual(self.path_md, out_table)

    def test_import_shogun_biom(self):
        shogun_table = ('#OTU ID\t1450\t2563\n'
                        'k__Archaea\t26\t25\n'
                        'k__Archaea;p__Crenarchaeota\t3\t5\n'
                        'k__Archaea;p__Crenarchaeota;c__Thermoprotei\t1\t25\n')

        exp_biom = Table(np.array([[26, 25],
                                   [3, 5],
                                   [1, 25]]),
                         ['k__Archaea',
                          'k__Archaea;p__Crenarchaeota',
                          'k__Archaea;p__Crenarchaeota;c__Thermoprotei'],
                         ['1450',
                          '2563'])

        obs_biom = import_shogun_biom(StringIO(shogun_table))
        self.assertEqual(exp_biom, obs_biom)

        tax_metadata = {'k__Archaea': {
                            'taxonomy': ['k__Archaea']},
                        'k__Archaea;p__Crenarchaeota': {
                            'taxonomy': ['k__Archaea',
                                         'p__Crenarchaeota']},
                        'k__Archaea;p__Crenarchaeota;c__Thermoprotei': {
                            'taxonomy': ['k__Archaea',
                                         'p__Crenarchaeota',
                                         'c__Thermoprotei']}}
        exp_biom_tax = Table(np.array([[26, 25],
                                       [3, 5],
                                       [1, 25]]),
                             ['k__Archaea',
                              'k__Archaea;p__Crenarchaeota',
                              'k__Archaea;p__Crenarchaeota;c__Thermoprotei'],
                             ['1450',
                              '2563'])
        exp_biom_tax.add_metadata(tax_metadata, axis='observation')
        obs_biom_tax = import_shogun_biom(
            StringIO(shogun_table), names_to_taxonomy=True)

        self.assertEqual(exp_biom_tax, obs_biom_tax)

        # test modules
        module_table = ('#MODULE ID\t1450\t2563\n'
                        'M00017\t26\t25\n'
                        'M00018\t3\t5\n')

        exp_m_biom = Table(np.array([[26, 25],
                                     [3, 5]]),
                           ['M00017', 'M00018'],
                           ['1450', '2563'])
        exp_m_biom.add_metadata(self.mod_md, axis='observation')
        obs_m_biom = import_shogun_biom(
            StringIO(module_table), annotation_table=StringIO(self.modules),
            annotation_type='module')

        self.assertEqual(exp_m_biom, obs_m_biom)

        # test pathways
        path_table = ('#PATHWAY ID\t1450\t2563\n'
                      '1.4.1  With NAD+ or NADP+ as acceptor\t26\t25\n'
                      '1.4.3  With oxygen as acceptor\t3\t5\n')

        exp_p_biom = Table(np.array([[26, 25],
                                     [3, 5]]),
                           ['1.4.1  With NAD+ or NADP+ as acceptor',
                            '1.4.3  With oxygen as acceptor'],
                           ['1450', '2563'])

        exp_p_biom.add_metadata(self.path_md, axis='observation')
        obs_p_biom = import_shogun_biom(
            StringIO(path_table), annotation_table=StringIO(self.pathways),
            annotation_type='pathway')

        self.assertEqual(exp_p_biom, obs_p_biom)

        # test enzymes
        enzyme_table = ('#KEGG ID\t1450\t2563\n'
                        'K00001\t26\t25\n'
                        'K00002\t3\t5\n'
                        'K00003\t1\t25\n')
        exp_e_biom = Table(np.array([[26, 25],
                                     [3, 5],
                                     [1, 25]]),
                           ['K00001',
                            'K00002',
                            'K00003'],
                           ['1450', '2563'])
        exp_e_biom.add_metadata(self.enz_md, axis='observation')
        obs_e_biom = import_shogun_biom(
            StringIO(enzyme_table), annotation_table=StringIO(self.enzymes),
            annotation_type='enzyme')

        self.assertEqual(exp_e_biom, obs_e_biom)

        # test empty
        empty_table = ('#KEGG ID\t1450\t2563\n')
        exp_empty_biom = Table(np.zeros((0, 2)),
                               [],
                               ['1450', '2563'])
        obs_empty_biom = import_shogun_biom(
            StringIO(empty_table), annotation_table=StringIO(self.enzymes),
            annotation_type='enzyme')

        self.assertEqual(exp_empty_biom, obs_empty_biom)

    def test_format_shogun_params(self):
        obs = _format_params(self.params, SHOGUN_PARAMS)
        exp = {
            'database': join(self.db_path, 'shogun'),
            'aligner': 'bowtie2',
            'threads': 5
        }

        self.assertEqual(obs, exp)

    def test_generate_shogun_align_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:

            exp_cmd = [
                ('shogun align --aligner bowtie2 --threads 5 '
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

    # Testing shogun with bowtie2
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
                'command': dumps(['qp-shogun', '0.1.0', 'Shogun']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 1 artifacts in total
        self.assertEqual(1, len(ainfo))
        ainfo = ainfo[0]
        self.assertEqual(ainfo.artifact_type, 'BIOM')
        exp = [(join(out_dir, 'otu_table.alignment.profile.biom'), 'biom')]
        self.assertCountEqual(ainfo.files, exp)

    def test_shogun_burst(self):
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
        self.params['Aligner tool'] = 'burst'
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shogun', '0.1.0', 'Shogun']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 1 artifacts in total
        self.assertEqual(1, len(ainfo))
        ainfo = ainfo[0]
        self.assertEqual(ainfo.artifact_type, 'BIOM')
        exp = [(join(out_dir, 'otu_table.alignment.profile.biom'), 'biom')]
        self.assertCountEqual(ainfo.files, exp)

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
                'command': dumps(['qp-shogun', '0.1.0', 'Shogun']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = shogun(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 1 artifacts in total
        self.assertEqual(1, len(ainfo))
        ainfo = ainfo[0]
        self.assertEqual(ainfo.artifact_type, 'BIOM')
        exp = [(join(out_dir, 'otu_table.alignment.profile.biom'), 'biom')]
        self.assertCountEqual(ainfo.files, exp)


if __name__ == '__main__':
    main()

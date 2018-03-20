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
from shutil import rmtree
from tempfile import TemporaryDirectory
from qp_shotgun import plugin
from tempfile import mkdtemp
from qp_shotgun.shogun.utils import (
    get_dbs, get_dbs_list, generate_shogun_dflt_params)
from qp_shotgun.shogun.shogun import (
    generate_shogun_align_commands, _format_params,
    generate_shogun_assign_taxonomy_commands, generate_fna_file,
    generate_shogun_functional_commands, generate_biom_conversion_commands)

SHOGUN_PARAMS = {
    'Database': 'database', 'Aligner tool': 'aligner',
    'Taxonomic Level': 'levels', 'Number of threads': 'threads'}


class ShogunTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.out_dir = out_dir
        self.db_path = os.environ["QC_SHOGUN_DB_DP"]
        self.params = {
            'Database': join(self.db_path, 'shogun'),
            'Aligner tool': 'bowtie2',
            'Taxonomic Level': 'all',
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
                'Taxonomy Level': 'all',
                'Number of threads': 1},
            'shogun_utree': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'utree',
                'Taxonomy Level': 'all',
                'Number of threads': 1},
            'shogun_burst': {
                'Database': join(self.db_path, 'shogun'),
                'Aligner tool': 'burst',
                'Taxonomy Level': 'all',
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
            'levels': 'all',
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
                 join(temp_dir, 'profile.functional.species.tsv'))
                ]
            profile_dir = join(temp_dir, 'profile.tsv')
            params = _format_params(self.params, SHOGUN_PARAMS)
            obs_cmd = generate_shogun_functional_commands(
                profile_dir, temp_dir, params, 'species')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_generate_biom_conversion_commands(self):
        out_dir = self.out_dir
        with TemporaryDirectory(dir=out_dir, prefix='shogun_') as temp_dir:
            exp_cmd = [
                ('biom convert -i %s '
                 '-o %s '
                 '--table-type="OTU table" '
                 '--process-obs-metadata taxonomy --to-hdf5') %
                (join(temp_dir, 'profile.tsv'),
                 join(out_dir, 'otu_table.species.biom'))
                ]
            profile_fp = join(temp_dir, 'profile.tsv')
            obs_cmd = generate_biom_conversion_commands(profile_fp, out_dir,
                                                        'species')

        self.assertEqual(obs_cmd, exp_cmd)

    def test_shogun(self):
        pass


if __name__ == '__main__':
    main()

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps
from os.path import exists, isdir, basename

from qiita_client.testing import PluginTestCase


from qp_shotgun.humann2 import plugin
from qp_shotgun.humann2.humann2 import (
    humann2, generate_humann2_analysis_commands)


class Humann2Tests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'nucleotide-database': 'chocophlan', 'protein-database': 'uniref',
            'bypass-prescreen': False, 'bypass-nucleotide-index': False,
            'bypass-translated-search': False,
            'bypass-nucleotide-search': False,
            'annotation-gene-index': 8, 'evalue': 1.0, 'search-mode': '',
            'metaphlan-options': '-t rel_ab', 'log-level': 'DEBUG',
            'remove-temp-output': False, 'threads': 1,
            'prescreen-threshold': 0.01, 'identity-threshold': 50.0,
            'translated-subject-coverage-threshold': 50.0,
            'translated-query-coverage-threshold': 90.0,
            'translated-alignment': 'diamond',
            'xipe': 'off', 'minpath': 'on', 'pick-frames': 'off',
            'gap-fill': 'off', 'minpath': 'on', 'output-format': 'biom',
            'output-max-decimals': 10, 'remove-stratified-output': False,
            'input-format': '', 'pathways': 'metacyc', 'memory-use': 'minimum'}
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_generate_humann2_analysis_commands_names_no_match(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        with self.assertRaises(ValueError):
            generate_humann2_analysis_commands(['a', 'b', 'c'], [], fp,
                                               'output', {})

    def test_generate_humann2_analysis_commands_fwd_rev_not_match(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        with self.assertRaises(ValueError):
            generate_humann2_analysis_commands(['s1', 's2', 's3'], ['a'], fp,
                                               'output', {})

    def test_generate_humann2_analysis_commands_only_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = [
            'humann2 --input "fastq/s1.fastq" --output "output/s1" '
            '--output-basename "SKB8.640193" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"',
            'humann2 --input "fastq/s2.fastq.gz" --output "output/s2" '
            '--output-basename "SKD8.640184" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"',
            'humann2 --input "fastq/s3.fastq" --output "output/s3" '
            '--output-basename "SKB7.640196" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"']
        obs = generate_humann2_analysis_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'], [],
            fp, 'output', self.params)
        self.assertEqual(obs, exp)

    def test_generate_humann2_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = [
            'humann2 --input "fastq/s1.fastq" --output "output/s1" '
            '--output-basename "SKB8.640193" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"',
            'humann2 --input "fastq/s1.R2.fastq" --output "output/s1.R2" '
            '--output-basename "SKB8.640193" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" '
            '--pathways "metacyc" --pick-frames "off" '
            '--translated-alignment "diamond" --log-level "DEBUG"',
            'humann2 --input "fastq/s2.fastq.gz" --output "output/s2" '
            '--output-basename "SKD8.640184" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" --protein-database '
            '"uniref" --threads "1" --pathways "metacyc" --pick-frames "off" '
            '--translated-alignment "diamond" --log-level "DEBUG"',
            'humann2 --input "fastq/s2.R2.fastq.gz" --output "output/s2.R2" '
            '--output-basename "SKD8.640184" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" '
            '--evalue "1.0" --minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"',
            'humann2 --input "fastq/s3.fastq" --output "output/s3" '
            '--output-basename "SKB7.640196" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" '
            '--evalue "1.0" --minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"',
            'humann2 --input "fastq/s3.R2.fastq" --output "output/s3.R2" '
            '--output-basename "SKB7.640196" --output-format biom '
            '--gap-fill "off" '
            '--identity-threshold "50.0" --output-format "biom" '
            '--metaphlan-options "-t rel_ab" '
            '--translated-query-coverage-threshold "90.0" '
            '--prescreen-threshold "0.01" '
            '--translated-subject-coverage-threshold "50.0" --evalue "1.0" '
            '--minpath "on" --output-max-decimals "10" '
            '--nucleotide-database "chocophlan" --memory-use "minimum" '
            '--xipe "off" --annotation-gene-index "8" '
            '--protein-database "uniref" --threads "1" --pathways "metacyc" '
            '--pick-frames "off" --translated-alignment "diamond" '
            '--log-level "DEBUG"']
        obs = generate_humann2_analysis_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output', self.params)
        self.assertEqual(obs, exp)

    def test_humann2(self):
        # generating filepaths
        fd, fp1 = mkstemp(prefix='demo_SKB7_', suffix='_seqs.fastq.gz')
        close(fd)
        self._clean_up_files.append(fp1)
        fd, fp2 = mkstemp(prefix='demo_SKB8_', suffix='_seqs.fastq.gz')
        close(fd)
        self._clean_up_files.append(fp2)
        copyfile('support_files/demo.fastq.gz', fp1)
        copyfile('support_files/demo.fastq.gz', fp2)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {
                'run_prefix': basename(fp1).replace('.fastq.gz', '')},
            'SKB8.640193': {
                'run_prefix': basename(fp2).replace('.fastq.gz', '')}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1, 'raw_forward_seqs'),
                (fp2, 'raw_forward_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "New test artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        # overwriting so the run uses the defaults
        self.params['nucleotide-database'] = ''
        self.params['protein-database'] = ''
        data = {'user': 'demo@microbio.me',
                'command': dumps(['HUMAnN2', '0.9.1', 'HUMAnN2']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        humann2(self.qclient, jid, self.params, out_dir)

        # TODO: test that the files are created properly


MAPPING_FILE = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\trun_prefix\t"
    "instrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts2\tIllumina MiSeq\tdesc3\n"
)

MAPPING_FILE_2 = (
    "#SampleID\tplatform\tbarcode\texperiment_design_description\t"
    "library_construction_protocol\tcenter_name\tprimer\t"
    "run_prefix\tinstrument_model\tDescription\n"
    "SKB7.640196\tILLUMINA\tA\tA\tA\tANL\tA\ts3\tIllumina MiSeq\tdesc1\n"
    "SKB8.640193\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc2\n"
    "SKD8.640184\tILLUMINA\tA\tA\tA\tANL\tA\ts1\tIllumina MiSeq\tdesc3\n"
)


if __name__ == '__main__':
    main()

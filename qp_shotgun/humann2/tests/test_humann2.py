# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from os import close, environ
from tempfile import mkstemp


from qiita_client import QiitaClient
# , ArtifactInfo


from qp_shotgun.humann2 import (
    get_sample_names_by_run_prefix, generate_humann2_analysis_commands)


CLIENT_ID = '19ndkO3oMKsoChjVVWluF7QkxHRfYhTKSFbAVt8IhK7gZgDaO4'
CLIENT_SECRET = ('J7FfQ7CQdOxuKhQAf1eoGgBAE81Ns8Gu3EKaWFm3IO2JKh'
                 'AmmCWZuabe0O5Mp28s1')


class Humann2Tests(TestCase):
    @classmethod
    def setUpClass(cls):
        server_cert = environ.get('QIITA_SERVER_CERT', None)
        cls.qclient = QiitaClient("https://localhost:21174", CLIENT_ID,
                                  CLIENT_SECRET, server_cert=server_cert)
        cls._clean_up_files = []

    @classmethod
    def tearDownClass(cls):
        cls.qclient.post('/apitest/reset/')

    def test_get_sample_names_by_run_prefix(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        obs = get_sample_names_by_run_prefix(fp)
        exp = {'s3': 'SKB7.640196', 's2': 'SKD8.640184', 's1': 'SKB8.640193'}
        self.assertEqual(obs, exp)

    def test_get_sample_names_by_run_prefix_error(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE_2)
        self._clean_up_files.append(fp)

        with self.assertRaises(ValueError):
            get_sample_names_by_run_prefix(fp)

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

        params = {"--nucleotide-database": "chocophlan",
                  "--protein-database": "uniref"}
        exp = [
            'humann2 --input fastq/s1.fastq --output output/s1 '
            '--output-basename SKB8.640193 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s2.fastq.gz --output output/s2 '
            '--output-basename SKD8.640184 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s3.fastq --output output/s3 '
            '--output-basename SKB7.640196 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan']
        obs = generate_humann2_analysis_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'], [],
            fp, 'output', params)
        self.assertEqual(obs, exp)

    def test_generate_humann2_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        params = {"--nucleotide-database": "chocophlan",
                  "--protein-database": "uniref"}
        exp = [
            'humann2 --input fastq/s1.fastq --output output/s1 '
            '--output-basename SKB8.640193 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s1.R2.fastq --output output/s1.R2 '
            '--output-basename SKB8.640193 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s2.fastq.gz --output output/s2 '
            '--output-basename SKD8.640184 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s2.R2.fastq.gz --output output/s2.R2 '
            '--output-basename SKD8.640184 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s3.fastq --output output/s3 '
            '--output-basename SKB7.640196 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan',
            'humann2 --input fastq/s3.R2.fastq --output output/s3.R2 '
            '--output-basename SKB7.640196 --output-format biom '
            '--protein-database uniref --nucleotide-database chocophlan']
        obs = generate_humann2_analysis_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output', params)
        self.assertEqual(obs, exp)

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

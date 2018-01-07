# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove, makedirs
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial
import os

from qiita_client.testing import PluginTestCase

from qp_shotgun import plugin
from qp_shotgun.qc_filter.qc_filter import (qc_filter,
                                            make_read_pairs_per_sample,
                                            _format_qc_filter_params,
                                            generate_qc_filter_commands)
import qp_shotgun.qc_filter as kd


class QC_FilterTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        self.params = {
                       'Bowtie2 database to filter': 'Human',
                       'Number of threads to be used': '4'
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_make_read_pairs_per_sample_match_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s3_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz',
                  './folder/s1_S009_L001_R2.fastq.gz']

        exp = [('s1', 'SKB8.640193', './folder/s1_S009_L001_R1.fastq.gz',
                './folder/s1_S009_L001_R2.fastq.gz'),
               ('s2', 'SKD8.640184', './folder/s2_S011_L001_R1.fastq.gz',
                './folder/s2_S011_L001_R2.fastq.gz'),
               ('s3', 'SKB7.640196', './folder/s3_S013_L001_R1.fastq.gz',
                './folder/s3_S013_L001_R2.fastq.gz')]

        obs = make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

        self.assertEqual(obs, exp)

    def test_make_read_pairs_per_sample_match_fwd_only(self):
        # do we want to support forward-only yet?
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = []

        exp = [('s1', 'SKB8.640193', './folder/s1_S009_L001_R1.fastq.gz',
                None),
               ('s2', 'SKD8.640184', './folder/s2_S011_L001_R1.fastq.gz',
                None),
               ('s3', 'SKB7.640196', './folder/s3_S013_L001_R1.fastq.gz',
                None)]

        obs = make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

        self.assertEqual(obs, exp)

    def test_make_read_pairs_per_sample_match_fwd_rev_diffnum(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s4_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz',
                  './folder/s1_S009_L001_R2.fastq.gz']

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_rev_notmatch(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s1_S009_L001_R1.fastq.gz']

        rev_fp = ['./folder/s3_S013_L001_R2.fastq.gz',
                  './folder/s2_S011_L001_R2.fastq.gz']

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_no_rp(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s4_S009_L001_R1.fastq.gz']

        rev_fp = []

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_make_read_pairs_per_sample_match_fwd_2match(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001_R1.fastq.gz',
                  './folder/s2_S011_L001_R1.fastq.gz',
                  './folder/s2_S009_L001_R1.fastq.gz']

        rev_fp = []

        with self.assertRaises(ValueError):
            make_read_pairs_per_sample(fwd_fp, rev_fp, fp)

    def test_format_qc_filter_params(self):
        db_path = os.environ["QC_FILTER_DB_DP"]
        obs = _format_qc_filter_params(self.params)
        exp = ('-p 4 -x %sHuman/phix') % db_path
        self.assertEqual(obs, exp)

    def test_generate_qc_trim_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        db_path = os.environ["QC_FILTER_DB_DP"]

        exp_cmd = [
            ('bowtie2 -p 4 -x %sHuman/phix --very-sensitive -1 fastq/s1.fastq.gz '
            '-2 fastq/s1.R2.fastq.gz | '
            'samtools view -f 12 -F 256 -b -o temp/SKB8.640193.unsorted.bam'

            'samtools sort -T temp/SKB8.640193 -@ 4 -n -o temp/SKB8.640193.bam '
            'temp/SKB8.640193.unsorted.bam'

            'bedtools bamtofastq -i temp/SKB8.640193.bam -fq '
            'temp/SKB8.640193.R1.trimmed.filtered.fastq -fq2 '
            'temp/SKB8.640193.R2.trimmed.filtered.fastq'

            'pigz -p 4 -c temp/SKB8.640193.R1.trimmed.filtered.fastq > '
            'output/SKB8.640193.R1.trimmed.filtered.fastq.gz'
            'pigz -p 4 -c temp/SKB8.640193.R2.trimmed.filtered.fastq > '
            'output/SKB8.640193.R2.trimmed.filtered.fastq.gz') % db_path
            ]

        exp_sample = [
            ('s1', 'SKB8.640193', 'fastq/s1.fastq.gz', 'fastq/s1.R2.fastq.gz')
            #('s2', 'SKD8.640184', 'fastq/s2.fastq.gz', 'fastq/s2.R2.fastq.gz'),
            #('s3', 'SKB7.640196', 'fastq/s3.fastq.gz', 'fastq/s3.R2.fastq.gz')
            ]

        obs_cmd, obs_sample = generate_qc_filter_commands(
            ['fastq/s1.fastq.gz'],
            ['fastq/s1.R2.fastq.gz'],
            fp, 'output', 'temp', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)


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

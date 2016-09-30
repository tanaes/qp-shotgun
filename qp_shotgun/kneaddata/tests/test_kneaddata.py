# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove
from os.path import exists, isdir
from shutil import rmtree
from tempfile import mkstemp

from qiita_client.testing import PluginTestCase

from qp_shotgun.kneaddata import plugin
from qp_shotgun.kneaddata.kneaddata import (make_read_pairs_per_sample,
                                            generate_kneaddata_commands)


class KneaddataTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {}
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

    def test_generate_humann2_analysis_commands_only_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = [
            'kneaddata --input "fastq/s1.fastq" --output "output/s1" '
            '--output-prefix "s1" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"',
            'kneaddata --input "fastq/s2.fastq" --output "output/s2" '
            '--output-prefix "s2" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"',
            'kneaddata --input "fastq/s3.fastq" --output "output/s3" '
            '--output-prefix "s3" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"']

        obs = generate_kneaddata_commands(
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
            'kneaddata --input "fastq/s1.fastq" --input "fastq/s1.R2.fastq"'
            ' --output "output/s1" '
            '--output-prefix "s1" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"',
            'kneaddata --input "fastq/s2.fastq" --input "fastq/s2.R2.fastq" '
            ' --output "output/s2" '
            '--output-prefix "s2" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"',
            'kneaddata --input "fastq/s3.fastq" --input "fastq/s3.R2.fastq" '
            ' --output "output/s3" '
            '--output-prefix "s3" --reference-db human_genome --threads 1 '
            '--processes 1 --quality-scores phred33 --run-fastqc-start '
            '--run-fastqc-end --log-level DEBUG --max-memory 500 '
            '--trimmomatic-options "ILLUMINACLIP:$trimmomatic/adapters/'
            'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36" --bowtie2-options "--very-sensitive"']

        obs = generate_kneaddata_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output', self.params)
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

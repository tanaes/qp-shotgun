# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial

from qiita_client.testing import PluginTestCase

from qp_shotgun import plugin
from qp_shotgun.kneaddata.kneaddata import (make_read_pairs_per_sample,
                                            generate_kneaddata_commands,
                                            _format_kneaddata_params,
                                            kneaddata)


class KneaddataTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'reference-db': 'default', 'bypass-trim': False, 'threads': 1,
            'processes': 1, 'quality-scores': 'phred33', 'run-bmtagger': False,
            'run-trf': False, 'run-fastqc-start': True,
            'run-fastqc-end': False, 'max-memory': '500m',
            'trimmomatic-options': (
                '"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"')
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_kneaddata_params(self):
        obs = _format_kneaddata_params(self.params)
        exp = ('--max-memory 500m --processes 1 --quality-scores phred33 '
               '--run-fastqc-start --threads 1 '
               '--trimmomatic-options "LEADING:3 '
               'TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"')

        self.assertEqual(obs, exp)

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

    def test_generate_kneaddata_analysis_commands_only_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp_cmd = [
            'kneaddata --input "fastq/s1.fastq" '
            '--output "output/s1" --output-prefix "s1" '
            '--max-memory 500m --processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 '
            '--trimmomatic-options "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36"',
            'kneaddata --input "fastq/s2.fastq.gz" '
            '--output "output/s2" --output-prefix "s2" '
            '--max-memory 500m --processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 '
            '--trimmomatic-options "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36"',
            'kneaddata --input "fastq/s3.fastq" '
            '--output "output/s3" --output-prefix "s3" '
            '--max-memory 500m --processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 '
            '--trimmomatic-options "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36"']

        exp_sample = [('s1', 'SKB8.640193', 'fastq/s1.fastq', None),
                      ('s2', 'SKD8.640184', 'fastq/s2.fastq.gz', None),
                      ('s3', 'SKB7.640196', 'fastq/s3.fastq', None)]

        obs_cmd, obs_sample = generate_kneaddata_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'], [],
            fp, 'output', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)

    def test_generate_kneaddata_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp_cmd = [
            'kneaddata --input "fastq/s1.fastq" --input "fastq/s1.R2.fastq" '
            '--output "output/s1" --output-prefix "s1" '
            '--max-memory 500m --processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 '
            '--trimmomatic-options "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36"',
            'kneaddata --input "fastq/s2.fastq.gz" '
            '--input "fastq/s2.R2.fastq.gz" --output "output/s2" '
            '--output-prefix "s2" --max-memory 500m '
            '--processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 --trimmomatic-options '
            '"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"',
            'kneaddata --input "fastq/s3.fastq" --input "fastq/s3.R2.fastq" '
            '--output "output/s3" --output-prefix "s3" '
            '--max-memory 500m --processes 1 --quality-scores phred33 '
            '--run-fastqc-start --threads 1 '
            '--trimmomatic-options "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 '
            'MINLEN:36"']

        exp_sample = [
            ('s1', 'SKB8.640193', 'fastq/s1.fastq', 'fastq/s1.R2.fastq'),
            ('s2', 'SKD8.640184', 'fastq/s2.fastq.gz', 'fastq/s2.R2.fastq.gz'),
            ('s3', 'SKB7.640196', 'fastq/s3.fastq', 'fastq/s3.R2.fastq')]

        obs_cmd, obs_sample = generate_kneaddata_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)

    def test_kneaddata(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'kd_test_1_R1.fastq.gz')
        fp1_2 = join(in_dir, 'kd_test_1_R2.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp1_2)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {
                'run_prefix': 'kd_test_1'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test KneadData artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shotgun', '0.0.1', 'KneadData 0.5.1']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = kneaddata(self.qclient, jid,
                                        self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)
        # we are expecting 3 artifacts in total
        self.assertEqual(3, len(ainfo))

        obs_fps = []
        for a in ainfo:
            self.assertEqual("per_sample_FASTQ", a.artifact_type)
            obs_fps.append(a.files)
        od = partial(join, out_dir, 'kd_test_1')
        exp_fps = [
            [(od('kd_test_1_paired_1.fastq'), 'preprocessed_fastq'),
             (od('kd_test_1_paired_2.fastq'), 'preprocessed_fastq')],
            [(od('kd_test_1_unmatched_1.fastq'), 'preprocessed_fastq')],
            [(od('kd_test_1_unmatched_2.fastq'), 'preprocessed_fastq')]]
        self.assertItemsEqual(exp_fps, obs_fps)


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

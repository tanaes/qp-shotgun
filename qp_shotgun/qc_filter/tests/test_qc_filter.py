# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove, makedirs
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import mkstemp, mkdtemp
from json import dumps
from functools import partial
import os

from qiita_client.testing import PluginTestCase

from qp_shotgun import plugin
from qp_shotgun.qc_filter.qc_filter import (make_read_pairs_per_sample,
                                            _format_qc_filter_params,
                                            generate_qc_filter_commands,
                                            qc_filter, _per_sample_ainfo)


class QC_FilterTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        self.params = {
                       'Bowtie2 database to filter': 'Human',
                       'Number of threads to be used': '1'
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_qc_filter_params(self):
        db_path = os.environ["QC_FILTER_DB_DP"]
        obs = _format_qc_filter_params(self.params)
        exp = ('-p 1 -x %sHuman/phix') % db_path
        self.assertEqual(obs, exp)

    def test_generate_qc_filter_analysis_commands_forward_reverse(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)
        db_path = os.environ["QC_FILTER_DB_DP"]

        exp_cmd = [
            ('bowtie2 -p 1 -x %sHuman/phix --very-sensitive '
             '-1 fastq/s1.fastq.gz -2 fastq/s1.R2.fastq.gz | '
             'samtools view -f 12 -F 256 -b -o temp/SKB8.640193.unsorted.bam; '

             'samtools sort -T temp/SKB8.640193 -@ 1 -n '
             '-o temp/SKB8.640193.bam temp/SKB8.640193.unsorted.bam; '

             'bedtools bamtofastq -i temp/SKB8.640193.bam -fq '
             'temp/SKB8.640193.R1.trimmed.filtered.fastq -fq2 '
             'temp/SKB8.640193.R2.trimmed.filtered.fastq; '

             'pigz -p 1 -c temp/SKB8.640193.R1.trimmed.filtered.fastq > '
             'output/SKB8.640193.R1.trimmed.filtered.fastq.gz; '
             'pigz -p 1 -c temp/SKB8.640193.R2.trimmed.filtered.fastq > '
             'output/SKB8.640193.R2.trimmed.filtered.fastq.gz;') % db_path
            ]

        exp_sample = [
            ('s1', 'SKB8.640193', 'fastq/s1.fastq.gz', 'fastq/s1.R2.fastq.gz')
            ]

        obs_cmd, obs_sample = generate_qc_filter_commands(
            ['fastq/s1.fastq.gz'],
            ['fastq/s1.R2.fastq.gz'],
            fp, 'output', 'temp', self.params)

        self.assertEqual(obs_cmd, exp_cmd)
        self.assertEqual(obs_sample, exp_sample)

    def test_qc_filter(self):
        # generating filepaths
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'kd_test_1_R1.fastq.gz')
        fp1_2 = join(in_dir, 'kd_test_1_R2.fastq.gz')
        fp2_1 = join(in_dir, 'kd_test_2_R1.fastq.gz')
        fp2_2 = join(in_dir, 'kd_test_2_R2.fastq.gz')
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp1_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp1_2)
        copyfile('support_files/kd_test_1_R1.fastq.gz', fp2_1)
        copyfile('support_files/kd_test_1_R2.fastq.gz', fp2_2)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'run_prefix': 'kd_test_1'},
            'SKB8.640193': {'run_prefix': 'kd_test_2'}
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
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test QC_Trim artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shotgun', '0.0.1', 'QC_Filter']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = qc_filter(self.qclient, jid,
                                        self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # we are expecting 3 artifacts in total
        self.assertEqual(1, len(ainfo))

        obs_fps = []
        for a in ainfo:
            self.assertEqual("per_sample_FASTQ", a.artifact_type)
            obs_fps.append(a.files)
        od = partial(join, out_dir)

        # ftype = 'per_sample_FASTQ'
        exp_fps = [
            [od('1.SKB7.640196.R1.trimmed.filtered.fastq.gz'),
             od('1.SKB7.640196.R2.trimmed.filtered.fastq.gz'),
             od('1.SKB8.640193.R1.trimmed.filtered.fastq.gz'),
             od('1.SKB8.640193.R2.trimmed.filtered.fastq.gz')]]
        self.assertEqual(exp_fps, obs_fps)

    def test_per_sample_ainfo_error(self):
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)
        makedirs(join(in_dir, 'sampleA'))
        makedirs(join(in_dir, 'sampleB'))

        # Paired-end
        with self.assertRaises(ValueError):
            _per_sample_ainfo(in_dir, (('sampleA', None, None, None),
                                       ('sampleB', None, None, None)), True)


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

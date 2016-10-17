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

from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo

from qp_shotgun.fastqc import plugin
from qp_shotgun.fastqc.fastqc import (format_fastqc_params,
                                      generate_fastqc_commands,
                                      _guess_fastqc_filename,
                                      _per_sample_ainfo,
                                      fastqc)


class FastQCTests(PluginTestCase):
    maxDiff = None

    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'extract': False, 'noextract': True, 'threads': 1, 'kmers': 7
        }
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_format_fastqc_params(self):
        obs = format_fastqc_params(self.params)
        exp = ('--kmers 7 --noextract --threads 1')

        self.assertEqual(obs, exp)

    def test__guess_fastqc_filename(self):
        obs = _guess_fastqc_filename('./folder/file1.R1.fastq.gz')
        exp = ('file1.R1_fastqc.html', 'file1.R1_fastqc.zip')

        self.assertEqual(obs,exp)

    def test__per_sample_ainfo(self):
        out_dir = './folder'
        samples = [('s1', 'SKB8.640193', './folder/s1_S009_L001_R1.fastq.gz',
                './folder/s1_S009_L001_R2.fastq.gz'),
               ('s2', 'SKD8.640184', './folder/s2_S011_L001_R1.fastq.gz',
                './folder/s2_S011_L001_R2.fastq.gz')]

        exp = [ArtifactInfo('FastQC html summary', 'html_summary',
                         [('./folder/s1/s1_S009_L001_R1_fastqc.html',
                          'html_summary')]),
               ArtifactInfo('FastQC data summary', 'zip_file',
                         [('./folder/s1/s1_S009_L001_R1_fastqc.zip',
                          'zip_file')]),
               ArtifactInfo('FastQC html summary', 'html_summary',
                         [('./folder/s1/s1_S009_L001_R2_fastqc.html',
                          'html_summary')]),
               ArtifactInfo('FastQC data summary', 'zip_file',
                         [('./folder/s1/s1_S009_L001_R2_fastqc.zip',
                          'zip_file')]),
               ArtifactInfo('FastQC html summary', 'html_summary',
                         [('./folder/s2/s2_S011_L001_R1_fastqc.html',
                          'html_summary')]),
               ArtifactInfo('FastQC data summary', 'zip_file',
                         [('./folder/s2/s2_S011_L001_R1_fastqc.zip',
                          'zip_file')]),
               ArtifactInfo('FastQC html summary', 'html_summary',
                         [('./folder/s2/s2_S011_L001_R2_fastqc.html',
                          'html_summary')]),
               ArtifactInfo('FastQC data summary', 'zip_file',
                         [('./folder/s2/s2_S011_L001_R2_fastqc.zip',
                          'zip_file')])]

        obs = _per_sample_ainfo(out_dir, samples)

        for i, a in enumerate(exp):
            self.assertEqual(exp[i],obs[i])

        #self.assertItemsEqual(obs, exp)

    def test_generate_fastqc_commands_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = ['mkdir -p output/s1; fastqc --outdir "output/s1" --kmers 7 --noextract --threads 1 '
               'fastq/s1.fastq fastq/s1.R2.fastq',
               'mkdir -p output/s2; fastqc --outdir "output/s2" --kmers 7 --noextract --threads 1 '
               'fastq/s2.fastq.gz fastq/s2.R2.fastq.gz',
               'mkdir -p output/s3; fastqc --outdir "output/s3" --kmers 7 --noextract --threads 1 '
               'fastq/s3.fastq fastq/s3.R2.fastq']

        obs, samp = generate_fastqc_commands(
            ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'],
            ['fastq/s1.R2.fastq', 'fastq/s2.R2.fastq.gz', 'fastq/s3.R2.fastq'],
            fp, 'output', self.params)

        self.assertEqual(obs, exp)

    def test_generate_fastqc_commands_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        exp = ['mkdir -p output/s1; fastqc --outdir "output/s1" --kmers 7 --noextract --threads 1 '
               'fastq/s1.fastq',
               'mkdir -p output/s2; fastqc --outdir "output/s2" --kmers 7 --noextract --threads 1 '
               'fastq/s2.fastq.gz',
               'mkdir -p output/s3; fastqc --outdir "output/s3" --kmers 7 --noextract --threads 1 '
               'fastq/s3.fastq']

        obs, samp = generate_fastqc_commands(
                ['fastq/s1.fastq', 'fastq/s2.fastq.gz', 'fastq/s3.fastq'], [],
                fp, 'output', self.params)

        self.assertEqual(obs, exp)

    def test_fastqc(self):
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
            'name': "New test artifact",
            'prep': pid}

        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['FastQC', '0.11.5', 'FastQC']),
                'status': 'running',
                'parameters': dumps(self.params)}

        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = fastqc(self.qclient, jid,
                                        self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)
        # we are expecting 16 artifacts per sample eventually
        # but for now just the four fastqs
        self.assertEqual(4, len(ainfo))

        obs_fps = []
        obs_arts = []
        for a in ainfo:
            obs_arts.append(a.artifact_type)
            obs_fps.append(a.files)
        self.assertEqual({'zip_file','html_summary'}, set(obs_arts))

        exp_fps = [[(join(out_dir, 'kd_test_1', 'kd_test_1_R1_fastqc.html'),
                    'html_summary')],
                   [(join(out_dir, 'kd_test_1', 'kd_test_1_R1_fastqc.zip'),
                    'zip_file')],
                   [(join(out_dir, 'kd_test_1', 'kd_test_1_R2_fastqc.html'),
                    'html_summary')],
                   [(join(out_dir, 'kd_test_1', 'kd_test_1_R2_fastqc.zip'),
                    'zip_file')]]

        self.assertItemsEqual(exp_fps, obs_fps)

        for f_a in exp_fps:
            assert exists(f_a[0][0])

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

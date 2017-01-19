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
from filecmp import cmp as fcmp
from os.path import exists, isdir, basename, join

from qiita_client.testing import PluginTestCase

from qp_shotgun import plugin
from qp_shotgun.humann2.humann2 import (
    make_read_sets_per_sample, make_single_fastq_gz,
    humann2, generate_humann2_analysis_commands)


class Humann2Tests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'read-set': 'fwd',
            'nucleotide-database': 'default', 'protein-database': 'default',
            'bypass-prescreen': False, 'bypass-nucleotide-index': False,
            'bypass-translated-search': False,
            'bypass-nucleotide-search': False,
            'annotation-gene-index': 8, 'evalue': 1.0,
            'metaphlan-options': '-t rel_ab', 'log-level': 'DEBUG',
            'remove-temp-output': False, 'threads': 1,
            'prescreen-threshold': 0.01, 'identity-threshold': 50.0,
            'translated-subject-coverage-threshold': 50.0,
            'translated-query-coverage-threshold': 90.0,
            'translated-alignment': 'diamond',
            'xipe': 'off', 'minpath': 'on', 'pick-frames': 'off',
            'gap-fill': 'off', 'output-format': 'biom',
            'output-max-decimals': 10, 'remove-stratified-output': False,
            'pathways': 'metacyc',
            'memory-use': 'minimum', 'remove-column-description-output': True}
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_make_read_sets_per_sample_match_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_p_fp = ['./folder/s3_paired_1.fastq.gz',
                    './folder/s2_paired_1.fastq.gz',
                    './folder/s1_paired_1.fastq.gz']

        rev_p_fp = ['./folder/s3_paired_2.fastq.gz',
                    './folder/s2_paired_2.fastq.gz',
                    './folder/s1_paired_2.fastq.gz']

        fwd_u_fp = ['./folder/s3_unmatched_1.fastq.gz',
                    './folder/s2_unmatched_1.fastq.gz',
                    './folder/s1_unmatched_1.fastq.gz']

        rev_u_fp = ['./folder/s3_unmatched_2.fastq.gz',
                    './folder/s2_unmatched_2.fastq.gz',
                    './folder/s1_unmatched_2.fastq.gz']

        files = fwd_p_fp + rev_p_fp + fwd_u_fp + rev_u_fp

        exp = [('s1', 'SKB8.640193',
                './folder/s1_paired_1.fastq.gz',
                './folder/s1_paired_2.fastq.gz',
                './folder/s1_unmatched_1.fastq.gz',
                './folder/s1_unmatched_2.fastq.gz',
                None),
               ('s2', 'SKD8.640184',
                './folder/s2_paired_1.fastq.gz',
                './folder/s2_paired_2.fastq.gz',
                './folder/s2_unmatched_1.fastq.gz',
                './folder/s2_unmatched_2.fastq.gz',
                None),
               ('s3', 'SKB7.640196',
                './folder/s3_paired_1.fastq.gz',
                './folder/s3_paired_2.fastq.gz',
                './folder/s3_unmatched_1.fastq.gz',
                './folder/s3_unmatched_2.fastq.gz',
                None)]

        obs = make_read_sets_per_sample(files, fp)

        self.assertEqual(obs, exp)

    def test_make_read_sets_per_sample_match_fwd_rev_extra_f(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_p_fp = ['./folder/s3_paired_1.fastq.gz',
                    './folder/s2_paired_1.fastq.gz',
                    './folder/s1_paired_1.fastq.gz']

        rev_p_fp = ['./folder/s3_paired_2.fastq.gz',
                    './folder/s2_paired_2.fastq.gz',
                    './folder/s1_paired_2.fastq.gz']

        fwd_u_fp = ['./folder/s3_unmatched_1.fastq.gz',
                    './folder/s2_unmatched_1.fastq.gz',
                    './folder/s1_unmatched_1.fastq.gz']

        rev_u_fp = ['./folder/s3_unmatched_2.fastq.gz',
                    './folder/s2_unmatched_2.fastq.gz',
                    './folder/s1_unmatched_2.fastq.gz']

        files = (fwd_p_fp + rev_p_fp + fwd_u_fp + rev_u_fp +
                 ['./random/file.txt'])

        exp = [('s1', 'SKB8.640193',
                './folder/s1_paired_1.fastq.gz',
                './folder/s1_paired_2.fastq.gz',
                './folder/s1_unmatched_1.fastq.gz',
                './folder/s1_unmatched_2.fastq.gz',
                None),
               ('s2', 'SKD8.640184',
                './folder/s2_paired_1.fastq.gz',
                './folder/s2_paired_2.fastq.gz',
                './folder/s2_unmatched_1.fastq.gz',
                './folder/s2_unmatched_2.fastq.gz',
                None),
               ('s3', 'SKB7.640196',
                './folder/s3_paired_1.fastq.gz',
                './folder/s3_paired_2.fastq.gz',
                './folder/s3_unmatched_1.fastq.gz',
                './folder/s3_unmatched_2.fastq.gz',
                None)]

        obs = make_read_sets_per_sample(files, fp)

        self.assertEqual(obs, exp)

    def test_make_read_sets_per_sample_match_fwd_only(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_fp = ['./folder/s3_S013_L001.fastq.gz',
                  './folder/s2_S011_L001.fastq.gz',
                  './folder/s1_S009_L001.fastq.gz']

        files = fwd_fp

        exp = [('s1', 'SKB8.640193',
                None,
                None,
                None,
                None,
                './folder/s1_S009_L001.fastq.gz'),
               ('s2', 'SKD8.640184',
                None,
                None,
                None,
                None,
                './folder/s2_S011_L001.fastq.gz'),
               ('s3', 'SKB7.640196',
                None,
                None,
                None,
                None,
                './folder/s3_S013_L001.fastq.gz')]

        obs = make_read_sets_per_sample(files, fp)

        self.assertEqual(obs, exp)

    def test_make_read_sets_per_sample_extra_fwd(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_p_fp = ['./folder/s3_paired_1.fastq.gz',
                    './folder/s2_paired_1.fastq.gz',
                    './folder/s1_paired_1.fastq.gz',
                    './folder/s0_paired_1.fastq.gz']

        rev_p_fp = ['./folder/s3_paired_2.fastq.gz',
                    './folder/s2_paired_2.fastq.gz',
                    './folder/s1_paired_2.fastq.gz']

        fwd_u_fp = ['./folder/s3_unmatched_1.fastq.gz',
                    './folder/s2_unmatched_1.fastq.gz',
                    './folder/s1_unmatched_1.fastq.gz']

        rev_u_fp = ['./folder/s3_unmatched_2.fastq.gz',
                    './folder/s2_unmatched_2.fastq.gz',
                    './folder/s1_unmatched_2.fastq.gz']

        files = fwd_p_fp + rev_p_fp + fwd_u_fp + rev_u_fp

        with self.assertRaises(ValueError):
            make_read_sets_per_sample(files, fp)

    def test_make_read_sets_per_sample_fwd_single(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_p_fp = ['./folder/s3_paired_1.fastq.gz',
                    './folder/s2_paired_1.fastq.gz',
                    './folder/s1_paired_1.fastq.gz',
                    './folder/s0_paired_1.fastq.gz']

        rev_p_fp = ['./folder/s3_paired_2.fastq.gz',
                    './folder/s2_paired_2.fastq.gz',
                    './folder/s1_paired_2.fastq.gz']

        fwd_u_fp = ['./folder/s3_unmatched_1.fastq.gz',
                    './folder/s2_unmatched_1.fastq.gz',
                    './folder/s1_unmatched_1.fastq.gz']

        rev_u_fp = ['./folder/s3_unmatched_2.fastq.gz',
                    './folder/s2_unmatched_2.fastq.gz',
                    './folder/s1_unmatched_2.fastq.gz']

        single = ['./folder/s3_single.fasta.gz']

        files = fwd_p_fp + rev_p_fp + fwd_u_fp + rev_u_fp + single

        with self.assertRaises(ValueError):
            make_read_sets_per_sample(files, fp)

    def test_make_read_sets_per_sample_no_fastq(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        files = ['not_a_seqfile.txt']

        with self.assertRaises(ValueError):
            make_read_sets_per_sample(files, fp)

    def test_make_read_sets_per_sample_unmatch_fwd_rev(self):
        fd, fp = mkstemp()
        close(fd)
        with open(fp, 'w') as f:
            f.write(MAPPING_FILE)
        self._clean_up_files.append(fp)

        fwd_p_fp = ['./folder/wrong_paired_1.fastq.gz',
                    './folder/s2_paired_1.fastq.gz',
                    './folder/s1_paired_1.fastq.gz']

        rev_p_fp = ['./folder/s3_paired_2.fastq.gz',
                    './folder/s2_paired_2.fastq.gz',
                    './folder/s1_paired_2.fastq.gz']

        fwd_u_fp = ['./folder/s3_unmatched_1.fastq.gz',
                    './folder/s2_unmatched_1.fastq.gz',
                    './folder/s1_unmatched_1.fastq.gz']

        rev_u_fp = ['./folder/s3_unmatched_2.fastq.gz',
                    './folder/s2_unmatched_2.fastq.gz',
                    './folder/s1_unmatched_2.fastq.gz']

        files = fwd_p_fp + rev_p_fp + fwd_u_fp + rev_u_fp

        with self.assertRaises(ValueError):
            make_read_sets_per_sample(files, fp)

    def test_make_single_fastq_gz_fwd_rev(self):
        self.params['read-set'] = 'fwd_rev'

        read_sets = [('s1', 'SKB8.640193',
                      './support_files/s1_paired_1.fastq.gz',
                      './support_files/s1_paired_2.fastq.gz',
                      './support_files/s1_unmatched_1.fastq.gz',
                      './support_files/s1_unmatched_2.fastq.gz',
                      None),
                     ('s2', 'SKD8.640184',
                      './support_files/s2_paired_1.fastq.gz',
                      './support_files/s2_paired_2.fastq.gz',
                      './support_files/s2_unmatched_1.fastq.gz',
                      './support_files/s2_unmatched_2.fastq.gz',
                      None)]

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        exp_out = [('s1', 'SKB8.640193',
                    join(out_dir, 's1.fastq.gz')),
                   ('s2', 'SKD8.640184',
                    join(out_dir, 's2.fastq.gz'))]

        obs_out = make_single_fastq_gz(read_sets, out_dir, True)

        self.assertEqual(exp_out, obs_out)

        self.assertTrue(fcmp(obs_out[0][2], './support_files/s1.all.fastq.gz'))
        self.assertTrue(fcmp(obs_out[1][2], './support_files/s2.all.fastq.gz'))

    def test_make_single_fastq_gz_fwd_rev_empty(self):
        self.params['read-set'] = 'fwd_rev'

        read_sets = [('s1', 'SKB8.640193',
                      './support_files/s1_paired_1.fastq.gz',
                      './support_files/s1_paired_2.fastq.gz',
                      './support_files/s1_unmatched_1.fastq.gz',
                      './support_files/s1_unmatched_2.fastq.gz',
                      None),
                     ('s2', 'SKD8.640184',
                      './support_files/empty.fastq.gz',
                      './support_files/empty.fastq.gz',
                      './support_files/empty.fastq.gz',
                      './support_files/empty.fastq.gz',
                      None)]

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        exp_out = [('s1', 'SKB8.640193',
                    join(out_dir, 's1.fastq.gz'))]

        obs_out = make_single_fastq_gz(read_sets, out_dir, True)

        self.assertEqual(exp_out, obs_out)

        self.assertTrue(fcmp(obs_out[0][2],
                        './support_files/s1.all.fastq.gz'))

    def test_make_single_fastq_gz_paired_fwd(self):
        self.params['read-set'] = 'fwd'

        read_sets = [('s1', 'SKB8.640193',
                      './support_files/s1_paired_1.fastq.gz',
                      './support_files/s1_paired_2.fastq.gz',
                      './support_files/s1_unmatched_1.fastq.gz',
                      './support_files/s1_unmatched_2.fastq.gz',
                      None),
                     ('s2', 'SKD8.640184',
                      './support_files/s2_paired_1.fastq.gz',
                      './support_files/s2_paired_2.fastq.gz',
                      './support_files/s2_unmatched_1.fastq.gz',
                      './support_files/s2_unmatched_2.fastq.gz',
                      None)]

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        exp_out = [('s1', 'SKB8.640193',
                    join(out_dir, 's1.fastq.gz')),
                   ('s2', 'SKD8.640184',
                    join(out_dir, 's2.fastq.gz'))]

        obs_out = make_single_fastq_gz(read_sets, out_dir, False)

        self.assertEqual(exp_out, obs_out)

        self.assertTrue(fcmp(obs_out[0][2], './support_files/s1.fwd.fastq.gz'))
        self.assertTrue(fcmp(obs_out[1][2], './support_files/s2.fwd.fastq.gz'))

    def test_make_single_fastq_gz_single_fwd(self):
        self.params['read-set'] = 'fwd'

        read_sets = [('s1', 'SKB8.640193',
                      None,
                      None,
                      None,
                      None,
                      './support_files/s1_single.fastq.gz'),
                     ('s2', 'SKD8.640184',
                      None,
                      None,
                      None,
                      None,
                      './support_files/s2_single.fastq.gz')]

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        exp_out = [('s1', 'SKB8.640193',
                    join(out_dir, 's1.fastq.gz')),
                   ('s2', 'SKD8.640184',
                    join(out_dir, 's2.fastq.gz'))]

        obs_out = make_single_fastq_gz(read_sets, out_dir, False)

        self.assertEqual(exp_out, obs_out)

        self.assertTrue(fcmp(obs_out[0][2],
                        './support_files/s1_single.fastq.gz'))
        self.assertTrue(fcmp(obs_out[1][2],
                        './support_files/s2_single.fastq.gz'))

    def test_generate_humann2_analysis_commands(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        self.params['read-set'] = 'fwd'

        exp = [
            'humann2 --input "./folder/s1.fastq.gz" --output "%s/s1" '
            '--output-basename "SKB8.640193" --output-format biom '
            '--annotation-gene-index "8" '
            '--evalue "1.0" '
            '--gap-fill "off" '
            '--identity-threshold "50.0" '
            '--log-level "DEBUG" '
            '--memory-use "minimum" '
            '--metaphlan-options "-t rel_ab" '
            '--minpath "on" '
            '--output-format "biom" '
            '--output-max-decimals "10" '
            '--pathways "metacyc" '
            '--pick-frames "off" '
            '--prescreen-threshold "0.01" '
            '--remove-column-description-output '
            '--threads "1" '
            '--translated-alignment "diamond" '
            '--translated-query-coverage-threshold "90.0" '
            '--translated-subject-coverage-threshold "50.0" '
            '--xipe "off"' % out_dir,
            'humann2 --input "./folder/s2.fastq.gz" --output "%s/s2" '
            '--output-basename "SKD8.640184" --output-format biom '
            '--annotation-gene-index "8" '
            '--evalue "1.0" '
            '--gap-fill "off" '
            '--identity-threshold "50.0" '
            '--log-level "DEBUG" '
            '--memory-use "minimum" '
            '--metaphlan-options "-t rel_ab" '
            '--minpath "on" '
            '--output-format "biom" '
            '--output-max-decimals "10" '
            '--pathways "metacyc" '
            '--pick-frames "off" '
            '--prescreen-threshold "0.01" '
            '--remove-column-description-output '
            '--threads "1" '
            '--translated-alignment "diamond" '
            '--translated-query-coverage-threshold "90.0" '
            '--translated-subject-coverage-threshold "50.0" '
            '--xipe "off"' % out_dir,
            'humann2 --input "./folder/s3.fastq.gz" --output "%s/s3" '
            '--output-basename "SKB7.640196" --output-format biom '
            '--annotation-gene-index "8" '
            '--evalue "1.0" '
            '--gap-fill "off" '
            '--identity-threshold "50.0" '
            '--log-level "DEBUG" '
            '--memory-use "minimum" '
            '--metaphlan-options "-t rel_ab" '
            '--minpath "on" '
            '--output-format "biom" '
            '--output-max-decimals "10" '
            '--pathways "metacyc" '
            '--pick-frames "off" '
            '--prescreen-threshold "0.01" '
            '--remove-column-description-output '
            '--threads "1" '
            '--translated-alignment "diamond" '
            '--translated-query-coverage-threshold "90.0" '
            '--translated-subject-coverage-threshold "50.0" '
            '--xipe "off"' % out_dir]

        params_fwd = dict(self.params)

        params_fwd.pop('read-set')

        combined_reads = [('s1', 'SKB8.640193', './folder/s1.fastq.gz'),
                          ('s2', 'SKD8.640184', './folder/s2.fastq.gz'),
                          ('s3', 'SKB7.640196', './folder/s3.fastq.gz')]

        obs = generate_humann2_analysis_commands(combined_reads, out_dir,
                                                 params_fwd)

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
                'run_prefix': basename(fp2).replace('.fastq.gz', '')}}
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([
                (fp1, 'preprocessed_fastq'),
                (fp2, 'preprocessed_fastq')]),
            'type': "per_sample_FASTQ",
            'name': "New test artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid
        # overwriting so the run uses the defaults
        self.params['nucleotide-database'] = ''
        self.params['protein-database'] = ''
        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-shotgun', '0.0.1', 'HUMAnN2 0.9.1']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = humann2(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)
        # we are expecting 6 artifacts
        self.assertEqual(12, len(ainfo))

        obs_fps = []
        for a in ainfo:
            self.assertEqual("BIOM", a.artifact_type)
            obs_fps.append(a.files)
        exp_fps = [
            [(join(out_dir, 'genefamilies.biom'), 'biom')],
            [(join(out_dir, 'pathcoverage.biom'), 'biom')],
            [(join(out_dir, 'pathabundance.biom'), 'biom')],
            [(join(out_dir, 'genefamilies_cpm.biom'), 'biom')],
            [(join(out_dir, 'pathcoverage_relab.biom'), 'biom')],
            [(join(out_dir, 'pathabundance_relab.biom'), 'biom')],
            [(join(out_dir, 'genefamilies_cpm_stratified.biom'), 'biom')],
            [(join(out_dir, 'pathcoverage_relab_stratified.biom'), 'biom')],
            [(join(out_dir, 'pathabundance_relab_stratified.biom'), 'biom')],
            [(join(out_dir, 'genefamilies_cpm_unstratified.biom'), 'biom')],
            [(join(out_dir, 'pathcoverage_relab_unstratified.biom'), 'biom')],
            [(join(out_dir, 'pathabundance_relab_unstratified.biom'), 'biom')]]
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

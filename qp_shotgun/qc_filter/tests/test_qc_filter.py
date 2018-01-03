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

from qiita_client.testing import PluginTestCase

from qp_shotgun import plugin
from qp_shotgun.qc_filter.qc_filter import (qc_filter)
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

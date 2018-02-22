# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from qiita_client.testing import PluginTestCase

from os import remove
from os.path import exists, isdir
from shutil import rmtree

from qp_shotgun import plugin
# from qp_shotgun.shogun.shogun import shogun


class ShogunTests(PluginTestCase):

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

    def test_shogun(self):
        pass


if __name__ == '__main__':
    main()

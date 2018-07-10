# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .trim import trim_cmd
from .filter import filter_cmd
from .shogun import shogun_cmd


# Initialize the plugin
plugin = QiitaPlugin(
    'qp-shogun', '0.1.1', 'Shogun analysis tools for shotgun data')

plugin.register_command(trim_cmd)
plugin.register_command(filter_cmd)
plugin.register_command(shogun_cmd)

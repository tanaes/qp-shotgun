# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .qc_trim import qc_trim_cmd
from .qc_filter import qc_filter_cmd
from .shogun import shogun_cmd


# Initialize the plugin
plugin = QiitaPlugin(
    'qp-shogun', '0.0.1', 'Shogun analysis tools for shotgun data')

plugin.register_command(qc_trim_cmd)
plugin.register_command(qc_filter_cmd)
plugin.register_command(shogun_cmd)

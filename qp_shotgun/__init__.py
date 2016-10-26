# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin

from .humann2 import humann2_cmd

# Initialize the plugin
plugin = QiitaPlugin(
    'HUMAnN2', '0.9.1', 'HUMAnN2 is the next generation of HUMAnN (HMP '
    'Unified Metabolic Analysis Network)')

plugin.register_command(humann2_cmd)

#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import sys
from setuptools import setup
from stat import S_IEXEC
from os import (chdir, getcwd, chmod, rename, stat)
from os.path import join
from glob import glob
from tempfile import mkdtemp
from urllib import FancyURLopener
from subprocess import Popen, PIPE
from shutil import rmtree, copy

__version__ = "0.1.0-dev"


# The next block has been taken from, with minor editions:
# https://github.com/biocore/qiime/blob/master/setup.py#L44
#

class URLOpener(FancyURLopener):
    def http_error_default(self, url, fp, errcode, errmsg, headers):
        raise IOError(
            'Could not download %s\nPlease ensure the URL is valid and that '
            'you have an active Internet connection.' % url)


def system_call(cmd, error_msg):
    """Call `cmd` and return whether it was successful or not.
    This function is taken and modified from qcli (previously
    `qcli_system_call`).
    """
    proc = Popen(cmd,
                 shell=True,
                 universal_newlines=True,
                 stdout=PIPE,
                 stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    success = return_value == 0
    if not success:
        status("Unable to %s:" % error_msg)
        status("  stdout:")
        for line in stdout.split('\n'):
            status("    " + line)
        status("  stderr:")
        for line in stderr.split('\n'):
            status("    " + line)
    return success


def status(msg):
    """Write message immediately to stdout."""
    sys.stdout.write(msg)
    sys.stdout.write('\n')
    sys.stdout.flush()


def download_file(URL, dest_dir, local_file, num_retries=4):
    """General file downloader

    Inputs:
    URL: string to download the file from
    dest_dir: directory where you want to download the file
    local_file: output filename of the download
    num_retries: number of times the function will try to download the file

    Output:
    return_code: exit status for the download 0 = success, 1 = fail
    """
    status('  Downloading %s...' % local_file)

    url_opener = URLOpener()
    localFP = join(dest_dir, local_file)
    tmpDownloadFP = '%s.part' % localFP

    return_code = 1
    while num_retries > 0:
        try:
            tmpLocalFP, headers = url_opener.retrieve(URL, tmpDownloadFP)
            rename(tmpDownloadFP, localFP)
            return_code = 0
        except IOError:
            if num_retries == 1:
                status('  Download of %s failed.' % URL)
            else:
                status('  Download failed. Trying again... %d tries remain.' %
                       (num_retries - 1))
            num_retries -= 1
        else:
            num_retries = 0
            status('  %s downloaded successfully.' % local_file)
    return return_code


def download_metaphlan2():
    """Download the metaphlan2 executable and mv to scripts directory"""
    status("Installing metaphlan2.py ...")

    cwd = getcwd()
    scripts = join(cwd, 'scripts')

    tempdir = mkdtemp()
    URL = ('https://bitbucket.org/biobakery/metaphlan2/get/default.zip')
    if download_file(URL, tempdir, 'metaphlan2.py'):
        status("Could not download SortMeRNA, so cannot install it.\n")
        return

    chdir(tempdir)
    try:
        if not system_call('unzip biobakery-metaphlan2-*.zip'):
            return

        chdir('biobakery-metaphlan2-*')
        # make the file an executable file
        fname = 'metaphlan2.py'
        chmod(fname, stat(fname).st_mode | S_IEXEC)

        copy('db_v20', scripts)
        copy('metaphlan2.py', scripts)

        status("metaphlan2.py installed.\n")
    except:
        status("metaphlan2.py could not be installed.\n")
    finally:
        # remove the source
        rmtree(tempdir)
        chdir(cwd)


def catch_install_errors(install_function, name):
    try:
        install_function()
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        exception_type, exception_value = sys.exc_info()[:2]
        status(
            "Skipping installation of %s due to failure while downloading, "
            "building, or installing:\n  %s: %s\n" %
            (name, exception_type.__name__, exception_value))

catch_install_errors(download_metaphlan2, 'metaphlan2.py')

#
# end of block
#

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""


with open('README.rst') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='qp-shotgun',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Qiita Plugin: Shotgun',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/biocore/qiita',
      test_suite='nose.collector',
      packages=['qp_shotgun', 'qp_shotgun/humann2'],
      package_data={'qp_shotgun': ['support_files/config_file.cfg']},
      scripts=glob('scripts/*'),
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'future', 'pandas >= 0.15', 'humann2',
                        'h5py >= 2.3.1', 'biom-format'],
      classifiers=classifiers
      )

# -*- coding: utf-8 -*-

import os
import re
import configuration
from configparser import ConfigParser
import hgvs

CONF_ROOT = os.environ.get('CONF_ROOT')


# Config Section Mapping function
def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print ("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


# Configure
Config = ConfigParser()
Config.read(os.path.join(CONF_ROOT, 'config.ini'))

# Set SeqRepo and UTA environment variables
HGVS_SEQREPO_DIR = os.environ.get('HGVS_SEQREPO_DIR')
UTA_DB_URL = os.environ.get('UTA_DB_URL')
if HGVS_SEQREPO_DIR is None:
    HGVS_SEQREPO_DIR = ConfigSectionMap("SeqRepo")['seqrepo_dir']
    os.environ['HGVS_SEQREPO_DIR'] = HGVS_SEQREPO_DIR
if UTA_DB_URL is None:
    UTA_DB_URL = ConfigSectionMap("UTA")['uta_url']
    os.environ['UTA_DB_URL'] = UTA_DB_URL
__version__ = ConfigSectionMap("variantFormatter")['version']
__released__ = ConfigSectionMap("variantFormatter")['release_date']

# Collect metadata objects
hgvs_version = hgvs.__version__,
__hgvs_version__ = str(hgvs_version[0])
__seqrepo_db__ = str(HGVS_SEQREPO_DIR.split('/')[-1])
__uta_schema__ = str(UTA_DB_URL.split('/')[-1])


# <LICENSE>
# Copyright (C) 2019  Peter Causey-Freeman, University of Leicester
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
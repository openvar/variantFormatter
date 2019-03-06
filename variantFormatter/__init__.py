# -*- coding: utf-8 -*-

import os
import re
import configuration
from configparser import ConfigParser
import hgvs
import backports

CONF_ROOT = os.environ.get('HOME')

try:
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
    Config.read(os.path.join(CONF_ROOT, '.VariantFormatter.conf'))
    __version__ = ConfigSectionMap("VariantFormatter")['version']

except backports.configparser.NoSectionError:
    pass
except backports.configparser.ParsingError:
    pass




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
import os
import VariantFormatter.variantformatter as vf
vfo = vf.initializeFormatter()

ROOT = os.path.dirname(os.path.abspath(__file__))

class NotSetError(Exception):
    pass

inputfile = os.path.join(ROOT, 'inputVariants.txt')
outputfile = os.path.join(ROOT, 'test_inputs_auto.py')
out = open(outputfile, 'a')

# Set number by referring to test_inputs_auto.py
number = 0
tx_set = 'refseq'

if number is None:
    raise NotSetError('Add start number from test_inputs_auto.py')
if tx_set is None:
    raise NotSetError('Add requested transcript set')

out.write('import VariantFormatter\n'
          'import VariantFormatter.variantformatter as vf\n'
          'vfo = vf.initializeFormatter()\n'
          'class TestVariantsAuto(object):\n\n'
          '\t@classmethod\n'
          '\tdef setup_class(cls):\n'
          '\t\tVariantFormatter.__version__\n\n')

with open(inputfile) as ins:
    for l in ins:
        number += 1
        testline = '\n\tdef test_variant%s(self):\n' % number
        testline += "\t\tvariant = '%s'\n\t\tresults = vf.FormatVariant(variant, 'GRCh37', vfo,  '%s', None)\n\t\tresults = results.stucture_data()\n\t\tprint(results)\n\n" % (l.strip(), tx_set)
        result = vf.FormatVariant(l.strip(), 'GRCh37', vfo,  tx_set, None)
        res = result.stucture_data()
        print(("Variant %s: %s" % (number, l.strip())))
        for key in list(res.keys()):
            testline += "\t\tassert '%s' in results.keys()\n" % key
            for secondkey, item in list(res[key].items()):
                if secondkey == "hgvs_t_and_p":
                    if res[key]["hgvs_t_and_p"] is None:
                        testline += "\t\tassert results['%s']['%s'] is None\n" % (key, "hgvs_t_and_p")
                        continue
                    else:
                        for thirdkey, itm in list(res[key]["hgvs_t_and_p"].items()):
                            testline += "\t\tassert '%s' in results['%s']['%s'].keys()\n" % (thirdkey, key, "hgvs_t_and_p")
                            for forthkey, it_m in list(res[key]["hgvs_t_and_p"][thirdkey].items()):
                                if 'gap' in forthkey:
                                    continue
                                elif it_m is None:
                                    testline += "\t\tassert results['%s']['%s']['%s']['%s'] is None\n" % (
                                        key, 'hgvs_t_and_p', thirdkey, forthkey)
                                else:
                                    testline += "\t\tassert results['%s']['%s']['%s']['%s'] == '%s'\n" % (
                                    key, 'hgvs_t_and_p', thirdkey, forthkey, it_m)
                        continue
                elif item is None:
                    testline += "\t\tassert results['%s']['%s'] is None\n" % (key, secondkey)
                else:
                    testline += "\t\tassert results['%s']['%s'] == '%s'\n" % (key, secondkey, item)
            testline += '\n'
        out.write(testline)


out.close()

# <LICENSE>
# Copyright (C) 2019  Peter Causey-Freeman, University of Manchester
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
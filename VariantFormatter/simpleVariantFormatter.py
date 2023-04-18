# -*- coding: utf-8 -*-

"""
This module creates an initialization object.
This object connects to the hgvs Python library and associated databases

The Initialization object is used by FormatVariant
The FormatVariant object contains all HGVS descriptions available for a given genomic variant, g_to_p
"""

import re
import collections
import VariantValidator
import VariantFormatter
import VariantFormatter.variantformatter as vf
vfo = VariantValidator.Validator()

# Collect metadata
metadata = vfo.my_config()
print(metadata)
metadata['variantformatter_version'] = VariantFormatter.__version__
sr_root, sr_version = metadata['vvseqrepo_db'].split('/')[-2:]
metadata['vvseqrepo_db'] = '/'.join([sr_root, sr_version])


# If called in a threaded environment you MUST set validator to a thread local
# VariantValidator instance, due to non thread-safe SQLite3 access via SeqRepo
def format(batch_input, genome_build, transcript_model=None, specify_transcripts=None,
           checkOnly=False, liftover=False, validator=vfo, testing=False):

    # Testing?
    if testing is True:
        validator.testing = True

    # Set select_transcripts == 'all' to None
    vfo.select_transcripts = specify_transcripts
    if specify_transcripts == 'all':
        specify_transcripts = None
    is_a_list = type(batch_input) is list
    if is_a_list is True:
        batch_list = batch_input
    else:
        batch_list = batch_input.split('|')
    # batch_vars = []
    formatted_variants = collections.OrderedDict()
    for variant in batch_list:
        # remove external whitespace
        variant = variant.strip()
        # Remove internal whitespace
        wsl = variant.split()
        variant = ''.join(wsl)
        formatted_variants[variant] = collections.OrderedDict()
        formatted_variants[variant]['errors'] = []
        # Set validation warning flag
        formatted_variants[variant]['flag'] = None
        format_these = []
        if re.match('chr[\w\d]+-', variant) or re.match('chr[\w\d]+:', variant) or re.match('[\w\d]+-', variant)\
                or re.match('[\w\d]+:', variant):
            pseudo_vcf = variant
            if re.search(':', pseudo_vcf):
                vcf_list = pseudo_vcf.split(':')
                delimiter = ':'
            else:
                vcf_list = pseudo_vcf.split('-')
                delimiter = '-'
            if len(vcf_list) != 4:
                error = ('%s is an unsupported format: For assistance, submit variant description '
                         'to https://rest.variantvalidator.org' % str(pseudo_vcf))
                formatted_variants[variant]['errors'].append(error)
                formatted_variants[variant]['flag'] = 'submission_warning'
                continue
            if ',' in str(vcf_list[-1]):
                alts = vcf_list[-1].split(',')
                for eachalt in alts:
                    base = vcf_list[:3]
                    base.append(eachalt)
                    pv = delimiter.join(base)
                    format_these.append(pv)
            else:
                try:
                    format_these.append(variant)
                except Exception:
                    error = ('%s is an unsupported format: For assistance, submit variant description '
                             'to https://rest.variantvalidator.org' % variant)
                    formatted_variants[variant]['errors'].append(error)
                    formatted_variants[variant]['flag'] = 'submission_warning'
                    continue

        else:
            format_these.append(variant)

        # Processing
        for needs_formatting in format_these:
            result = vf.FormatVariant(needs_formatting, genome_build, validator,  transcript_model, specify_transcripts,
                                      checkOnly, liftover)
            res = result.stucture_data()
            formatted_variants[variant]['flag'] = result.warning_level
            formatted_variants[variant][needs_formatting] = res[needs_formatting]

    # Add metadata
    formatted_variants['metadata'] = metadata
    return formatted_variants


# <LICENSE>
# Copyright (C) 2016-2023 VariantValidator Contributors
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

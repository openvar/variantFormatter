# -*- coding: utf-8 -*-

"""
This module creates an initialization object.
This object connects to the hgvs Python library and associated databases

The Initialization object is used by FormatVariant
The FormatVariant object contains all HGVS descriptions available for a given genomic variant, g_to_p
"""

import re
import json
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

    # Format specify transcripts key options
    if specify_transcripts == '["all"]':
        specify_transcripts = "all"
    if specify_transcripts == '["raw"]':
        specify_transcripts = "raw"
    if specify_transcripts == '["mane"]':
        specify_transcripts = "mane"
    if specify_transcripts == '["mane_select"]':
        specify_transcripts = "mane_select"
    if specify_transcripts == '["select"]':
        specify_transcripts = "select"

    # Set select_transcripts == 'all' to None
    vfo.select_transcripts = specify_transcripts
    if specify_transcripts == 'all':
        specify_transcripts = None
    is_a_list = type(batch_input) is list
    if is_a_list is True:
        batch_list = batch_input
    else:
        try:
            batch_list = json.loads(batch_input)
        except json.decoder.JSONDecodeError:
            batch_list = [batch_input]

    # batch_vars = []
    formatted_variants = collections.OrderedDict()
    for variant in batch_list:
        # Set a bypass variable
        bypass = False
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
        # specially exclude LRGs as they do not have a "." separated
        # version number, we don't handle them (yet?) but they are not VCF
        if not variant.startswith('LRG') and (
                re.match('chr[\w\d]+-', variant) or
                re.match('chr[\w\d]+:', variant) or
                re.match('[\w\d]+-', variant) or
                re.match('[\w\d]+:', variant)):
            pseudo_vcf = variant

            if re.search(':', pseudo_vcf):
                vcf_list = pseudo_vcf.split(':')
                delimiter = ':'
            else:
                vcf_list = pseudo_vcf.split('-')
                delimiter = '-'
            if len(vcf_list) != 4:
                # Is it a hybrid format that VV can handle?
                try:
                    result = vfo.validate(variant, genome_build, "check_only").format_as_dict(test=True)
                except Exception:
                    pass
                try:
                    if "NC_" in result["intergenic_variant_1"]["primary_assembly_loci"][
                                       genome_build.lower()]["hgvs_genomic_description"]:

                        format_these.append(result["intergenic_variant_1"]["primary_assembly_loci"][
                                            genome_build.lower()]["hgvs_genomic_description"])
                        error = (f"{pseudo_vcf} is not HGVS compliant because a valid reference sequence has not been "
                                 f"provided. Updating to "
                                 f"{result['intergenic_variant_1']['primary_assembly_loci'][genome_build.lower()]['hgvs_genomic_description']}")
                        formatted_variants[variant]['errors'].append(error)
                        bypass = True
                    else:
                        error = ('%s is an unsupported format: For assistance, submit variant description '
                                 'to https://rest.variantvalidator.org' % str(pseudo_vcf))
                        formatted_variants[variant]['errors'].append(error)
                        formatted_variants[variant]['flag'] = 'submission_warning'
                        continue
                except KeyError:
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
                if bypass is True:
                    pass
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

        for needs_formatting in format_these:

            try:
                result = vf.FormatVariant(needs_formatting, genome_build, validator,  transcript_model,
                                          specify_transcripts, checkOnly, liftover)
            except Exception as e:
                import traceback
                traceback.print_exc()

            res = result.stucture_data()
            formatted_variants[variant]['flag'] = result.warning_level
            formatted_variants[variant][needs_formatting] = res[needs_formatting]

    # Add metadata
    formatted_variants['metadata'] = metadata
    return formatted_variants


# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
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

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
from VariantValidator.modules import vcf_to_pvcf


# ---------------------------------------------------------------------
# Lazy global state (legacy functional API only)
# ---------------------------------------------------------------------
_GLOBAL_VFO = None
_METADATA = None


def _get_global_validator():
    """
    Lazily create the legacy global VariantValidator instance.
    """
    global _GLOBAL_VFO
    if _GLOBAL_VFO is None:
        _GLOBAL_VFO = VariantValidator.Validator()
    return _GLOBAL_VFO


def _get_metadata():
    """
    Lazily compute metadata (requires a validator).
    """
    global _METADATA
    if _METADATA is None:
        vfo = _get_global_validator()
        _METADATA = vfo.my_config()
        _METADATA['variantformatter_version'] = VariantFormatter.__version__
        sr_root, sr_version = _METADATA['vvseqrepo_db'].split('/')[-2:]
        _METADATA['vvseqrepo_db'] = '/'.join([sr_root, sr_version])
    return _METADATA


# ---------------------------------------------------------------------
# Internal shared implementation
# ---------------------------------------------------------------------
def _format_impl(batch_input, genome_build, transcript_model=None,
                 specify_transcripts=None, checkOnly=False, liftover=False,
                 validator=None, testing=None):

    # Testing?
    validator.testing = bool(testing)

    # Normalise specify_transcripts
    mapping = {
        '["all"]': "all",
        '["raw"]': "raw",
        '["mane"]': "mane",
        '["mane_select"]': "mane_select",
        '["select"]': "select",
    }
    specify_transcripts = mapping.get(specify_transcripts, specify_transcripts)

    vfo = validator
    vfo.select_transcripts = specify_transcripts

    if specify_transcripts == 'all':
        specify_transcripts = None

    if isinstance(batch_input, list):
        batch_list = batch_input
    else:
        try:
            batch_list = json.loads(batch_input)
        except json.decoder.JSONDecodeError:
            batch_list = [batch_input]

    formatted_variants = collections.OrderedDict()

    for variant in batch_list:
        bypass = False
        variant = variant.strip()
        vcf_processing_warnings = []

        # VCF handling
        if "\t" in variant and not re.search(r"[gcrnmo]\.", variant):
            try:
                variant = vcf_to_pvcf.vcf_to_shorthand(variant)
            except vcf_to_pvcf.VcfConversionError:
                pass
            else:
                vcf_processing_warnings.append(
                    f"VcfConversionWarning: VCF line identified and converted to {variant}"
                )
                vcf_data = re.split(r'[-:]', variant)
                if (
                    re.search(r"\d+", vcf_data[2]) and
                    (re.search("del", vcf_data[3], re.IGNORECASE) or
                     re.search("inv", vcf_data[3], re.IGNORECASE))
                ):
                    if not re.search(r"[gatcnmo]\.", str(vcf_data)):
                        variant = (
                            f"{vcf_data[0]}:{vcf_data[1]}_"
                            f"{vcf_data[2]}{vcf_data[3].lower()}"
                        )
                        vcf_processing_warnings.append(
                            f"VcfConversionWarning: CNV identified, and mapped to {variant}"
                        )

        variant = ''.join(variant.split())
        formatted_variants[variant] = collections.OrderedDict()
        formatted_variants[variant]['errors'] = []
        formatted_variants[variant]['flag'] = None

        format_these = []

        if not variant.startswith('LRG') and (
            re.match(r'chr[\w\d]+[-:]', variant) or
            re.match(r'[\w\d]+[-:]', variant)
        ):
            pseudo_vcf = variant
            delimiter = ':' if ':' in pseudo_vcf else '-'
            vcf_list = pseudo_vcf.split(delimiter)

            if len(vcf_list) != 4:
                try:
                    result = vfo.validate(
                        variant, genome_build, "check_only"
                    ).format_as_dict(test=True)
                    hgvs = result["intergenic_variant_1"][
                        "primary_assembly_loci"
                    ][genome_build.lower()]["hgvs_genomic_description"]

                    if "NC_" in hgvs:
                        format_these.append(hgvs)
                        formatted_variants[variant]['errors'].append(
                            f"{pseudo_vcf} is not HGVS compliant because a valid "
                            f"reference sequence has not been provided. "
                            f"Updating to {hgvs}"
                        )
                        bypass = True
                    else:
                        raise KeyError
                except Exception:
                    formatted_variants[variant]['errors'].append(
                        f"{pseudo_vcf} is an unsupported format: "
                        "For assistance, submit variant description "
                        "to https://rest.variantvalidator.org"
                    )
                    formatted_variants[variant]['flag'] = 'submission_warning'
                    continue

            if ',' in vcf_list[-1]:
                for alt in vcf_list[-1].split(','):
                    format_these.append(
                        delimiter.join(vcf_list[:3] + [alt])
                    )
            elif not bypass:
                format_these.append(variant)
        else:
            format_these.append(variant)

        for needs_formatting in format_these:
            result = vf.FormatVariant(
                needs_formatting, genome_build, vfo,
                transcript_model, specify_transcripts,
                checkOnly, liftover
            )
            res = result.stucture_data()
            formatted_variants[variant]['flag'] = result.warning_level
            formatted_variants[variant][needs_formatting] = res[needs_formatting]

            if vcf_processing_warnings:
                formatted_variants[variant][needs_formatting][
                    'genomic_variant_warnings'
                ] = vcf_processing_warnings

    formatted_variants['metadata'] = _get_metadata()
    return formatted_variants


# ---------------------------------------------------------------------
# Legacy functional API (backwards compatible)
# ---------------------------------------------------------------------
def format(batch_input, genome_build, transcript_model=None,
           specify_transcripts=None, checkOnly=False,
           liftover=False, validator=None, testing=None):
    """
    Legacy functional API.
    Not thread-safe unless validator is thread-local.
    """
    if validator is None:
        validator = _get_global_validator()

    return _format_impl(
        batch_input, genome_build,
        transcript_model=transcript_model,
        specify_transcripts=specify_transcripts,
        checkOnly=checkOnly,
        liftover=liftover,
        validator=validator,
        testing=testing
    )


# ---------------------------------------------------------------------
# Object-oriented API (safe, poolable)
# ---------------------------------------------------------------------
class SimpleVariantFormatter:
    """
    Object-oriented formatter.

    Each instance owns its own VariantValidator.Validator.
    Safe for pooling and concurrent use (one request per instance).
    """

    def __init__(self, *, testing=False):
        self.validator = VariantValidator.Validator()
        self.testing = testing

    def format(self, batch_input, genome_build, transcript_model=None,
               specify_transcripts=None, checkOnly=False, liftover=False):
        return _format_impl(
            batch_input, genome_build,
            transcript_model=transcript_model,
            specify_transcripts=specify_transcripts,
            checkOnly=checkOnly,
            liftover=liftover,
            validator=self.validator,
            testing=self.testing
        )
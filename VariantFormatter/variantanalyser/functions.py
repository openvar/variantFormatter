# -*- coding: utf-8 -*-
"""
functions.py

Module containing VariantFormatter sub-functions. The majoirty of these functions require
hgvs Python package top-level functions or sub-functions contained in uta.py and
seqfetcher.py
"""

# IMPORT REQUIRED PYTHON MODULES
import re
import os
import copy
import warnings
import hgvs
import hgvs.exceptions
from hgvs.exceptions import HGVSDataNotAvailableError
import supportFunctions as links
from Bio.Seq import Seq



# Needs functions from VariantFormatter - directory above, unless in a single directory
try:
    from VariantFormatter import supportedChromosomeBuilds as supported_chromosome_builds
except ImportError:
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0, parentdir)
    import supportedChromosomeBuilds as supported_chromosome_builds
except AttributeError:
    parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.sys.path.insert(0, parentdir)
    import supportedChromosomeBuilds as supported_chromosome_builds


"""
Enhanced transcript to genome position mapping function using evm
Deals with mapping from transcript positions that do not exist in the genomic sequence
i.e. the stated position aligns to a genomic gap!
Trys to ensure that a genomic position is always returned even if the c. or n. transcript
will not map to the specified genome build primary assembly.
Deals with transcript mapping to several genomic assemblies
Order 
Map to a single NC_ for the specified genome build primary assembly
Map to a single NC_ for an alternate genome build primary assembly
Map to an NT_ from the specified genome build
Map to an NT_ from an alternative genome build
Map to an NW_ from the specified genome build
Map to an NW_ from an alternative genome buildRequires parsed c. or n. object
returns parsed hgvs g. object
"""


def myevm_t_to_g(hgvs_c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm, utilise_gap_code):

    # store the input
    stored_hgvs_c = copy.deepcopy(hgvs_c)
    expand_out = 'false'

    if utilise_gap_code is True and (
            hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type == 'delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins' or hgvs_c.posedit.edit.type == 'inv'):

        # if NM_ need the n. position
        if re.match('NM_', str(hgvs_c.ac)):
            hgvs_c = no_norm_evm.c_to_n(hgvs_c)

        # Check for intronic
        try:
            hn.normalize(hgvs_c)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('intronic variant', error):
                pass
            elif re.search('Length implied by coordinates must equal sequence deletion length', error) and re.match(
                    'NR_', hgvs_c.ac):
                hgvs_c.posedit.pos.end.base = hgvs_c.posedit.pos.start.base + len(hgvs_c.posedit.edit.ref) - 1

        # Check again before continuing
        if re.search('\d+\+', str(hgvs_c.posedit.pos)) or re.search('\d+\-', str(hgvs_c.posedit.pos)) or re.search(
                '\*\d+\+', str(hgvs_c.posedit.pos)) or re.search('\*\d+\-', str(hgvs_c.posedit.pos)):
            pass

        else:
            try:
                # For non-intronic sequence
                hgvs_t = copy.deepcopy(hgvs_c)
                if hgvs_t.posedit.edit.type == 'inv':
                    inv_alt = revcomp(hgvs_t.posedit.edit.ref)
                    t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                        hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t_delins = hp.parse_hgvs_variant(t_delins)
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    inv_alt = pre_base + inv_alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(
                        end) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                elif hgvs_c.posedit.edit.type == 'dup':
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base
                    ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        (hgvs_t.posedit.pos.start.base + len(ref)) - 2) + 'del' + ref + 'ins' + alt
                    hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
                elif hgvs_c.posedit.edit.type == 'ins':
                    ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                           hgvs_t.posedit.pos.end.base + 1)
                    ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(
                        hgvs_t.posedit.pos.start.base - 1) + '_' + str(
                        hgvs_t.posedit.pos.end.base + 1) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                else:
                    if str(hgvs_t.posedit.edit.alt) == 'None':
                        hgvs_t.posedit.edit.alt = ''
                    pre_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 2,
                                            hgvs_t.posedit.pos.start.base - 1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.end.base,
                                             hgvs_t.posedit.pos.end.base + 1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(
                        hgvs_t.posedit.edit)
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)
                hgvs_c = copy.deepcopy(hgvs_t)

                # Set expanded out test to true
                expand_out = 'true'

            except Exception:
                hgvs_c = hgvs_c

        if re.match('NM_', str(hgvs_c.ac)):
            try:
                hgvs_c = no_norm_evm.n_to_c(hgvs_c)
            except hgvs.exceptions.HGVSError as e:
                hgvs_c = copy.deepcopy(stored_hgvs_c)

        # Ensure the altered c. variant has not crossed intro exon boundaries
        hgvs_check_boundaries = copy.deepcopy(hgvs_c)
        try:
            h_variant = hn.normalize(hgvs_check_boundaries)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if re.search('spanning the exon-intron boundary', error):
                hgvs_c = copy.deepcopy(stored_hgvs_c)
        # Catch identity at the exon/intron boundary by trying to normalize ref only
        if hgvs_check_boundaries.posedit.edit.type == 'identity':
            reform_ident = str(hgvs_c).split(':')[0]
            reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(
                hgvs_c.posedit.edit.ref)  # + 'ins' + str(hgvs_c.posedit.edit.alt)
            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
            try:
                hn.normalize(hgvs_reform_ident)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if re.search('spanning the exon-intron boundary', error) or re.search(
                        'Normalization of intronic variants', error):
                    hgvs_c = copy.deepcopy(stored_hgvs_c)
    try:
        hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
        hn.normalize(hgvs_genomic)  # Check the validity of the mapping
        # This will fail on multiple refs for NC_
    except hgvs.exceptions.HGVSError as e:
        # Recover all available mapping options from UTA
        mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac)

        if mapping_options == []:
            raise HGVSDataNotAvailableError(
                "No alignment data between the specified transcript reference sequence and any GRCh37 and GRCh38 genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) are available.")

        # Capture errors from attempted mappings
        attempted_mapping_error = ''

        for option in mapping_options:
            if re.match('blat', option[2]):
                continue
            if re.match('NC_', option[1]):
                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                if chr_num != 'false':
                    try:
                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                        break
                    except Exception as e:
                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                            1] + '~'
                        print e
                        continue

        # If not mapped, raise error
        try:
            hn.normalize(hgvs_genomic)
        except:
            for option in mapping_options:
                if re.match('blat', option[2]):
                    continue
                if re.match('NC_', option[1]):
                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                    if chr_num == 'false':
                        try:
                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                            break
                        except Exception as e:
                            if re.search(option[1], attempted_mapping_error):
                                pass
                            else:
                                attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                                    1] + '~'
                            print e
                            continue
            try:
                hn.normalize(hgvs_genomic)
            except:
                for option in mapping_options:
                    if re.match('blat', option[2]):
                        continue
                    if re.match('NT_', option[1]):
                        chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]), primary_assembly)
                        if chr_num != 'false':
                            try:
                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[
                                    1] + '~'
                                print e
                                continue
                try:
                    hn.normalize(hgvs_genomic)
                except:
                    for option in mapping_options:
                        if re.match('blat', option[2]):
                            continue
                        if re.match('NT_', option[1]):
                            chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                        primary_assembly)
                            if chr_num == 'false':
                                try:
                                    hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                    break
                                except Exception as e:
                                    attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                              option[
                                                                  1] + '~'
                                    print e
                                    continue
                    try:
                        hn.normalize(hgvs_genomic)
                    except:
                        for option in mapping_options:
                            if re.match('blat', option[2]):
                                continue
                            if re.match('NW_', option[1]):
                                chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                            primary_assembly)
                                if chr_num != 'false':
                                    try:
                                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                        break
                                    except Exception as e:
                                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                                  option[1] + '~'
                                        print e
                                        continue
                        try:
                            hn.normalize(hgvs_genomic)
                        except:
                            for option in mapping_options:
                                if re.match('blat', option[2]):
                                    continue
                                if re.match('NW_', option[1]):
                                    chr_num = supported_chromosome_builds.supported_for_mapping(str(option[1]),
                                                                                                primary_assembly)
                                    if chr_num == 'false':
                                        try:
                                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                            break
                                        except Exception as e:
                                            attempted_mapping_error = attempted_mapping_error + str(
                                                e) + "/" + hgvs_c.ac + "/" + \
                                                                      option[1] + '~'
                                            print e
                                            continue

                            # Only a RefSeqGene available
                            try:
                                hn.normalize(hgvs_genomic)
                            except:
                                for option in mapping_options:
                                    if re.match('blat', option[2]):
                                        continue
                                    if re.match('NG_', option[1]):
                                        try:
                                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                            break
                                        except Exception as e:
                                            attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + \
                                                                      option[1] + '~'
                                            print e
                                            continue

    # If not mapped, raise error
    try:
        hgvs_genomic
    except Exception:
        raise HGVSDataNotAvailableError(attempted_mapping_error)

    if hgvs_c.posedit.edit.type == 'identity' and hgvs_genomic.posedit.edit.type == 'delins' and hgvs_genomic.posedit.edit.alt == '' and expand_out != 'true':
        hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.ref
    if hgvs_genomic.posedit.edit.type == 'ins' and utilise_gap_code is True:
        try:
            hgvs_genomic = hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                   hgvs_genomic.posedit.pos.end.base)
                hgvs_genomic.posedit.edit.ref = ref
                hgvs_genomic.posedit.edit.alt = ref[0:1] + hgvs_genomic.posedit.edit.alt + ref[-1:]
                hgvs_genomic = hn.normalize(hgvs_genomic)
            if error == 'base start position must be <= end position':
                start = hgvs_genomic.posedit.pos.start.base
                end = hgvs_genomic.posedit.pos.end.base
                hgvs_genomic.posedit.pos.start.base = end
                hgvs_genomic.posedit.pos.end.base = start
                hgvs_genomic = hn.normalize(hgvs_genomic)

    # Statements required to reformat the stored_hgvs_c into a useable synonym
    if (stored_hgvs_c.posedit.edit.ref == '' or stored_hgvs_c.posedit.edit.ref is None) and expand_out != 'false':
        if stored_hgvs_c.type == 'c':
            stored_hgvs_n = vm.c_to_n(stored_hgvs_c)
        else:
            stored_hgvs_n = stored_hgvs_c
        stored_ref = sf.fetch_seq(str(stored_hgvs_n.ac), stored_hgvs_n.posedit.pos.start.base - 1,
                                  stored_hgvs_n.posedit.pos.end.base)
        stored_hgvs_c.posedit.edit.ref = stored_ref

    if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
        if hgvs_genomic.posedit.edit.type == 'ins':
            stored_ref = sf.fetch_seq(str(hgvs_genomic.ac), hgvs_genomic.posedit.pos.start.base - 1,
                                      hgvs_genomic.posedit.pos.end.base)
            stored_alt = stored_ref[:1] + hgvs_genomic.posedit.edit.alt + stored_ref[-1:]
            hgvs_genomic.posedit.edit.ref = stored_ref
            hgvs_genomic.posedit.edit.alt = stored_alt

    # First look for variants mapping to the flanks of gaps
    # either in the gap or on the flank but not fully within the gap
    if expand_out == 'true':

        nr_genomic = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

        try:
            hn.normalize(nr_genomic)
        except hgvs.exceptions.HGVSInvalidVariantError as e:
            error_type_1 = str(e)
            if re.match('Length implied by coordinates must equal sequence deletion length', str(e)) or str(
                    e) == 'base start position must be <= end position':
                # Effectively, this code is designed to handle variants that are directly proximal to
                # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed bases due to
                # the deletion length being > the specified range.

                # Warn of variant location wrt the gap
                if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                    try:
                        hn.normalize(genomic_gap_variant)
                    # Still a problem
                    except hgvs.exceptions.HGVSInvalidVariantError as e:
                        if 'base start position must be <= end position' in str(e) and \
                                'Length implied by coordinates must equal' in error_type_1:
                            make_gen_var = copy.copy(nr_genomic)
                            make_gen_var.posedit.edit.ref = sf.fetch_seq(nr_genomic.ac,
                                                                         nr_genomic.posedit.pos.start.base - 1,
                                                                         nr_genomic.posedit.pos.end.base)
                            genomic_gap_variant = make_gen_var

                            error_type_1 = None
                    else:
                        genomic_gap_variant = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)

                if error_type_1 == 'base start position must be <= end position':
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)

                # Logic
                # We have checked that the variant does not cross boundaries, or is intronic
                # So is likely mapping to a genomic gap
                try:
                    hn.normalize(genomic_gap_variant)
                except Exception as e:
                    if str(e) == 'base start position must be <= end position':
                        # This will only happen when the variant is fully within the gap
                        gap_start = genomic_gap_variant.posedit.pos.end.base
                        gap_end = genomic_gap_variant.posedit.pos.start.base
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                    if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        # This will only happen if the variant is flanking the gap but is
                        # not inside the gap
                        gap_start = genomic_gap_variant.posedit.pos.start.base - 1
                        gap_end = genomic_gap_variant.posedit.pos.end.base + 1
                        genomic_gap_variant.posedit.pos.start.base = gap_start
                        genomic_gap_variant.posedit.pos.end.base = gap_end
                        genomic_gap_variant.posedit.edit.ref = ''
                        stored_hgvs_c = copy.deepcopy(hgvs_c)

                    # Remove alt
                    try:
                        genomic_gap_variant.posedit.edit.alt = ''
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            pass

                    # Should be a delins so will normalize statically and replace the reference bases
                    genomic_gap_variant = hn.normalize(genomic_gap_variant)
                    # Static map to c. and static normalize
                    transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                    stored_transcript_gap_variant = transcript_gap_variant

                    if not re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                        try:
                            transcript_gap_variant = hn.normalize(transcript_gap_variant)
                        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                            if ' Unsupported normalization of variants spanning the UTR-exon boundary' in str(e):
                                pass

                    # if NM_ need the n. position
                    if re.match('NM_', str(hgvs_c.ac)):
                        transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                        transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                    else:
                        transcript_gap_n = transcript_gap_variant
                        transcript_gap_alt_n = stored_hgvs_c

                    # Ensure an ALT exists
                    try:
                        if transcript_gap_alt_n.posedit.edit.alt is None:
                            transcript_gap_alt_n.posedit.edit.alt = 'X'
                    except Exception as e:
                        if str(e) == "'Dup' object has no attribute 'alt'":
                            transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                                transcript_gap_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                            transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                                transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                                transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                            transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                    # Split the reference and replacing alt sequence into a dictionary
                    reference_bases = list(transcript_gap_n.posedit.edit.ref)
                    if transcript_gap_alt_n.posedit.edit.alt is not None:
                        alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                    else:
                        # Deletions with no ins
                        pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                        alternate_bases = []
                        for base in pre_alternate_bases:
                            alternate_bases.append('X')

                    # Create the dictionaries
                    ref_start = transcript_gap_n.posedit.pos.start.base
                    alt_start = transcript_gap_alt_n.posedit.pos.start.base
                    ref_base_dict = {}
                    for base in reference_bases:
                        ref_base_dict[ref_start] = str(base)
                        ref_start = ref_start + 1

                    alt_base_dict = {}

                    # Note, all variants will be forced into the format delete insert
                    # Deleted bases in the ALT will be substituted for X
                    for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                     transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                        if int == alt_start:
                            alt_base_dict[int] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[int] = 'X'

                            # Generate the alt sequence
                    alternate_sequence_bases = []
                    for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1,
                                     1):
                        if int in alt_base_dict.keys():
                            alternate_sequence_bases.append(alt_base_dict[int])
                        else:
                            alternate_sequence_bases.append(ref_base_dict[int])
                    alternate_sequence = ''.join(alternate_sequence_bases)
                    alternate_sequence = alternate_sequence.replace('X', '')

                    # Update variant, map to genome using vm and normalize
                    transcript_gap_n.posedit.edit.alt = alternate_sequence

                    try:
                        transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                    except:
                        transcript_gap_variant = transcript_gap_n

                    try:
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)
                    except Exception as e:
                        if str(e) == "base start position must be <= end position":
                            # Expansion out is required to map back to the genomic position
                            pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                    transcript_gap_n.posedit.pos.start.base - 1)
                            post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                     transcript_gap_n.posedit.pos.end.base + 1)
                            transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                            transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                            transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                            transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                            try:
                                transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                            except:
                                transcript_gap_variant = transcript_gap_n
                            hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                            hgvs_genomic = hn.normalize(hgvs_genomic)

                    # Bypass the next bit of gap code
                    expand_out = 'false'

            else:
                pass
        # No map to the flank of a gap or within the gap
        else:
            pass

    # CASCADING STATEMENTS WHICH CAPTURE t to g MAPPING OPTIONS
    # Remove identity bases
    if hgvs_c == stored_hgvs_c:
        expand_out = 'false'
    elif expand_out == 'false' or utilise_gap_code is False:
        pass
    # Correct expansion ref + 2
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) == (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        hgvs_genomic.posedit.pos.start.base = hgvs_genomic.posedit.pos.start.base + 1
        hgvs_genomic.posedit.pos.end.base = hgvs_genomic.posedit.pos.end.base - 1
        hgvs_genomic.posedit.edit.ref = hgvs_genomic.posedit.edit.ref[1:-1]
        if hgvs_genomic.posedit.edit.alt is not None:
            hgvs_genomic.posedit.edit.alt = hgvs_genomic.posedit.edit.alt[1:-1]
    elif expand_out == 'true' and (
            len(hgvs_genomic.posedit.edit.ref) != (len(stored_hgvs_c.posedit.edit.ref) + 2)):  # >= 3:
        if expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) == 2:
            gn = hn.normalize(hgvs_genomic)
            pass

        # Likely if the start or end position aligns to a gap in the genomic sequence
        # Logic
        # We have checked that the variant does not cross boundaries, or is intronic
        # So is likely mapping to a genomic gap
        elif expand_out == 'true' and len(hgvs_genomic.posedit.edit.ref) <= 1:
            # Incorrect expansion, likely < ref + 2
            genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
            try:
                hn.normalize(genomic_gap_variant)
            except Exception as e:
                if str(e) == 'base start position must be <= end position':
                    gap_start = genomic_gap_variant.posedit.pos.end.base
                    gap_end = genomic_gap_variant.posedit.pos.start.base
                    genomic_gap_variant.posedit.pos.start.base = gap_start
                    genomic_gap_variant.posedit.pos.end.base = gap_end
                # Remove alt
                try:
                    genomic_gap_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        pass
                # Should be a delins so will normalize statically and replace the reference bases
                genomic_gap_variant = hn.normalize(genomic_gap_variant)
                # Static map to c. and static normalize
                transcript_gap_variant = vm.g_to_t(genomic_gap_variant, hgvs_c.ac)
                stored_transcript_gap_variant = transcript_gap_variant
                transcript_gap_variant = hn.normalize(transcript_gap_variant)
                # if NM_ need the n. position
                if re.match('NM_', str(hgvs_c.ac)):
                    transcript_gap_n = no_norm_evm.c_to_n(transcript_gap_variant)
                    transcript_gap_alt_n = no_norm_evm.c_to_n(stored_hgvs_c)
                else:
                    transcript_gap_n = transcript_gap_variant
                    transcript_gap_alt_n = stored_hgvs_c

                # Ensure an ALT exists
                try:
                    if transcript_gap_alt_n.posedit.edit.alt is None:
                        transcript_gap_alt_n.posedit.edit.alt = 'X'
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(
                            transcript_gap_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' + transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref
                        transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)
                        transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(
                            transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(
                            transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' + transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref
                        transcript_gap_alt_n = hp.parse_hgvs_variant(transcript_gap_alt_n_delins_from_dup)

                # Split the reference and replacing alt sequence into a dictionary
                reference_bases = list(transcript_gap_n.posedit.edit.ref)
                if transcript_gap_alt_n.posedit.edit.alt is not None:
                    alternate_bases = list(transcript_gap_alt_n.posedit.edit.alt)
                else:
                    # Deletions with no ins
                    pre_alternate_bases = list(transcript_gap_alt_n.posedit.edit.ref)
                    alternate_bases = []
                    for base in pre_alternate_bases:
                        alternate_bases.append('X')

                # Create the dictionaries
                ref_start = transcript_gap_n.posedit.pos.start.base
                alt_start = transcript_gap_alt_n.posedit.pos.start.base
                ref_base_dict = {}
                for base in reference_bases:
                    ref_base_dict[ref_start] = str(base)
                    ref_start = ref_start + 1

                alt_base_dict = {}

                # Note, all variants will be forced into the format delete insert
                # Deleted bases in the ALT will be substituted for X
                for int in range(transcript_gap_alt_n.posedit.pos.start.base,
                                 transcript_gap_alt_n.posedit.pos.end.base + 1, 1):
                    if int == alt_start:
                        alt_base_dict[int] = str(''.join(alternate_bases))
                    else:
                        alt_base_dict[int] = 'X'

                # Generate the alt sequence
                alternate_sequence_bases = []
                for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base + 1, 1):
                    if int in alt_base_dict.keys():
                        alternate_sequence_bases.append(alt_base_dict[int])
                    else:
                        alternate_sequence_bases.append(ref_base_dict[int])
                alternate_sequence = ''.join(alternate_sequence_bases)
                alternate_sequence = alternate_sequence.replace('X', '')

                # Update variant, map to genome using vm and normalize
                transcript_gap_n.posedit.edit.alt = alternate_sequence

                try:
                    transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                except:
                    transcript_gap_variant = transcript_gap_n

                try:
                    hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                    hgvs_genomic = hn.normalize(hgvs_genomic)
                except Exception as e:
                    if str(e) == "base start position must be <= end position":
                        # Expansion out is required to map back to the genomic position
                        pre_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.start.base - 2,
                                                transcript_gap_n.posedit.pos.start.base - 1)
                        post_base = sf.fetch_seq(transcript_gap_n.ac, transcript_gap_n.posedit.pos.end.base,
                                                 transcript_gap_n.posedit.pos.end.base + 1)
                        transcript_gap_n.posedit.pos.start.base = transcript_gap_n.posedit.pos.start.base - 1
                        transcript_gap_n.posedit.pos.end.base = transcript_gap_n.posedit.pos.end.base + 1
                        transcript_gap_n.posedit.edit.ref = pre_base + transcript_gap_n.posedit.edit.ref + post_base
                        transcript_gap_n.posedit.edit.alt = pre_base + transcript_gap_n.posedit.edit.alt + post_base
                        try:
                            transcript_gap_variant = vm.n_to_c(transcript_gap_n)
                        except:
                            transcript_gap_variant = transcript_gap_n
                        hgvs_genomic = vm.t_to_g(transcript_gap_variant, hgvs_genomic.ac)
                        hgvs_genomic = hn.normalize(hgvs_genomic)

    # Ins variants map badly - Especially between c. exon/exon boundary
    if hgvs_c.posedit.edit.type == 'ins' and hgvs_c.posedit.pos.start.offset == 0 and hgvs_c.posedit.pos.end.offset == 0:
        try:
            hn.normalize(hgvs_genomic)
        except hgvs.exceptions.HGVSError as e:
            error = str(e)
            if error == 'insertion length must be 1':
                if hgvs_c.type == 'c':
                    hgvs_t = vm.c_to_n(hgvs_c)
                else:
                    hgvs_t = copy.copy(hgvs_c)
                ins_ref = sf.fetch_seq(str(hgvs_t.ac), hgvs_t.posedit.pos.start.base - 1, hgvs_t.posedit.pos.end.base)
                ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(
                    hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                try:
                    hgvs_c = vm.n_to_c(hgvs_t)
                except Exception:
                    hgvs_c = copy.copy(hgvs_t)
                try:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                except Exception as e:
                    error = str(e)

    return hgvs_genomic

"""
Function which takes a NORMALIZED hgvs Python transcript variant and maps to a specified protein reference sequence. A protein
level hgvs python object is returned.

Note the function currently assumes that the transcript description is correctly normalized having come from the 
previous g_to_t function
"""


def myc_to_p(hgvs_transcript, evm, hdp, hp, hn, vm, sf, re_to_p):
    # Create dictionary to store the information
    hgvs_transcript_to_hgvs_protein = {'error': '', 'hgvs_protein': '', 'ref_residues': ''}

    # Collect the associated protein
    if hgvs_transcript.type == 'c':
        associated_protein_accession = hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)
        # This method sometimes fails
        if str(associated_protein_accession) == 'None':
            cod = str(hgvs_transcript)
            cod = cod.replace('inv', 'del')
            cod = hp.parse_hgvs_variant(cod)
            p = evm.c_to_p(cod)
            associated_protein_accession = p.ac
    else:
        pass

        # Check for non-coding transcripts
    if hgvs_transcript.type == 'c':
        # Handle non inversions with simple c_to_p mapping

        if (hgvs_transcript.posedit.edit.type != 'inv') and (hgvs_transcript.posedit.edit.type != 'delins') and (
                re_to_p is False):
            # Does the edit affect the start codon?
            if ((
                        hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0) or (
                        hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                    and not re.search('\*', str(
                hgvs_transcript.posedit.pos)):
                residue_one = sf.fetch_seq(associated_protein_accession, start_i=1-1,end_i=1)
                threed_residue_one = links.one_to_three(residue_one)
                r_one_report = '(%s1?)' % threed_residue_one
                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                    type='p', posedit=r_one_report)
            else:
                try:
                    hgvs_protein = evm.c_to_p(hgvs_transcript)
                except IndexError as e:
                    error = str(e)
                    if re.search('string index out of range', error) and re.search('dup', str(hgvs_transcript)):
                        hgvs_ins = hp.parse_hgvs_variant(str(hgvs_transcript))
                        hgvs_ins = hn.normalize(hgvs_ins)
                        inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + str(
                            hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
                        hgvs_transcript = hp.parse_hgvs_variant(inst)
                        hgvs_protein = evm.c_to_p(hgvs_transcript)

            try:
                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                return hgvs_transcript_to_hgvs_protein
            except UnboundLocalError:
                hgvs_transcript_to_hgvs_protein = myc_to_p(hgvs_transcript, evm, hdp, hp, hn, vm, sf, re_to_p=True)
                return hgvs_transcript_to_hgvs_protein

        else:
            # Additional code required to process inversions
            # Note, this code was developed for VariantValidator and is not native to the biocommons hgvs Python package
            # Convert positions to n. position
            hgvs_naughty = vm.c_to_n(hgvs_transcript)

            # Collect the deleted sequence using fetch_seq
            del_seq = sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base - 1,
                                   end_i=hgvs_naughty.posedit.pos.end.base)

            # Make the inverted sequence
            my_seq = Seq(del_seq)

            if hgvs_transcript.posedit.edit.type == 'inv':
                inv_seq = my_seq.reverse_complement()
            else:
                inv_seq = hgvs_transcript.posedit.edit.alt
                if inv_seq is None:
                    inv_seq = ''

            # Look for p. delins or del
            not_delins = True
            if hgvs_transcript.posedit.edit.type != 'inv':
                try:
                    shifts = evm.c_to_p(hgvs_transcript)
                    if re.search('del', shifts.posedit.edit.type):
                        not_delins = False
                except Exception:
                    not_delins = False
            else:
                not_delins = False

            # Use inv delins code?
            if not_delins == False:
                # Collect the associated protein
                associated_protein_accession = hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)

                # Intronic inversions are marked as uncertain i.e. p.?
                if re.search('\d+\-', str(hgvs_transcript.posedit.pos)) or re.search('\d+\+', str(
                        hgvs_transcript.posedit.pos)) or re.search('\*', str(hgvs_transcript.posedit.pos)) or re.search(
                        '[cn].\-', str(hgvs_transcript)):
                    if ((
                                hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                        or
                        (
                                hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                            and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                        residue_one = sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                        threed_residue_one = links.one_to_three(residue_one)
                        r_one_report = '(%s1?)' % threed_residue_one # was (MET1?)
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                            type='p', posedit=r_one_report)
                    else:
                        # Make the variant
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                            posedit='?')
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                    return hgvs_transcript_to_hgvs_protein
                else:
                    # Need to obtain the cds_start
                    inf = hdp.get_tx_identity_info(hgvs_transcript.ac)
                    cds_start = inf[3]

                    # Extract the reference coding sequence from SeqRepo
                    try:
                        ref_seq = sf.fetch_seq(str(hgvs_naughty.ac))
                    except Exception as e:
                        error = str(e)
                        hgvs_transcript_to_hgvs_protein['error'] = error
                        return hgvs_transcript_to_hgvs_protein

                    # Create the variant coding sequence
                    var_seq = links.n_inversion(ref_seq, del_seq, inv_seq,
                                                hgvs_naughty.posedit.pos.start.base,
                                                hgvs_naughty.posedit.pos.end.base)
                    # Translate the reference and variant proteins
                    prot_ref_seq = links.translate(ref_seq, cds_start)

                    try:
                        prot_var_seq = links.translate(var_seq, cds_start)
                    except IndexError:
                        hgvs_transcript_to_hgvs_protein[
                            'error'] = 'Cannot identify an in-frame Termination codon in the variant mRNA sequence'
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                            posedit='?')
                        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                        return hgvs_transcript_to_hgvs_protein

                    if prot_ref_seq == 'error':
                        error = 'Unable to generate protein variant description'
                        hgvs_transcript_to_hgvs_protein['error'] = error
                        return hgvs_transcript_to_hgvs_protein
                    elif prot_var_seq == 'error':
                        # Does the edit affect the start codon?
                        if ((
                                    hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                            or
                            (
                                    hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                                and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                            residue_one = sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                            threed_residue_one # was (MET1?) = links.one_to_three(residue_one)
                            r_one_report = '(%s1?)' % threed_residue_one # was (MET1?)
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p', posedit=r_one_report)

                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein
                        else:
                            error = 'Unable to generate protein variant description'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein
                    else:
                        # Gather the required information regarding variant interval and sequences
                        if hgvs_transcript.posedit.edit.type != 'delins':
                            pro_inv_info = links.pro_inv_info(prot_ref_seq, prot_var_seq)
                        else:
                            pro_inv_info = links.pro_delins_info(prot_ref_seq, prot_var_seq)

                        # Error has occurred
                        if pro_inv_info['error'] == 'true':
                            error = 'Translation error occurred, please contact admin'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein

                        # The Nucleotide variant has not affected the protein sequence i.e. synonymous
                        elif pro_inv_info['variant'] != 'true':
                            # Make the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p',
                                                                                posedit='=')
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein

                        else:
                            # Early termination i.e. stop gained
                            # if pro_inv_info['terminate'] == 'true':
                            #     end = 'Ter' + str(pro_inv_info['ter_pos'])
                            #     pro_inv_info['prot_ins_seq'].replace('*', end)

                            # Complete variant description
                            # Recode the single letter del and ins sequences into three letter amino acid codes
                            del_thr = links.one_to_three(pro_inv_info['prot_del_seq'])
                            ins_thr = links.one_to_three(pro_inv_info['prot_ins_seq'])

                            # Write the HGVS position and edit
                            del_len = len(del_thr)
                            from_aa = del_thr[0:3]
                            to_aa = del_thr[del_len - 3:]

                            # Handle a range of amino acids
                            if pro_inv_info['edit_start'] != pro_inv_info['edit_end']:
                                if len(ins_thr) > 0:
                                    if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'delins' + ins_thr + '?)'
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'
                                else:
                                    if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'del?)'
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                            pro_inv_info['edit_end']) + 'del)'
                            else:
                                # Handle extended proteins i.e. stop_lost
                                if del_thr == 'Ter' and (len(ins_thr) > len(del_thr)):
                                    # Nucleotide variant range aligns to the Termination codon
                                    if ins_thr[-3:] == 'Ter':
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                            ins_thr[:3]) + 'ext' + str(ins_thr[-3:]) + str((len(ins_thr) / 3) - 1) + ')'
                                    # Nucleotide variant range spans the Termination codon
                                    else:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                            ins_thr[:3]) + 'ext?)'

                                # Nucleotide variation has not affected the length of the protein thus substitution or del
                                else:
                                    if len(ins_thr) == 3:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + ins_thr + ')'
                                    elif len(ins_thr) == 0:
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'del)'
                                    else:
                                        posedit = '(' + from_aa + str(
                                            pro_inv_info['edit_start']) + 'delins' + ins_thr + ')'

                            # Complete the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                type='p',
                                                                                posedit=posedit)

                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein

            else:
                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = shifts

            # Return
            return hgvs_transcript_to_hgvs_protein


    # Handle non-coding transcript and non transcript descriptions
    elif hgvs_transcript.type == 'n':
        # non-coding transcripts
        hgvs_protein = copy.deepcopy(hgvs_transcript)
        hgvs_protein.ac = 'Non-coding '
        hgvs_protein.posedit = ''
        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
        return hgvs_transcript_to_hgvs_protein
    else:
        hgvs_transcript_to_hgvs_protein['error'] = 'Unable to map %s to %s' % (
            hgvs_transcript.ac, associated_protein_accession)
        return hgvs_transcript_to_hgvs_protein


"""
Returns exon information for a given transcript
e.g. how the exons align to the genomic reference
see hgvs.dataproviders.uta.py for details
"""


def tx_exons(tx_ac, alt_ac, alt_aln_method, hdp):
    # Interface with the UTA database via get_tx_exons in uta.py
    try:
        tx_exons = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
    except hgvs.exceptions.HGVSError as e:
        tx_exons = 'hgvs Exception: ' + str(e)
        return tx_exons
    try:
        tx_exons[0]['alt_strand']
    except TypeError:
        tx_exons = 'error'
        return tx_exons
    # If on the reverse strand, reverse the order of elements
    if tx_exons[0]['alt_strand'] == -1:
        tx_exons = tx_exons[::-1]
        return tx_exons
    else:
        return tx_exons

"""
Simple reverse complement function for nucleotide sequences
"""

def revcomp(bases):
    l2 = []
    l = list(bases)
    element = 0
    for base in l:
        element = element + 1
        if base == 'G':
            l2.append('C')
        if base == 'C':
            l2.append('G')
        if base == 'A':
            l2.append('T')
        if base == 'T':
            l2.append('A')
    revcomp = ''.join(l2)
    revcomp = revcomp[::-1]
    return revcomp
    
# <LICENSE>
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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
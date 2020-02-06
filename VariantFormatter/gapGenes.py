# -*- coding: utf-8 -*-

# Import python modules
# from distutils.version import StrictVersion
import re
import copy
import vvhgvs.exceptions
import vvhgvs.assemblymapper
import vvhgvs.variantmapper
import VariantFormatter

# VV
import VariantValidator
import VariantValidator.modules.seq_data as seq_data
import VariantValidator.modules.hgvs_utils as va_H2V

"""
Internal function that returns True if gene symbol is blacklisted and False if not
"""
def gap_black_list(symbol):
    gap_check = seq_data.gap_black_list(symbol)
    return gap_check



"""
Head function for g_to_t mapping
Runs quick tests to see if gap compensation is required

RefSeq only
hgvs <= 1.1.3
"""


def compensate_g_to_t(hgvs_tx, hgvs_genomic, un_norm_hgvs_genomic, vm,
                        hn, reverse_normalizer, primary_assembly, hdp, hp, sf, hgvs_version, vfo):

    # Not required in these instances
    if re.match('ENST', hgvs_tx.ac): # or (StrictVersion(str(hgvs_version)) > StrictVersion('1.1.3') is True):
        
        # Push to absolute position
        normalized_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn, reverse_normalizer,
                                             hdp, vm, vfo)
        hgvs_tx_returns = [normalized_tx, False, None, None, None]

    else:
        gene_symbol = hdp.get_tx_identity_info(hgvs_tx.ac)[6]
        # Check the blacklist
        gap_compensation = gap_black_list(gene_symbol)
        # print('Is gapped = ' + str(gap_compensation))
        if gap_compensation is False:
            normalized_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn, 
                                                reverse_normalizer, hdp, vm, vfo)
            hgvs_tx_returns = [normalized_tx, False, None, None, None]
        else:
            # At this stage, we know that:
            # the gene is on the gap list
            # the hgvs version is <= 1.1.3
            # The requested transcript set is RefSeq
            gap_compensated_tx = g_to_t_compensation_code(hgvs_tx, hgvs_genomic,
                                                                un_norm_hgvs_genomic, 
                                                                vm, hn, 
                                                                reverse_normalizer,
                                                                primary_assembly,
                                                                hdp, hp, sf, vfo)
            # except Exception as e:
            #     import sys
            #     import traceback
            #     exc_type, exc_value, last_traceback = sys.exc_info()
            #     te = traceback.format_exc()
            #     tbk = [str(exc_type), str(exc_value), str(te)]
            #     er = str('\n'.join(tbk))

            if gap_compensated_tx[1] is False:
                refresh_hgvs_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn,
                                                reverse_normalizer, hdp, vm, vfo)
                gap_compensated_tx[0] = refresh_hgvs_tx                                         
            hgvs_tx_returns = gap_compensated_tx
                        
    
    hgvs_tx_dict = {'hgvs_transcript': hgvs_tx_returns[0],
                    'position_lock': hgvs_tx_returns[1],
                    'gapped_alignment_warning': hgvs_tx_returns[3],
                    'corrective_action': hgvs_tx_returns[2],
                    'gap_position': hgvs_tx_returns[-1],
                    'transcript_accession': hgvs_tx_returns[0].ac
                    }
    return hgvs_tx_dict


"""
Fully normalizes the hgvs_tx variant from the hgvs_genomic variant

Is only activated if the g_to_t_compensation_code IS NOT USED!
"""


def fully_normalize(hgvs_tx, hgvs_genomic, hn, reverse_normalizer, hdp, vm, vfo):
    
    # set required variables
    tx_id = hgvs_tx.ac
    if re.match('ENST', tx_id):
        alt_aln_method = 'genebuild'
    else:
        alt_aln_method = 'splign'
    rhn = reverse_normalizer

    # Obtain the orientation of the transcript wrt selected genomic accession
    exon_alignments = vfo.tx_exons(tx_id, hgvs_genomic.ac, alt_aln_method)
    orientation = int(exon_alignments[0]['alt_strand'])
    # Normalize the genomic variant 5 prime if antisense or 3 prime if sense
    if orientation == -1:
        hgvs_genomic = rhn.normalize(hgvs_genomic)
    else:
        hgvs_genomic = hn.normalize(hgvs_genomic)

    # print('Try 1')
    try:
        hgvs_tx = vm.g_to_t(hgvs_genomic, hgvs_tx.ac)
    except vvhgvs.exceptions.HGVSError as e:
        # print('Error area 1')
        # print(e)
        pass
    # print('Try 2')
    try:
        hgvs_tx = hn.normalize(hgvs_tx)
    except vvhgvs.exceptions.HGVSError as e:
        # print('Error area 2')
        # print(e)
        pass
        
    return hgvs_tx


"""
Gap compensation code from genome to transcript

Source is VariantValidator

Requires an un-normalized genomic variant for stashing if available

Also requres a hgvs_genomic and tx id
"""
                        
def g_to_t_compensation_code(hgvs_tx, hgvs_genomic, un_norm_hgvs_genomic, vm, hn, 
                                reverse_normalizer, primary_assembly, hdp, hp, sf, vfo):

    # Ensure hgvs_tx is good to go
    try:
        hn.normalize(hgvs_tx)
    except vvhgvs.exceptions.HGVSInvalidVariantError as e:
        if 'insertion length must be 1' in str(e):
            hgvs_tx_anew = '%s:%s.%sdelins%s' % (hgvs_tx.ac, hgvs_tx.type, str(hgvs_tx.posedit.pos), hgvs_tx.posedit.edit.alt)
            hgvs_tx = hp.parse_hgvs_variant(hgvs_tx_anew)
    except vvhgvs.exceptions.HGVSUnsupportedOperationError:
        pass

    """
    Gap aware projection from g. to c.
    """

    # Create mappers: has to be inline because requires primary_assembly and alt_aln_method
    alt_aln_method = 'splign'
    no_norm_evm = vvhgvs.assemblymapper.AssemblyMapper(hdp,
                                                     assembly_name=primary_assembly,
                                                     alt_aln_method=alt_aln_method, # Only RefSeq should be here!!!
                                                     normalize=False,
                                                     replace_reference=True
                                                     )
    evm = vvhgvs.assemblymapper.AssemblyMapper(hdp,
                                                     assembly_name=primary_assembly,
                                                     alt_aln_method=alt_aln_method, # Only RefSeq should be here!!!
                                                     normalize=True,
                                                     replace_reference=True
                                                     )

    nr_vm = vvhgvs.variantmapper.VariantMapper(hdp, replace_reference=False)

    utilise_gap_code = True
    automap = ''
    
    # Set variables for problem specific warnings
    gapped_alignment_warning = None
    corrective_action_taken = None
    gapped_transcripts = ''
    auto_info = None
    
    # Set variables
    stash_input = hgvs_genomic
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Create a pseudo VCF so that normalization can be applied and a delins can be generated
    vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly, reverse_normalizer, sf)
    chr = vcf_dict['chr']
    pos = vcf_dict['pos']
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']

    # Generate an end position
    end = str(int(pos) + len(ref) - 1)
    pos = str(pos)

    # take a look at the input genomic variant for potential base salvage
    stash_ac = vcf_dict['chr']
    stash_pos = int(vcf_dict['pos'])
    stash_ref = vcf_dict['ref']
    stash_alt = vcf_dict['alt']
    stash_end = end
    # Re-Analyse genomic positions
    if re.match('NG_', str(stash_input)):
        c = hgvs_tx
        try:
            c.posedit.edit.ref = c.posedit.edit.ref.upper()
        except Exception:
            pass
        try:
            c.posedit.edit.alt = c.posedit.edit.alt.upper()
        except Exception:
            pass
        stash_input = vfo.myevm_t_to_g(c, no_norm_evm, primary_assembly, hn)
    if re.match('NC_', str(stash_input)) or re.match('NT_', str(stash_input)) or re.match('NW_',
                                                                                          str(
                                                                                                  stash_input)):
        try:
            hgvs_stash = hp.parse_hgvs_variant(stash_input)
        except:
            hgvs_stash = stash_input
        try:
            hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
        except Exception:
            pass
        try:
            hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()
        except Exception:
            pass

        stash_ac = hgvs_stash.ac
        # MAKE A NO NORM HGVS2VCF
        stash_dict = va_H2V.pos_lock_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer, sf)

        stash_ac = hgvs_stash.ac
        stash_pos = int(stash_dict['pos'])
        stash_ref = stash_dict['ref']
        stash_alt = stash_dict['alt']
        # Generate an end position
        stash_end = str(stash_pos + len(stash_ref) - 1)

    # Store a not real deletion insertion
    stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
        hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
    stash_hgvs_not_delins = hp.parse_hgvs_variant(
        stash_ac + ':' + hgvs_genomic_5pr.type + '.' + str(
            stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)

    # Set non-valid caution to false
    non_valid_caution = 'false'

    # Store the current hgvs:c. description
    saved_hgvs_coding = hgvs_tx

    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
    ori = vfo.tx_exons(saved_hgvs_coding.ac, hgvs_genomic.ac, alt_aln_method)
    orientation = int(ori[0]['alt_strand'])
    intronic_variant = 'false'

    if orientation == -1:
        # position genomic at its most 5 prime position
        try:
            query_genomic = reverse_normalizer.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
        # Map to the transcript ant test for movement
        try:
            hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
        except vvhgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = str(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
            pass
        else:
            hgvs_seek_var = saved_hgvs_coding

    elif orientation != -1:
        # position genomic at its most 3 prime position
        try:
            query_genomic = hn.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
    # Map to the transcript and test for movement
    try:
        hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
    except vvhgvs.exceptions.HGVSError as e:
        hgvs_seek_var = saved_hgvs_coding
    else:
        #seek_var = str(hgvs_seek_var)
        #seek_ac = str(hgvs_seek_var.ac)
        if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
            pass
        else:
            hgvs_seek_var = saved_hgvs_coding

    try:
        intron_test = hn.normalize(hgvs_seek_var)
    except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
        error = str(e)
        if re.match('Normalization of intronic variants is not supported', error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
            if re.match(
                    'Unsupported normalization of variants spanning the exon-intron boundary',
                    error):
                intronic_variant = 'hard_fail'
            else:
                # Double check to see whether the variant is actually intronic?
                for exon in ori:
                    genomic_start = int(exon['alt_start_i'])
                    genomic_end = int(exon['alt_end_i'])
                    if (
                            hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                            hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                        intronic_variant = 'false'
                        break
                    else:
                        intronic_variant = 'true'

    if intronic_variant != 'hard_fail':
        if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
            hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-',
                                                     str(hgvs_seek_var.posedit.pos)):
            # Double check to see whether the variant is actually intronic?
            for exon in ori:
                genomic_start = int(exon['alt_start_i'])
                genomic_end = int(exon['alt_end_i'])
                if (
                        hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                        hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                    intronic_variant = 'false'
                    break
                else:
                    intronic_variant = 'true'

    if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
            hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
        hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
        # Double check to see whether the variant is actually intronic?
        for exon in ori:
            genomic_start = int(exon['alt_start_i'])
            genomic_end = int(exon['alt_end_i'])
            if (
                    hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                    hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                intronic_variant = 'false'
                break
            else:
                intronic_variant = 'true'

    # If exonic, process
    if intronic_variant != 'true':
        # map form reverse normalized g. to c.
        hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)

        # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
        disparity_deletion_in = ['false', 'false']
        if stored_hgvs_not_delins != '':
            # Refresh hgvs_not_delins from stored_hgvs_not_delins
            hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
            # This test will only occur in dup of single base, insertion or substitution
            if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                     hgvs_genomic_5pr.posedit.edit.type):
                    # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                    plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                    plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                    plussed_hgvs_not_delins.posedit.edit.ref = ''
                    transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                            str(saved_hgvs_coding.ac))
                    if ((
                            transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                            hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                        elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        elif re.search('ins',
                                       str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                            'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                    else:
                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                        elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        elif re.search('ins',
                                       str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                            'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                else:
                    pass
            else:
                pass

            try:
                tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
            except vvhgvs.exceptions.HGVSInvalidIntervalError:
                tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_genomic_5pr, saved_hgvs_coding.ac)
            except vvhgvs.exceptions.HGVSError:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    tx_hgvs_not_delins = saved_hgvs_coding

            # Create normalized version of tx_hgvs_not_delins
            rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
            # Check for +ve base and adjust
            if (re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                           str(
                                                                                               rn_tx_hgvs_not_delins.posedit.pos.start))) and (
                    re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                '\-', str(rn_tx_hgvs_not_delins.posedit.pos.end))):
                # Remove offsetting to span the gap
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                try:
                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                except:
                    pass
            elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                # move tx end base to next available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = vfo.myevm_t_to_g(test_tx_var, no_norm_evm, primary_assembly, hn)

                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
            elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                # move tx start base to previous available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = vfo.myevm_t_to_g(test_tx_var, no_norm_evm, primary_assembly, hn)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0

            # Check for -ve base and adjust
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-',
                                                                                           str(
                                                                                               rn_tx_hgvs_not_delins.posedit.pos.start)):
                # Remove offsetting to span the gap
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                try:
                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                except:
                    pass
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                # move tx end base back to next available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                # Delete the ref
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                # Add the additional base to the ALT
                start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = vfo.myevm_t_to_g(test_tx_var, no_norm_evm, primary_assembly, hn)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                # move tx start base to previous available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = vfo.myevm_t_to_g(test_tx_var, no_norm_evm, primary_assembly, hn)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
            else:
                pass

            # Logic
            if len(hgvs_not_delins.posedit.edit.ref) < len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                    hgvs_not_delins.posedit.edit.ref)
                disparity_deletion_in = ['chromosome', gap_length]
            elif len(hgvs_not_delins.posedit.edit.ref) > len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref)
                disparity_deletion_in = ['transcript', gap_length]
            else:
                hgvs_stash_t = vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac)
                if len(stash_hgvs_not_delins.posedit.edit.ref) > len(
                        hgvs_stash_t.posedit.edit.ref):
                    try:
                        hn.normalize(hgvs_stash_t)
                    except:
                        pass
                    else:
                        gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                            hgvs_stash_t.posedit.edit.ref)
                        disparity_deletion_in = ['transcript', gap_length]
                        try:
                            tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                        except:
                            tx_hgvs_not_delins = hgvs_stash_t
                        hgvs_not_delins = stash_hgvs_not_delins
                elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                    disparity_deletion_in = ['transcript', 'Requires Analysis']
                    try:
                        tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                    except:
                        tx_hgvs_not_delins = hgvs_stash_t
                    hgvs_not_delins = stash_hgvs_not_delins
                    hgvs_genomic_5pr = stash_hgvs_not_delins
                else:
                    pass

        # Final sanity checks
        try:
            vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
        except Exception as e:
            if str(e) == 'start or end or both are beyond the bounds of transcript record':
                hgvs_not_delins = saved_hgvs_coding
                disparity_deletion_in = ['false', 'false']
        try:
            hn.normalize(tx_hgvs_not_delins)
        except vvhgvs.exceptions.HGVSUnsupportedOperationError as e:
            error = str(e)
            if re.match('Normalization of intronic variants is not supported',
                        error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
                if re.match(
                        'Unsupported normalization of variants spanning the exon-intron boundary',
                        error):
                    hgvs_not_delins = saved_hgvs_coding
                    disparity_deletion_in = ['false', 'false']
                elif re.match('Normalization of intronic variants is not supported', error):
                    # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                    disparity_deletion_in = ['transcript', 'Requires Analysis']

        # Pre-processing of tx_hgvs_not_delins
        try:
            if tx_hgvs_not_delins.posedit.edit.alt is None:
                tx_hgvs_not_delins.posedit.edit.alt = ''
        except Exception as e:
            if str(e) == "'Dup' object has no attribute 'alt'":
                tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                    tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                    tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                tx_hgvs_not_delins = hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

        # GAP IN THE TRANSCRIPT DISPARITY DETECTED
        if disparity_deletion_in[0] == 'transcript':
            if disparity_deletion_in[1] == 'Requires Analysis':
                analyse_gap = copy.deepcopy(tx_hgvs_not_delins)
                try:
                    analyse_gap = vm.n_to_c(analyse_gap)
                except vvhgvs.exceptions.HGVSError:
                    pass
                analyse_gap.posedit.pos.start.offset = 0
                analyse_gap.posedit.pos.end.offset = 0
                try:
                    analyse_gap.posedit.edit.ref = ''
                except AttributeError:
                    pass
                try:
                    analyse_gap.posedit.edit.alt = ''
                except AttributeError:
                    pass
                my_g_gap_v = vm.t_to_g(analyse_gap, hgvs_genomic_5pr.ac)
                g_gap = my_g_gap_v.posedit.pos.end.base - my_g_gap_v.posedit.pos.start.base + 1
                my_t_gap_v = vm.g_to_t(my_g_gap_v, analyse_gap.ac)
                t_gap = my_t_gap_v.posedit.pos.end.base - my_t_gap_v.posedit.pos.start.base
                disparity_deletion_in[1] = str(g_gap - t_gap - 1)

            gap_position = ''
            gapped_alignment_warning = 'The displayed variants may be artefacts of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly

            # ANY VARIANT WHOLLY WITHIN THE GAP
            if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                        str(
                                                                                            tx_hgvs_not_delins.posedit.pos.start))) and (
                    re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-',
                                                                                          str(
                                                                                              tx_hgvs_not_delins.posedit.pos.end))):
                gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                # Copy the current variant
                tx_gap_fill_variant = copy.deepcopy(tx_hgvs_not_delins)
                try:
                    if tx_gap_fill_variant.posedit.edit.alt is None:
                        tx_gap_fill_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        tx_gap_fill_variant_delins_from_dup = tx_gap_fill_variant.ac + ':' + tx_gap_fill_variant.type + '.' + str(
                            tx_gap_fill_variant.posedit.pos.start) + '_' + str(
                            tx_gap_fill_variant.posedit.pos.end) + 'del' + tx_gap_fill_variant.posedit.edit.ref + 'ins' + tx_gap_fill_variant.posedit.edit.ref + tx_gap_fill_variant.posedit.edit.ref
                        tx_gap_fill_variant = hp.parse_hgvs_variant(
                            tx_gap_fill_variant_delins_from_dup)

                # Identify which half of the NOT-intron the start position of the variant is in
                if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                    tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                    tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                    tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                    tx_gap_fill_variant.posedit.edit.alt = ''
                    tx_gap_fill_variant.posedit.edit.ref = ''
                elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                    tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                    tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                    tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                    tx_gap_fill_variant.posedit.edit.alt = ''
                    tx_gap_fill_variant.posedit.edit.ref = ''

                try:
                    tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                except:
                    pass
                genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                     reverse_normalized_hgvs_genomic.ac)
                genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                try:
                    c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                except Exception:
                    c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
                                                         hgvs_genomic_5pr.ac)

                # Ensure an ALT exists
                try:
                    if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                        genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        genomic_gap_fill_variant_delins_from_dup = genomic_gap_fill_variant.ac + ':' + genomic_gap_fill_variant.type + '.' + str(
                            genomic_gap_fill_variant.posedit.pos.start.base) + '_' + str(
                            genomic_gap_fill_variant.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant.posedit.edit.ref + 'ins' + genomic_gap_fill_variant.posedit.edit.ref + genomic_gap_fill_variant.posedit.edit.ref
                        genomic_gap_fill_variant = hp.parse_hgvs_variant(
                            genomic_gap_fill_variant_delins_from_dup)
                        genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                            genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                            genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                        genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                            genomic_gap_fill_variant_alt_delins_from_dup)

                # Correct insertion alts
                if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                    append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                              genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                              genomic_gap_fill_variant_alt.posedit.pos.end.base)
                    genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                        0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                    append_ref[1]

                # Split the reference and replacing alt sequence into a dictionary
                reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                    alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
                else:
                    # Deletions with no ins
                    pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
                    alternate_bases = []
                    for base in pre_alternate_bases:
                        alternate_bases.append('X')

                # Create the dictionaries
                ref_start = genomic_gap_fill_variant.posedit.pos.start.base
                alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
                ref_base_dict = {}
                for base in reference_bases:
                    ref_base_dict[ref_start] = str(base)
                    ref_start = ref_start + 1

                alt_base_dict = {}

                # NEED TO SEARCH FOR RANGE = and replace with interval_range
                # Need to search for int and replace with integer

                # Note, all variants will be forced into the format delete insert
                # Deleted bases in the ALT will be substituted for X
                for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                     genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                    if integer == alt_start:
                        alt_base_dict[integer] = str(''.join(alternate_bases))
                    else:
                        alt_base_dict[integer] = 'X'

                # Generate the alt sequence
                alternate_sequence_bases = []
                for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                     genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                    if integer in alt_base_dict.keys():
                        alternate_sequence_bases.append(alt_base_dict[integer])
                    else:
                        alternate_sequence_bases.append(ref_base_dict[integer])
                alternate_sequence = ''.join(alternate_sequence_bases)
                alternate_sequence = alternate_sequence.replace('X', '')

                # Add the new alt to the gap fill variant and generate transcript variant
                genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
                                                   tx_gap_fill_variant.ac)

                # Set warning
                gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                disparity_deletion_in[1] = [gap_size]
                auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                    stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                    tx_hgvs_not_delins.ac)
                non_valid_caution = 'true'

                # Alignment position
                for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                if re.match('NM_', str(for_location_c)):
                    for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
                    gps = for_location_c.posedit.pos.start.base - 1
                    gpe = for_location_c.posedit.pos.start.base
                else:
                    gps = for_location_c.posedit.pos.start.base
                    gpe = for_location_c.posedit.pos.start.base + 1
                gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                # auto_info = '%s' % (gap_position)
            else:
                if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                    # In this instance, we have identified a transcript gap but the n. version of
                    # the transcript variant but do not have a position which actually hits the gap,
                    # so the variant likely spans the gap, and is not picked up by an offset.
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                    g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                    ng2 = hn.normalize(g2)
                    g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                                len(g3.posedit.edit.ref) - 1)
                    try:
                        c2 = vm.g_to_t(g3, c1.ac)
                        if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                            pass
                        else:
                            tx_hgvs_not_delins = c2
                            try:
                                tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                            except vvhgvs.exceptions.HGVSError:
                                pass
                    except vvhgvs.exceptions.HGVSInvalidVariantError:
                        pass

                if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                        '\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c2 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c2 = tx_hgvs_not_delins
                    c1 = copy.deepcopy(c2)
                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                    c1.posedit.pos.start.offset = 0
                    c1.posedit.pos.end = c2.posedit.pos.start
                    c1.posedit.edit.ref = ''
                    c1.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                        gps = for_location_c.posedit.pos.start.base
                        gpe = for_location_c.posedit.pos.start.base + 1
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    # auto_info = '%s' % (gap_position)
                elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                        '\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                    auto_info = 'Genome position ' + str(
                        stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.end.base) + ' aligns within a gap in transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ''
                    c2.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.end.base
                    gpe = for_location_c.posedit.pos.end.base + 1
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    # auto_info = '%s' % (gap_position)
                elif re.search('\-',
                               str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                    '\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c2 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c2 = tx_hgvs_not_delins
                    c1 = copy.deepcopy(c2)
                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                    c1.posedit.pos.start.offset = 0
                    c1.posedit.pos.end = c2.posedit.pos.start
                    c1.posedit.edit.ref = ''
                    c1.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.start.base - 1
                    gpe = for_location_c.posedit.pos.start.base
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    # auto_info = '%s' % (gap_position)
                elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                        '\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                    auto_info = 'Genome position ' + str(
                        stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.end.base) + ' aligns within a gap in transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ''
                    c2.posedit.edit.alt = ''
                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                    c2 = vm.g_to_t(g2, c2.ac)
                    reference = c1.posedit.edit.ref + c2.posedit.edit.ref[1:]
                    alternate = c1.posedit.edit.alt + c2.posedit.edit.ref[1:]
                    c3 = copy.deepcopy(c1)
                    c3.posedit.pos.end = c2.posedit.pos.end
                    c3.posedit.edit.ref = ''  # reference
                    c3.posedit.edit.alt = alternate
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.end.base - 1
                    gpe = for_location_c.posedit.pos.end.base
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    # auto_info = '%s' % (gap_position)
                else:
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    hgvs_refreshed_variant = tx_hgvs_not_delins
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

        # GAP IN THE CHROMOSOME
        elif disparity_deletion_in[0] == 'chromosome':
            # Set warning variables
            gap_position = ''
            gapped_alignment_warning = 'The displayed variants may be artefacts of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
            hgvs_refreshed_variant = tx_hgvs_not_delins
            # Warn
            auto_info = str(hgvs_refreshed_variant.ac) + ':c.' + str(
                hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(disparity_deletion_in[
                                                                             1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                hgvs_genomic.ac)
            gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
        else:
            # Try the push
            hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
            stash_ac = hgvs_stash.ac
            # Make a hard left and hard right not delins g.
            stash_dict_right = va_H2V.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, sf)
            stash_pos_right = int(stash_dict_right['pos'])
            stash_ref_right = stash_dict_right['ref']
            stash_alt_right = stash_dict_right['alt']
            stash_end_right = str(stash_pos_right + len(stash_ref_right) - 1)
            stash_hgvs_not_delins_right = hp.parse_hgvs_variant(stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos_right) + '_' + stash_end_right + 'del' + stash_ref_right + 'ins' + stash_alt_right)
            stash_dict_left = va_H2V.hard_left_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer, sf)
            stash_pos_left = int(stash_dict_left['pos'])
            stash_ref_left = stash_dict_left['ref']
            stash_alt_left = stash_dict_left['alt']
            stash_end_left = str(stash_pos_left + len(stash_ref_left) - 1)
            stash_hgvs_not_delins_left = hp.parse_hgvs_variant(stash_ac + ':' + hgvs_stash.type + '.' + str(
                    stash_pos_left) + '_' + stash_end_left + 'del' + stash_ref_left + 'ins' + stash_alt_left)
            # Map in-situ to the transcript left and right
            try:
                tx_hard_right = vm.g_to_t(stash_hgvs_not_delins_right, saved_hgvs_coding.ac)
            except Exception as e:
                tx_hard_right = saved_hgvs_coding
            else:
                normalize_stash_right = hn.normalize(stash_hgvs_not_delins_right)
                if str(normalize_stash_right.posedit) == str(stash_hgvs_not_delins.posedit):
                    tx_hard_right = saved_hgvs_coding
            try:
                tx_hard_left = vm.g_to_t(stash_hgvs_not_delins_left, saved_hgvs_coding.ac)
            except Exception as e:
                tx_hard_left = saved_hgvs_coding
            else:
                normalize_stash_left = hn.normalize(stash_hgvs_not_delins_left)
                if str(normalize_stash_left.posedit) == str(stash_hgvs_not_delins.posedit):
                    tx_hard_left = saved_hgvs_coding
            # The Logic - Currently limited to genome gaps
            if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                    tx_hard_right.posedit.edit.ref):
                tx_hard_right = hn.normalize(tx_hard_right)
                gap_position = ''
                gapped_alignment_warning = 'The displayed variants may be artefacts of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
                hgvs_refreshed_variant = tx_hard_right
                gapped_transcripts = gapped_transcripts + str(tx_hard_right.ac) + ' '
            elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < len(
                    tx_hard_left.posedit.edit.ref):
                tx_hard_left = hn.normalize(tx_hard_left)
                gap_position = ''
                gapped_alignment_warning = 'The displayed variants may be artefacts of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
                hgvs_refreshed_variant = tx_hard_left
                gapped_transcripts = gapped_transcripts + str(tx_hard_left.ac) + ' '
            else:
                # Keep the same by re-setting rel_var
                hgvs_refreshed_variant = saved_hgvs_coding

        # Edit the output
        if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
                hgvs_refreshed_variant.type)):
            hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant)
        else:
            pass
        
        # Set pos_lock
        pos_lock = True
        
        try:
            hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
            if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                    hgvs_refreshed_variant.posedit.edit.alt[-1]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
            elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[0] == \
                    hgvs_refreshed_variant.posedit.edit.alt[0]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          1:]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          1:]
                hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
                
        except Exception as e:
            error = str(e)
            # Ensure the final variant is not intronic nor does it cross exon boundaries
            if re.match('Normalization of intronic variants is not supported',
                        error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
                hgvs_refreshed_variant = saved_hgvs_coding
                corrective_action_taken = None
                gapped_alignment_warning = None
                auto_info = None
                pos_lock = False
            else:
                pass

    # Otherwise these variants need to be set
    else:
        corrective_action_taken = None
        gapped_alignment_warning = None
        auto_info = None
        pos_lock = False
        hgvs_refreshed_variant = saved_hgvs_coding

    # Warn the user that the g. description is not valid
    if gapped_alignment_warning is not None:
        if disparity_deletion_in[0] == 'transcript':
            corrective_action_taken = 'Automap has deleted ' + str(
                disparity_deletion_in[1]) + ' bp from chromosomal reference sequence ' + str(
                hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence'
        if disparity_deletion_in[0] == 'chromosome':
            corrective_action_taken = 'Automap has added ' + str(
                disparity_deletion_in[1]) + ' bp to chromosomal reference sequence ' + str(
                hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence '

    # Add additional data to the front of automap
    if auto_info is not None:
        automap = auto_info + '\n' + automap
    
    # Make the return
    gap_report = [hgvs_refreshed_variant, pos_lock, corrective_action_taken, 
                    gapped_alignment_warning, auto_info]

    return gap_report
    
    
# <LICENSE>
# Copyright (C) 2019 VariantValidator Contributors
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

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
import sys
import copy
import warnings


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


def myevm_t_to_g(hgvs_c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm):

    # store the input
    stored_hgvs_c = copy.deepcopy(hgvs_c)
    expand_out = 'false'
    utilise_gap_code = True
    
    # Gap gene black list
    try:
        gene_symbol = dbControls.data.get_gene_symbol_from_transcriptID(hgvs_c.ac)
    except Exception:
        utilise_gap_code = False
    else:
        # If the gene symbol is not in the list, the value False will be returned
        utilise_gap_code = gap_genes.gap_black_list(gene_symbol)
    # Warn gap code in use
    warnings.warn("gap_compensation_myevm = " + str(utilise_gap_code))
        
    if utilise_gap_code is True and (hgvs_c.posedit.edit.type == 'identity' or hgvs_c.posedit.edit.type == 'del' or hgvs_c.posedit.edit.type =='delins' or hgvs_c.posedit.edit.type == 'dup' or hgvs_c.posedit.edit.type == 'sub' or hgvs_c.posedit.edit.type == 'ins' or hgvs_c.posedit.edit.type == 'inv'):

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
                    t_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(hgvs_t.posedit.pos.end.base) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t_delins = hp.parse_hgvs_variant(t_delins)
                    pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    inv_alt = pre_base + inv_alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + 'del' + hgvs_t.posedit.edit.ref + 'ins' + inv_alt
                    hgvs_t = hp.parse_hgvs_variant(hgvs_str)                    
                elif hgvs_c.posedit.edit.type == 'dup':
                    pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
                    alt = pre_base + hgvs_t.posedit.edit.ref + hgvs_t.posedit.edit.ref + post_base              
                    ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    dup_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base - 1) + '_' + str((hgvs_t.posedit.pos.start.base + len(ref)) -2) + 'del' + ref + 'ins' + alt
                    hgvs_t = hp.parse_hgvs_variant(dup_to_delins)
                elif hgvs_c.posedit.edit.type == 'ins': 
                    ins_ref = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.end.base+1)
                    ins_alt = ins_ref[:2] + hgvs_t.posedit.edit.alt + ins_ref[-2:]
                    ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base - 1) + '_' + str(hgvs_t.posedit.pos.end.base +1 ) + 'del' + ins_ref + 'ins' + ins_alt
                    hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                else:   
                    if str(hgvs_t.posedit.edit.alt) == 'None':
                        hgvs_t.posedit.edit.alt = ''    
                    pre_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-2,hgvs_t.posedit.pos.start.base-1)
                    post_base = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.end.base,hgvs_t.posedit.pos.end.base+1)
                    hgvs_t.posedit.edit.ref = pre_base + hgvs_t.posedit.edit.ref + post_base
                    hgvs_t.posedit.edit.alt = pre_base + hgvs_t.posedit.edit.alt + post_base
                    hgvs_t.posedit.pos.start.base = hgvs_t.posedit.pos.start.base - 1
                    start = hgvs_t.posedit.pos.start.base
                    hgvs_t.posedit.pos.start.base = start + 1
                    hgvs_t.posedit.pos.end.base = hgvs_t.posedit.pos.end.base + 1
                    end = hgvs_t.posedit.pos.end.base
                    hgvs_t.posedit.pos.start.base = start
                    hgvs_t.posedit.pos.end.base = end
                    hgvs_str = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(start) + '_' + str(end) + str(hgvs_t.posedit.edit)
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
            reform_ident = reform_ident + ':' + stored_hgvs_c.type + '.' + str(hgvs_c.posedit.pos) + 'del' + str(hgvs_c.posedit.edit.ref)# + 'ins' + str(hgvs_c.posedit.edit.alt)
            hgvs_reform_ident = hp.parse_hgvs_variant(reform_ident)
            try:
                hn.normalize(hgvs_reform_ident)
            except hgvs.exceptions.HGVSError as e:
                error = str(e)
                if re.search('spanning the exon-intron boundary', error) or re.search('Normalization of intronic variants', error):
                    hgvs_c = copy.deepcopy(stored_hgvs_c)
    try:
        hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
        hn.normalize(hgvs_genomic) # Check the validity of the mapping
        # This will fail on multiple refs for NC_
    except hgvs.exceptions.HGVSError as e:
        # Recover all available mapping options from UTA
        mapping_options = hdp.get_tx_mapping_options(hgvs_c.ac) 
        
        if mapping_options == []:
            raise HGVSDataNotAvailableError("No alignment data between the specified transcript reference sequence and any GRCh37 and GRCh38 genomic reference sequences (including alternate chromosome assemblies, patches and RefSeqGenes) are available.")
        
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
                        attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
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
                    try:
                        hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                        break
                    except Exception as e:
                        if re.search(option[1], attempted_mapping_error):
                            pass
                        else:
                            attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
                        print e
                        continue
            try:
                hn.normalize(hgvs_genomic)
            except:
                for option in mapping_options:
                    if re.match('blat', option[2]):
                        continue
                    if re.match('NT_', option[1]):
                        try:
                            hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                            break
                        except Exception as e:
                            attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
                            print e
                            continue
                try:
                    hn.normalize(hgvs_genomic)
                except:                         
                    for option in mapping_options:
                        if re.match('blat', option[2]):
                            continue
                        if re.match('NW_', option[1]):
                            try:
                                hgvs_genomic = vm.t_to_g(hgvs_c, str(option[1]))
                                break
                            except Exception as e:
                                attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
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
                                    attempted_mapping_error = attempted_mapping_error + str(e) + "/" + hgvs_c.ac + "/" + option[1] + '~'
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
                ref = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.end.base)
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
        stored_ref = sf.fetch_seq(str(stored_hgvs_n.ac),stored_hgvs_n.posedit.pos.start.base-1,stored_hgvs_n.posedit.pos.end.base)
        stored_hgvs_c.posedit.edit.ref = stored_ref

    if (hgvs_genomic.posedit.edit.ref == '' or hgvs_genomic.posedit.edit.ref is None) and expand_out != 'false':
        if hgvs_genomic.posedit.edit.type == 'ins':
            stored_ref = sf.fetch_seq(str(hgvs_genomic.ac),hgvs_genomic.posedit.pos.start.base-1,hgvs_genomic.posedit.pos.end.base) 
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
            if re.match('Length implied by coordinates must equal sequence deletion length', str(e)) or str(e) == 'base start position must be <= end position':
                # Effectively, this code is designed to handle variants that are directly proximal to 
                # gap BOUNDARIES, but in some cases the replace reference function of hgvs mapping has removed bases due to 
                # the deletion length being > the specified range.
        
                # Warn of variant location wrt the gap
                if re.match('Length implied by coordinates must equal sequence deletion length', str(e)):
                    warnings.warn('Variant is proximal to the flank of a genomic gap')
                    genomic_gap_variant = vm.t_to_g(stored_hgvs_c, hgvs_genomic.ac)
                    try:
                        hn.normalize(genomic_gap_variant)
                    except Exception:
                        pass
                    else:
                        genomic_gap_variant = nr_vm.t_to_g(hgvs_c, hgvs_genomic.ac)
                    
                if str(e) == 'base start position must be <= end position':
                    warnings.warn('Variant is fully within a genomic gap')
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
                        warnings.warn('Variant is on the flank of a genomic gap but not within the gap')
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
                            transcript_gap_n_delins_from_dup = transcript_gap_n.ac + ':' + transcript_gap_n.type + '.' + str(transcript_gap_n.posedit.pos.start.base) + '_' + str(transcript_gap_n.posedit.pos.end.base) + 'del' + transcript_gap_n.posedit.edit.ref + 'ins' +  transcript_gap_n.posedit.edit.ref + transcript_gap_n.posedit.edit.ref                           
                            transcript_gap_n = hp.parse_hgvs_variant(transcript_gap_n_delins_from_dup)                          
                            transcript_gap_alt_n_delins_from_dup = transcript_gap_alt_n.ac + ':' + transcript_gap_alt_n.type + '.' + str(transcript_gap_alt_n.posedit.pos.start.base) + '_' + str(transcript_gap_alt_n.posedit.pos.end.base) + 'del' + transcript_gap_alt_n.posedit.edit.ref + 'ins' +  transcript_gap_alt_n.posedit.edit.ref + transcript_gap_alt_n.posedit.edit.ref                           
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
                    for int in range(transcript_gap_alt_n.posedit.pos.start.base, transcript_gap_alt_n.posedit.pos.end.base+1, 1):
                        if int == alt_start:
                            alt_base_dict[int] = str(''.join(alternate_bases))
                        else:
                            alt_base_dict[int] = 'X'        

                    # Generate the alt sequence
                    alternate_sequence_bases = []
                    for int in range(transcript_gap_n.posedit.pos.start.base, transcript_gap_n.posedit.pos.end.base+1, 1):
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
                            pre_base = sf.fetch_seq(transcript_gap_n.ac,transcript_gap_n.posedit.pos.start.base-2,transcript_gap_n.posedit.pos.start.base-1)
                            post_base = sf.fetch_seq(transcript_gap_n.ac,transcript_gap_n.posedit.pos.end.base,transcript_gap_n.posedit.pos.end.base+1)
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
                ins_ref = sf.fetch_seq(str(hgvs_t.ac),hgvs_t.posedit.pos.start.base-1,hgvs_t.posedit.pos.end.base)
                ins_alt = ins_ref[:1] + hgvs_t.posedit.edit.alt + ins_ref[-1:]
                ins_to_delins = hgvs_t.ac + ':' + hgvs_t.type + '.' + str(hgvs_t.posedit.pos.start.base) + '_' + str(hgvs_t.posedit.pos.end.base) + 'del' + ins_ref + 'ins' + ins_alt
                hgvs_t = hp.parse_hgvs_variant(ins_to_delins)
                try:
                    hgvs_c = vm.n_to_c(hgvs_t)
                except Exception:
                    hgvs_c = copy.copy(hgvs_t)
                try:
                    hgvs_genomic = no_norm_evm.t_to_g(hgvs_c)
                except Exception as e:
                    error = str(e)
                    warnings.warn('Ins mapping error in myt_to_g ' + error) 

    return hgvs_genomic
    
    
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
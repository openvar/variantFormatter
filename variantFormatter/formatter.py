# -*- coding: utf-8 -*-

# Import Python modules
import re
import os
from configparser import ConfigParser
import copy

# Import hgvs modules
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import hgvs.exceptions
import hgvs.normalizer
import hgvs.validator
import hgvs.sequencevariant
# import hgvs.location

# import other required modules
from Bio.Seq import Seq

# Import variantFormatter modules
import supportedChromosomeBuilds as chr_dict
import supportFunctions as links
import hgvs2vcf
import gapGenes

# Ensure configuration is on the OS
if os.environ.get('CONF_ROOT') is None:
    import configuration

    CONF_ROOT = os.environ.get('CONF_ROOT')
else:
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
hgvs_version = hgvs.__version__,
hgvs_version = str(hgvs_version[0])

# Pre compile variables
hgvs.global_config.formatting.max_ref_length = 1000000
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect(pooling=True)
vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True)
vr = hgvs.validator.Validator(hdp)
sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
splign_normalizer = hgvs.normalizer.Normalizer(hdp,
                                               cross_boundaries=False,
                                               shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                               alt_aln_method='splign'
                                               )
                                               
genebuild_normalizer = hgvs.normalizer.Normalizer(hdp,
                                                  cross_boundaries=False,
                                                  shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                                  alt_aln_method='genebuild'
                                                  )
                                                  
reverse_splign_normalizer = hgvs.normalizer.Normalizer(hdp,
                                                       cross_boundaries=False,
                                                       shuffle_direction=5,
                                                       alt_aln_method='splign'
                                                       )
                                                       
reverse_genebuild_normalizer = hgvs.normalizer.Normalizer(hdp,
                                                          cross_boundaries=False,
                                                          shuffle_direction=5,
                                                          alt_aln_method='genebuild'
                                                          )
                                                          

"""
Check configuration
"""


def my_config():
    try:
        HGVS_SEQREPO_DIR = os.environ.get('HGVS_SEQREPO_DIR')
    except Exception:
        HGVS_SEQREPO_DIR = 'Unspecified'
    try:
        UTA_DB_URL = os.environ.get('UTA_DB_URL')
    except Exception:
        UTA_DB_URL = 'Unspecified'
    VERSION = __version__,
    VERSION = str(VERSION[0])
    RELEASE = __released__,
    RELEASE = str(RELEASE[0])
    locate = {
        'SeqRepo_Directory': HGVS_SEQREPO_DIR,
        'UTA_URL': UTA_DB_URL,
        'variantFormatter_Version': VERSION,
        'variantFormatter_release_date' : RELEASE,
        'hgvs_version': hgvs_version
    }
    return locate




"""
Simple function, parses an HGVS string into a hgvs (.py) object
NOTE: if the string is not valid, hgvs (py) parser will throw up an error
This function may be best used in a try/except statement
See hgvs (py) documentation for error types
"""


def parse(hgvs_string):
    hgvs_object = hp.parse_hgvs_variant(hgvs_string)
    return hgvs_object
    

"""
Function which takes a vcf string in the format Chr:Pos:Ref:Alt or Chr-Pos-Ref-Alt and the required genome build
and generates a hgvs python object. The object is validated to ensure the stated ref appears at the stated
positions in the RefSeq RefSeq reference sequence

Supported genome builds are GRCh38, GRCh37, hg19, hg38
"""


def vcf2hgvs_genomic(pseudo_vcf, genome_build):
    # Dictionary to store the output
    vcf_to_hgvs_genomic = {'error': '', 'hgvs_genomic': '', 'ref_bases': '', 'un_normalized_hgvs_genomic': ''}
    # Simple check of the format
    if not (re.search(':', pseudo_vcf) or re.search('-', pseudo_vcf)):
        vcf_to_hgvs_genomic['error'] = '%s is an unsupported variant description format' % pseudo_vcf
    else:
        if re.search(':', pseudo_vcf):
            vcf_list = pseudo_vcf.split(':')
        else:
            vcf_list = pseudo_vcf.split('-')
        if len(vcf_list) != 4:
            vcf_to_hgvs_genomic['error'] = '%s is an unsupported variant description format' % pseudo_vcf
        else:
            # extract the variant data
            chrom = vcf_list[0]
            pos = vcf_list[1]
            ref = vcf_list[2]
            alt = vcf_list[3]
            # assemble the HGVS genomic description
            ac = chr_dict.to_accession(chrom, genome_build)
            if ac is None:
                vcf_to_hgvs_genomic[
                    'error'] = 'chromosome ID %s is not associated with genome build %s' % (
                    ac, genome_build)
            else:
                # Create a genomic HGVS string then parse into an hgvs variant object
                gen_var = '%s:g.%s_%sdel%sins%s' % (ac, str(pos), str(int(pos) + (len(ref) - 1)), ref, alt)
                hgvs_genomic = hp.parse_hgvs_variant(gen_var)
                # validate the reference
                try:
                    vr.validate(hgvs_genomic)
                except hgvs.exceptions.HGVSError as e:
                    vcf_to_hgvs_genomic['error'] = str(e)
                else:
                    # Normalize the description to its most 3 prime position and apply the correct HGVS grammer
                    vcf_to_hgvs_genomic['un_normalized_hgvs_genomic'] = hgvs_genomic
                    splign_normalizer.normalize(hgvs_genomic)
                    hgvs_genomic = splign_normalizer.normalize(hgvs_genomic)
                    vcf_to_hgvs_genomic['hgvs_genomic'] = hgvs_genomic
                    vcf_to_hgvs_genomic['ref_bases'] = sf.fetch_seq(hgvs_genomic.ac, start_i=hgvs_genomic.posedit.pos.start.base - 1, end_i=hgvs_genomic.posedit.pos.end.base)
    return vcf_to_hgvs_genomic


"""
Function takes a genomic HGVS description and returns the component parts of a VCF
"""


def hgvs_genomic2vcf(hgvs_genomic, genome_build):
    vcf_dictionary = hgvs2vcf.report_hgvs2vcf(hgvs_genomic, genome_build, hdp, reverse_splign_normalizer, sf)
    return vcf_dictionary


"""
Function takes a parsed hgvs_genomic object and normalizes it
The format is intended to reflect the vcf2hgvs output
The intention is to also collect the un-normalized variant
"""


def format_hgvs_genomic(hgvs_genomic):
    format_hgvs_genomic = {'error': '', 'hgvs_genomic': '', 'ref_bases': '', 'un_normalized_hgvs_genomic': ''}
    
    # Simple check of the format
    try:
        vr.validate(hgvs_genomic)
    except hgvs.exceptions.HGVSError as e:
        format_hgvs_genomic['error'] = str(e)
    else:
        # Normalize the description to its most 3 prime position and apply the correct HGVS grammer
        format_hgvs_genomic['un_normalized_hgvs_genomic'] = hgvs_genomic
        splign_normalizer.normalize(hgvs_genomic)
        hgvs_genomic = splign_normalizer.normalize(hgvs_genomic)
        format_hgvs_genomic['hgvs_genomic'] = hgvs_genomic
        format_hgvs_genomic['ref_bases'] = sf.fetch_seq(hgvs_genomic.ac, start_i=hgvs_genomic.posedit.pos.start.base - 1, end_i=hgvs_genomic.posedit.pos.end.base)
        return format_hgvs_genomic
    

"""
Function which takes a hgvs Python genomic variant and maps to a specified transcript reference sequence. A transcript
level hgvs python object is returned. 
"""


def hgvs_genomic2hgvs_transcript(hgvs_genomic, tx_id):
    # Create dictionary to store the information
    hgvs_genomic_to_hgvs_transcript = {'error': '', 'hgvs_transcript': '', 'ref_bases': ''}
    if re.match('ENST', tx_id):
        alt_aln_method = 'genebuild'
        hn = genebuild_normalizer
        rhn = reverse_genebuild_normalizer
    elif re.match('NM', tx_id) or re.match('NR', tx_id):
        alt_aln_method = 'splign'
        hn = splign_normalizer
        rhn = reverse_splign_normalizer
    # Obtain the orientation of the transcript wrt selected genomic accession
    try:
        exon_alignments = hdp.get_tx_exons(tx_id, hgvs_genomic.ac, alt_aln_method)
        orientation = int(exon_alignments[0]['alt_strand'])

    # May want to tighten up this exception during testing
    except Exception as e:
    	print str(e)
        hgvs_genomic_to_hgvs_transcript['error'] = 'No alignment data available for transcript %s and chromosome %s' % (
            tx_id, hgvs_genomic.ac)
    else:
        # Normalize the genomic variant 5 prime if antisense or 3 prime if sense
        if orientation == -1:
            hgvs_genomic = rhn.normalize(hgvs_genomic)
        else:
            hgvs_genomic = hn.normalize(hgvs_genomic)
        # directly map from the normalized genomic variant to the transcript variant
        try:
            hgvs_tx = vm.g_to_t(hgvs_genomic, tx_id, alt_aln_method)
        except hgvs.exceptions.HGVSError as e:
            hgvs_genomic_to_hgvs_transcript['error'] = str(e)
        else:
            # Ensure complete normalization
            try:
            	hgvs_tx = hn.normalize(hgvs_tx)
            except hgvs.exceptions.HGVSError as e: 
            	pass # No restorative action required. Suspected intronic variant
            else:
				hgvs_genomic_to_hgvs_transcript['hgvs_transcript'] = hgvs_tx
				try:
					hgvs_genomic_to_hgvs_transcript['ref_bases'] = hgvs_tx.posedit.edit.ref
				except hgvs.exceptions.HGVSError:
					if hgvs_tx.type == 'c':
						hgvs_tx = vm.c_to_n(hgvs_t)
					if hgvs_tx.posedit.pos.start.offset == 0 and hgvs_tx.posedit.pos.end.offset == 0:
						hgvs_genomic_to_hgvs_transcript['ref_bases'] = sf.fetch_seq(hgvs_tx.ac, start_i=hgvs_tx.posedit.pos.start.base - 1, end_i=hgvs_tx.posedit.pos.end.base)
					else:
						hgvs_genomic_to_hgvs_transcript['ref_bases'] = ''
                    
    return hgvs_genomic_to_hgvs_transcript


"""
Function which takes a hgvs Python transctipt variant and maps to a specified protein reference sequence. A protein
level hgvs python object is returned.

Note the function currently assumes that the transcript description is correctly normalized having come from the 
previous g_to_t function
"""


def hgvs_transcript2hgvs_protein(hgvs_transcript, associated_protein_accession, re_to_p):
    # Create dictionary to store the information
    hgvs_transcript_to_hgvs_protein = {'error': '', 'hgvs_protein': ''} #, 'hgvs_protein_slc': '', 'ref_residues': ''}
    # Check for non-coding transcripts
    if hgvs_transcript.type == 'c':
        # Handle non inversions with simple c_to_p mapping
                
        if (hgvs_transcript.posedit.edit.type != 'inv') and (hgvs_transcript.posedit.edit.type != 'delins') and (re_to_p is False):
            # Does the edit affect the start codon?
            if ((hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0) or (
                    hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                    and not re.search('\*', str(
                hgvs_transcript.posedit.pos)):
                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                    type='p', posedit='(Met1?)')
            else:
                try:
                    hgvs_protein = vm.c_to_p(hgvs_transcript, associated_protein_accession)
                except IndexError as e:
                    error = str(e)
                    if re.search('string index out of range', error) and re.search('dup', str(hgvs_transcript)):
                        hgvs_ins = hp.parse_hgvs_variant(str(hgvs_transcript))
                        hgvs_ins = hn.normalize(hgvs_ins)
                        inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + str(hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
                        hgvs_transcript = hp.parse_hgvs_variant(inst)
                        hgvs_protein = vm.c_to_p(hgvs_transcript, associated_protein_accession)

            try:
                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                return hgvs_transcript_to_hgvs_protein
            except UnboundLocalError:
                hgvs_transcript_to_hgvs_protein = myc_to_p(hgvs_transcript, associated_protein_accession, re_to_p = True)
                return hgvs_transcript_to_hgvs_protein
                
        else:           
            # Additional code required to process inversions
            # Note, this code was developed for VariantValidator and is not native to the biocommons hgvs Python package
            # Convert positions to n. position
            hgvs_naughty = vm.c_to_n(hgvs_transcript)

            # Collect the deleted sequence using fetch_seq
            del_seq = sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base - 1, end_i=hgvs_naughty.posedit.pos.end.base)

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
                    shifts = vm.c_to_p(hgvs_transcript, associated_protein_accession)
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
                if re.search('\d+\-', str(hgvs_transcript.posedit.pos)) or re.search('\d+\+', str(hgvs_transcript.posedit.pos)) or re.search('\*', str(hgvs_transcript.posedit.pos)) or re.search('[cn].\-', str(hgvs_transcript)):
                    if ((
                                hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                        or
                        (hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                            and not re.search('\*', str(hgvs_transcript.posedit.pos)):                      
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                            posedit='(Met1?)')
                    else:                   
                        # Make the variant
                        hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p', posedit='?')
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
                        hgvs_transcript_to_hgvs_protein['error'] = 'Cannot identify an in-frame Termination codon in the variant mRNA sequence'
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
                            (hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                                and not re.search('\*', str(hgvs_transcript.posedit.pos)):                      
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                posedit='(Met1?)')

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
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
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
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'delins' + ins_thr + ')'

                            # Complete the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
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
Return all aligned transcripts for a given genomic hgvs object
"""


def fetch_aligned_transcripts(hgvs_genomic):
    enst_list = hdp.get_tx_for_region(hgvs_genomic.ac, 'genebuild', hgvs_genomic.posedit.pos.start.base - 1,
                                      hgvs_genomic.posedit.pos.end.base)
    tx_list = enst_list
    refseq_list = hdp.get_tx_for_region(hgvs_genomic.ac, 'splign', hgvs_genomic.posedit.pos.start.base - 1,
                                        hgvs_genomic.posedit.pos.end.base)
    tx_list.append(refseq_list)
    return tx_list

"""
Return the relevant protein sequence for a given transcript reference sequence
"""


def fetch_encoded_protein(tx_ac):
    pro_ac = hdp.get_pro_ac_for_tx_ac(tx_ac)
    return pro_ac
    

"""
format protein description into single letter aa code
"""

def single_letter_protein(hgvs_protein):
    hgvs_protein_slc = hgvs_protein.format({'p_3_letter': False})
    return hgvs_protein_slc
    
    
"""
format nucleotide descriptions to not display reference base
"""

def remove_reference(hgvs_nucleotide):
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless


"""
check RefSeq hgvs_tx descriptions for gaps
"""


def gap_checker(hgvs_transcript, hgvs_genomic, un_norm_hgvs_genomic):
	
    tx_id = hgvs_transcript.ac
    if re.match('ENST', tx_id):
        alt_aln_method = 'genebuild'
        hn = genebuild_normalizer
        rhn = reverse_genebuild_normalizer
    elif re.match('NM', tx_id) or re.match('NR', tx_id):
        alt_aln_method = 'splign'
        hn = splign_normalizer
        rhn = reverse_splign_normalizer

	checked = gapGenes.compensate_g_to_t(hgvs_transcript, hgvs_genomic, 
										un_norm_hgvs_genomic, vm, hn, rhn, 
										hdp, hp, sf, hgvs_version)

	return checked
	

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
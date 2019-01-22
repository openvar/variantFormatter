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

# import other required modules
from Bio.Seq import Seq

# Import variantFormatter modules
import supportedChromosomeBuilds as chr_dict
import hgvs2vcf
import gapGenes
import variantanalyser.functions as va_func

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
    vcf_dictionary = hgvs2vcf.report_hgvs2vcf(hgvs_genomic, genome_build, reverse_splign_normalizer, sf)
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


def hgvs_transcript2hgvs_protein(hgvs_transcript, genome_build):
    
    # Configure evm and normalizers (easier to maintain as now matches VV)
    if 'ENST' in hgvs_transcript.ac:
        alt_aln_method = 'genebuild'
        hn = genebuild_normalizer
        rhn = reverse_genebuild_normalizer
    else:
        alt_aln_method = 'splign'
        hn = splign_normalizer
        rhn = reverse_splign_normalizer 
    
    
    # Create evm
    evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                             assembly_name=genome_build,
                                             alt_aln_method=alt_aln_method, # Only RefSeq should be here!!!
                                             normalize=True,
                                             replace_reference=True
                                             )    
    
    # Create dictionary to store the information    
    hgvs_transcript_to_hgvs_protein = va_func.myc_to_p(hgvs_transcript, evm, hdp, hp, rhn, vm, sf, re_to_p=False)
    hgvs_transcript_to_hgvs_protein = hgvs_transcript_to_hgvs_protein['hgvs_protein']
    
    return hgvs_transcript_to_hgvs_protein

"""
Return all aligned transcripts for a given genomic hgvs object
"""


def fetch_aligned_transcripts(hgvs_genomic, transcript_model):
    tx_list = []
    
    if transcript_model == 'ensembl' or transcript_model == 'all':
        enst_list = hdp.get_tx_for_region(hgvs_genomic.ac, 'genebuild', hgvs_genomic.posedit.pos.start.base - 1,
                                          hgvs_genomic.posedit.pos.end.base)
        
        # Remove transcripts with incomplete identifiers
        cp_enst_list = []
        for et in enst_list:
            if re.search('\d+\.\d+', et[0]):
                cp_enst_list.append(et)
    
        tx_list = tx_list + cp_enst_list
                          
    if transcript_model == 'refseq' or transcript_model == 'all':   
        refseq_list = hdp.get_tx_for_region(hgvs_genomic.ac, 'splign', hgvs_genomic.posedit.pos.start.base - 1,
                                            hgvs_genomic.posedit.pos.end.base)
        tx_list = tx_list + refseq_list
        
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


def gap_checker(hgvs_transcript, hgvs_genomic, un_norm_hgvs_genomic, genome_build):
    
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
                                        genome_build, hdp, hp, sf, hgvs_version)

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
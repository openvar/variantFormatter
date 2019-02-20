# -*- coding: utf-8 -*-

# import modules
import formatter
import re
import collections
import VariantFormatter
 
# Custom Exceptions
class vcf2hgvsError(Exception):
    pass
class hgvs2VcfError(Exception):
    pass
class variableError(Exception):
    pass

    
# Create variantformatter object
class GenomicDescriptions(object):
    """
    Object contains genomic level sequence variant descriptions in the pseudo vcf (p_vcf)
    genomic hgvs (g_hgvs) and un-normalized g_hgvs. The reference bases are the HGVS
    description reference nucleotide sequence corresponding to the specified range 
    """
    # Initialise and add initialisation data to the object
    def __init__(self, p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases):
        self.p_vcf = p_vcf
        removed_ref_hgvs_g = formatter.remove_reference(g_hgvs)
        self.g_hgvs = removed_ref_hgvs_g
        self.un_norm_hgvs = un_norm_hgvs
        self.g_hgvs_ref = hgvs_ref_bases

    
# Create variantformatter object
class FormatVariant(object):
    """
    Returns and object containing the:
    Configuration metadata i.e. input variant_description, genome_build, transcript_model
    and specified transcripts;
    GenomicDescriptions object;
    transcript and protein level descriptions - includes gap checking!    
    """
    # Initialise and add initialisation data to the object
    def __init__(self, variant_description, genome_build, transcript_model=None, specify_transcripts=None):
        
        self.variant_description = variant_description
        
        if genome_build not in ['GRCh37', 'GRCh38']:
            raise variableError("genome_build must be one of: 'GRCh37'; 'GRCh38'")
        else:
            self.genome_build = genome_build
        
        if transcript_model is None:
            transcript_model = 'all' 
        elif transcript_model not in ['ensembl', 'refseq', 'all']:
            raise variableError("transcript_model must be one of: 'ensembl'; 'refseq'; 'all'")
        else:
            self.transcript_model = transcript_model
        
        self.specify_transcripts = specify_transcripts
                
        # vcf2hgvs route
        if re.match('N[CTW]_', self.variant_description):
            try:
                hgvs_genomic = formatter.parse(self.variant_description)
                vcf_dictionary = formatter.hgvs_genomic2vcf(hgvs_genomic, self.genome_build)
                vcf_list = [vcf_dictionary['grc_chr'], vcf_dictionary['pos'], vcf_dictionary['ref'], vcf_dictionary['alt']]
                p_vcf = ':'.join(vcf_list)
            except Exception as e:
                raise hgvs2VcfError(str(e)) 
            try:
                genomic_level = formatter.vcf2hgvs_genomic(p_vcf, self.genome_build)
            except Exception as e:
                raise vcf2hgvsError(str(e))
            else:
                if genomic_level['error'] != '':
                    raise vcf2hgvsError(str(genomic_level['error']))
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']
            
        # hgvs2vcf route
        elif re.match('chr\d+\-', self.variant_description) or re.match('chr\d+:', self.variant_description) or re.match('\d+\-', self.variant_description)  or re.match('\d+:', self.variant_description):        
            try:
                genomic_level = formatter.vcf2hgvs_genomic(self.variant_description, self.genome_build)
            except Exception as e:
                raise vcf2hgvsError(str(e))
            else:
                if genomic_level['error'] != '':
                    raise vcf2hgvsError(str(genomic_level['error']))
                p_vcf = self.variant_description
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']
        
        # Not recognidsed
        else:
            raise variableError('Variant description ' + self.variant_description + ' is not in a supported format') 
        
        # Create genomic_descriptions object
        gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases)
        self.genomic_descriptions = gds
        
        
        # Add transcript and protein data
        prelim_transcript_descriptions = {}
        transcript_list = []
        
        # Transcripts specified
        if self.specify_transcripts is not None:
            transcript_list = str(self.specify_transcripts).split('|')
        
        # No transcripts specified
        else:
            transcript_list = formatter.fetch_aligned_transcripts(g_hgvs, self.transcript_model)
            
        # Create transcript level descriptions
        for tx_alignment_data in transcript_list:
            tx_id = tx_alignment_data[0]
            hgvs_transcript_dict = formatter.hgvs_genomic2hgvs_transcript(g_hgvs, tx_id)
           
            # Gap checking            
            try:
                am_i_gapped = formatter.gap_checker(hgvs_transcript_dict['hgvs_transcript'], g_hgvs, un_norm_hgvs, self.genome_build)
            except Exception as e:

                am_i_gapped = {'hgvs_transcript': None,
                                'position_lock': False,
                                'gapped_alignment_warning': None,
                                'corrective_action': None,
                                'gap_position':None,
                                'transcript_accession': tx_id,
                                'error': hgvs_transcript_dict['error']  
                                }               
                
                
                # add to dictionary
                am_i_gapped['hgvs_protein_tlc'] = None
                am_i_gapped['hgvs_protein_slc'] = None                          
            
            else:   
                                                        
                # map to Protein
                if am_i_gapped['hgvs_transcript'].type == 'c':
                    hgvs_protein_tlc = formatter.hgvs_transcript2hgvs_protein(am_i_gapped['hgvs_transcript'], self.genome_build)
                    hgvs_protein_slc = formatter.single_letter_protein(hgvs_protein_tlc)
                if am_i_gapped['hgvs_transcript'].type == 'n':
                    hgvs_protein_tlc = 'non-coding'
                    hgvs_protein_slc = 'non-coding'                 
            
                # add to dictionary
                am_i_gapped['hgvs_protein_tlc'] = str(hgvs_protein_tlc)
                am_i_gapped['hgvs_protein_slc'] = str(hgvs_protein_slc)
                am_i_gapped['error'] = hgvs_transcript_dict['error']
                        
                # Remove ref bases
                removed_ref_tx = formatter.remove_reference(am_i_gapped['hgvs_transcript'])
                am_i_gapped['hgvs_transcript'] = str(removed_ref_tx)
            
            # Order the tx_p output
            order_my_tp = collections.OrderedDict()
            order_my_tp['t_hgvs'] = am_i_gapped['hgvs_transcript']
            order_my_tp['t_hgvs_ref'] = hgvs_transcript_dict['ref_bases']
            order_my_tp['p_hgvs_tlc'] = am_i_gapped['hgvs_protein_tlc']
            order_my_tp['p_hgvs_slc'] = am_i_gapped['hgvs_protein_slc']
            order_my_tp['gapped_alignment_warning'] = am_i_gapped['gapped_alignment_warning']
            order_my_tp['gap_position'] = am_i_gapped['gap_position']
            order_my_tp['corrective_action'] = am_i_gapped['corrective_action']            
            order_my_tp['position_lock'] = am_i_gapped['position_lock']

            
            # add to output dictionary keyed by tx_ac
            prelim_transcript_descriptions[tx_id] = order_my_tp
            
        self.t_and_p_descriptions = prelim_transcript_descriptions
    
    # Create ordered output
    def stucture_data(self):        
    	bring_order = collections.OrderedDict()
    	
    	# Add the data to the ordered dictionary structure
    	bring_order['submitted_variant'] = str(self.variant_description)
    	bring_order['p_vcf'] = str(self.genomic_descriptions.p_vcf)
    	bring_order['g_hgvs'] = str(self.genomic_descriptions.g_hgvs)
    	bring_order['g_hgvs_ref'] = str(self.genomic_descriptions.g_hgvs_ref)
    	bring_order['hgvs_t_and_p'] = self.t_and_p_descriptions
    	return bring_order
    	
    def collect_metadata(self):	
    	meta = collections.OrderedDict()
    	meta['api_version'] = VariantFormatter.__version__
    	meta['api_released'] = VariantFormatter.__released__
    	meta['hgvs_version'] = VariantFormatter.__hgvs_version__
    	meta['uta_schema'] = VariantFormatter.__uta_schema__
    	meta['seqrepo_db'] = VariantFormatter.__seqrepo_db__
    	return meta
                
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
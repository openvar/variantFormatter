# -*- coding: utf-8 -*-

"""
This module creates an initialization object.
This object connects to the hgvs Python library and soosciated databses

The Initialization object is used by FormatVariant
The FormatVariant object contains all HGVS descriptions available for a given genomic variant, g_to_p
"""

# import modules
import os
import re
import VariantFormatter
import formatter
import collections
from configparser import ConfigParser, RawConfigParser
import hgvs
import hgvs.dataproviders.uta
 
# Custom Exceptions
class vcf2hgvsError(Exception):
    pass
class hgvs2VcfError(Exception):
    pass
class variableError(Exception):
    pass


# Configure databases  
class initializeFormatter(object):
    """
    Sets paths to the UTA and SeqRepo databses
    """
    def __init__(self):

        # First load from the configuration file, if it exists.
        configName=".VariantFormatter.conf"
        configRoot=os.environ.get('HOME')

        # Now configpath points to the config file itself.
        configPath=os.path.join(configRoot,configName)
        # Does the file exist?
        if not os.path.exists(configPath):
            self.createConfig(configPath)

        # Load the configuration file.
        config=RawConfigParser(allow_no_value=True)
        with open(configPath) as file:
            config.read_file(file)

        # Set up versions
        __version__ = VariantFormatter.__version__
        self.version=__version__
        self.released_version = VariantFormatter._is_released_version

        # Handle databases
        if config["seqrepo"]["location"]!=None:
            self.seqrepoPath=config["seqrepo"]["location"]
            self.seqrepoVersion=str(self.seqrepoPath).split('/')[-1]
            os.environ['HGVS_SEQREPO_DIR']=self.seqrepoPath
        else:
            valueerror = "\nThe seqrepo database location has not been set in %s" % configPath
            example = "example: location = %s/seqrepo/2018-08-21\n" % configRoot
            valueerror = valueerror + '\n' + example
            raise ValueError(valueerror)

        if config["uta"]["location"]!=None:
            self.utaPath=config["uta"]["location"]
            self.utaVersion=str(self.utaPath).split('/')[-1]
            os.environ['UTA_DB_URL']=self.utaPath
        else:
            valueerror = "\nThe uta database location has not been set in %s" % configPath
            example = "example: location = postgresql://uta_admin:uta_admin@127.0.0.1/uta/uta_20180821\n"
            valueerror = valueerror + '\n' + example
            raise ValueError(valueerror)
        
        # Set up HGVS
        # Configure hgvs package global settings
        self.hgvsVersion=hgvs.__version__
        hgvs.global_config.uta.pool_max = 25
        hgvs.global_config.formatting.max_ref_length = 1000000
        self.hdp = hgvs.dataproviders.uta.connect(pooling=True)
        self.hp = hgvs.parser.Parser()
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=True)
        self.vr = hgvs.validator.Validator(self.hdp)
        self.sf = hgvs.dataproviders.seqfetcher.SeqFetcher()
        self.splign_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                               cross_boundaries=False,
                                               shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                               alt_aln_method='splign'
                                               )
                                               
        self.genebuild_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                  cross_boundaries=False,
                                                  shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                                  alt_aln_method='genebuild'
                                                  )
                                                  
        self.reverse_splign_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                       cross_boundaries=False,
                                                       shuffle_direction=5,
                                                       alt_aln_method='splign'
                                                       )
                                                       
        self.reverse_genebuild_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                          cross_boundaries=False,
                                                          shuffle_direction=5,
                                                          alt_aln_method='genebuild'
                                                          )
                                                          
    def myConfig(self):
        '''
        #Returns configuration:
        #version, hgvs version, uta schema, seqrepo db.
        '''
        meta = collections.OrderedDict()
        meta['VariantFormatter_version'] = self.version
        meta['biocommons.hgvs_version'] = self.hgvsVersion
        meta['biocommons.uta_schema'] = self.utaVersion
        meta['biocommons.seqrepo_db'] = self.seqrepoVersion
        return meta

    def createConfig(self,outPath):
        '''
        # This function reads from the default configuration file stored in the same folder as this module,
        # and transfers it to outPath.
        # Outpath should include a filename.
        '''
        lines=[]
        inPath=os.path.join(os.path.dirname(os.path.realpath(__file__)),"configuration/config.ini")
        with open(inPath) as file:
            for l in file:
                lines.append(l)
        with open(outPath, "w") as file:
            for l in lines:
                file.write(l)                                                                 
                        

# Create variantformatter object
class GenomicDescriptions(object):
    """
    Object contains genomic level sequence variant descriptions in the pseudo vcf (p_vcf)
    genomic hgvs (g_hgvs) and un-normalized g_hgvs. The reference bases are the HGVS
    description reference nucleotide sequence corresponding to the specified range 
    """
    # Initialise and add initialisation data to the object
    def __init__(self, p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error):
        if p_vcf == "None":
            p_vcf = None
        if g_hgvs == "None":
            g_hgvs = None
        if un_norm_hgvs == "None":
            un_norm_hgvs = None
        if hgvs_ref_bases == "None":
            hgvs_ref_bases = None
        if gen_error == "None":
            gen_error = None

        # Create object
        self.p_vcf = p_vcf
        try:
            removed_ref_hgvs_g = formatter.remove_reference(g_hgvs)
        except AttributeError:
            removed_ref_hgvs_g = None
        self.g_hgvs = removed_ref_hgvs_g
        self.un_norm_hgvs = un_norm_hgvs
        self.g_hgvs_ref = hgvs_ref_bases
        self.gen_error = gen_error

    

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
    def __init__(self, variant_description, genome_build, vfo, transcript_model=None, specify_transcripts=None, checkOnly=False):
        
        self.variant_description = variant_description
        self.varForm = vfo
        
        if genome_build not in ['GRCh37', 'GRCh38']:
            p_vcf = None
            g_hgvs = None
            hgvs_ref_bases = None
            un_norm_hgvs = None
            gen_error = "genome_build must be one of: 'GRCh37'; 'GRCh38'"
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error)
            self.genomic_descriptions = gds
            return
        else:
            self.genome_build = genome_build
        
        if transcript_model is None:
            transcript_model = 'all' 
        elif transcript_model not in ['ensembl', 'refseq', 'all']:
            p_vcf = None
            g_hgvs = None
            hgvs_ref_bases = None
            un_norm_hgvs = None
            gen_error = "transcript_model must be one of: 'ensembl'; 'refseq'; 'all'"
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error)
            self.genomic_descriptions = gds
            return
        else:
            self.transcript_model = transcript_model
        
        self.specify_transcripts = specify_transcripts
                
        # hgvs2vcf route
        if re.match('N[CTW]_', self.variant_description):
            try:
                hgvs_genomic = formatter.parse(self.variant_description, self.varForm)
                vcf_dictionary = formatter.hgvs_genomic2vcf(hgvs_genomic, self.genome_build, self.varForm)
                vcf_list = [vcf_dictionary['grc_chr'], vcf_dictionary['pos'], vcf_dictionary['ref'], vcf_dictionary['alt']]
                p_vcf = ':'.join(vcf_list)
            except Exception as e:
                raise hgvs2VcfError(str(e))
            try:
                genomic_level = formatter.vcf2hgvs_genomic(p_vcf, self.genome_build, self.varForm)
            except Exception as e:
                raise vcf2hgvsError(str(e))
            else:
                if genomic_level['error'] != '':
                    p_vcf = None
                    g_hgvs = None
                    hgvs_ref_bases = None
                    un_norm_hgvs = None
                    gen_error = genomic_level['error']
                    gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error)
                    self.genomic_descriptions = gds
                    return
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']
            
        # vcf2hgvs route
        elif re.match('chr[\w\d]+\-', self.variant_description) or re.match('chr[\w\d]+:', self.variant_description) or re.match('[\w\d]+\-', self.variant_description)  or re.match('[\w\d]+:', self.variant_description):
            try:
                genomic_level = formatter.vcf2hgvs_genomic(self.variant_description, self.genome_build, self.varForm)
            except Exception as e:
                raise vcf2hgvsError(str(e))
                return
            else:
                if genomic_level['error'] != '':
                    p_vcf = None
                    g_hgvs = None
                    hgvs_ref_bases = None
                    un_norm_hgvs = None
                    gen_error = genomic_level['error']
                    gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error)
                    self.genomic_descriptions = gds
                    return
                p_vcf = self.variant_description
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']
        
        # Not recognidsed
        else:
            raise variableError('Variant description ' + self.variant_description + ' is not in a supported format') 
        
        # Create genomic_descriptions object
        gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error=None)
        self.genomic_descriptions = gds

        # Return on checkonly
        if checkOnly is True:
            return

        # Add transcript and protein data
        prelim_transcript_descriptions = {}
        transcript_list = []
        
        # Transcripts specified        
        if self.specify_transcripts is not None:
            trans_list = str(self.specify_transcripts).split('|')
            for tx in trans_list:
                transcript_list.append([tx,''])

        # No transcripts specified
        else:
            transcript_list = formatter.fetch_aligned_transcripts(g_hgvs, self.transcript_model, self.varForm)
            
        # Create transcript level descriptions
        for tx_alignment_data in transcript_list:
            tx_id = tx_alignment_data[0]
            hgvs_transcript_dict = formatter.hgvs_genomic2hgvs_transcript(g_hgvs, tx_id, self.varForm)
            
            # Gap checking            
            try:
                am_i_gapped = formatter.gap_checker(hgvs_transcript_dict['hgvs_transcript'], g_hgvs, un_norm_hgvs, self.genome_build, self.varForm)
            except Exception as e:
                if hgvs_transcript_dict['error'] == '':
                    hgvs_transcript_dict['error'] = None

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
                if hgvs_transcript_dict['error'] == '':
                    hgvs_transcript_dict['error'] = None

                # map to Protein
                if am_i_gapped['hgvs_transcript'].type == 'c':
                    hgvs_protein_tlc = formatter.hgvs_transcript2hgvs_protein(am_i_gapped['hgvs_transcript'], self.genome_build, self.varForm)
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
            if am_i_gapped['error'] == '':
                am_i_gapped['error'] = None
            order_my_tp = collections.OrderedDict()
            order_my_tp['t_hgvs'] = am_i_gapped['hgvs_transcript']
            order_my_tp['p_hgvs_tlc'] = am_i_gapped['hgvs_protein_tlc']
            order_my_tp['p_hgvs_slc'] = am_i_gapped['hgvs_protein_slc']
            order_my_tp['gapped_alignment_warning'] = am_i_gapped['gapped_alignment_warning']
            order_my_tp['gap_statement'] = am_i_gapped['gap_position']
            order_my_tp['transcript_variant_error'] = am_i_gapped['error']
            
            # add to output dictionary keyed by tx_ac
            prelim_transcript_descriptions[tx_id] = order_my_tp
            
        self.t_and_p_descriptions = prelim_transcript_descriptions
    
    # Create ordered output
    def stucture_data(self):        
        bring_order = collections.OrderedDict()
        
        # Add the data to the ordered dictionary structure
        bring_order['p_vcf'] = self.genomic_descriptions.p_vcf
        bring_order['g_hgvs'] = self.genomic_descriptions.g_hgvs # Is the removed ref version!
        bring_order['genomic_variant_error'] = self.genomic_descriptions.gen_error
        try:
            if self.t_and_p_descriptions == {}:
                bring_order['hgvs_t_and_p'] = None
            else:
                bring_order['hgvs_t_and_p'] = self.t_and_p_descriptions
        except AttributeError:
            bring_order['hgvs_t_and_p'] = None
        brought_order = {str(self.variant_description): bring_order}
        return brought_order
        
    def collect_metadata(self):
        meta = collections.OrderedDict()
        meta['api_version'] = self.varForm.version
        meta['hgvs_version'] = self.varForm.hgvsVersion
        meta['uta_schema'] = self.varForm.utaVersion
        meta['seqrepo_db'] = self.varForm.seqrepoVersion
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
# -*- coding: utf-8 -*-

"""
This module creates an initialization object.
This object connects to the hgvs Python library and soosciated databses

The Initialization object is used by FormatVariant
The FormatVariant object contains all HGVS descriptions available for a given genomic variant, g_to_p
"""

# import modules
import re
import collections
import copy
import vvhgvs.exceptions
import VariantFormatter.formatter as formatter
import VariantValidator.modules.liftover as lo


# Custom Exceptions
class vcf2hgvsError(Exception):
    pass


class hgvs2VcfError(Exception):
    pass


class variableError(Exception):
    pass


# Create Genomic Descriptions object
class GenomicDescriptions(object):
    """
    Object contains genomic level sequence variant descriptions in the pseudo vcf (p_vcf)
    genomic hgvs (g_hgvs) and un-normalized g_hgvs. The reference bases are the HGVS
    description reference nucleotide sequence corresponding to the specified range 
    """
    # Initialise and add initialisation data to the object
    def __init__(self, p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build, variant_description):
        if p_vcf == "None":
            p_vcf = None
        try:
            if g_hgvs == "None":
                g_hgvs = None
            elif ('NC_012920.1' in g_hgvs.ac or 'NC_001807.4' in g_hgvs.ac) and ":g." in variant_description:
                gen_error = "The given reference sequence (%s) does not match the DNA type (g). For %s, " \
                            "please use (m). " \
                            "For g. variants, please use a linear genomic reference sequence" % (g_hgvs.ac, g_hgvs.ac)
        except AttributeError:
            if g_hgvs == "None":
                g_hgvs = None
        if un_norm_hgvs == "None":
            un_norm_hgvs = None
        if hgvs_ref_bases == "None":
            hgvs_ref_bases = None
        if gen_error == "None":
            gen_error = None

        # Warn incorrect m. accession for hg19
        try:
            if ("NC_012920.1" in str(g_hgvs)) and "hg19" in genome_build:
                gen_error = "NC_012920.1 is not associated with genome build hg19, instead use genome build GRCh37"
            elif ("NC_001807.4" in str(g_hgvs)) and "GRCh37" in genome_build:
                gen_error = "NC_001807.4 is not associated with genome build GRCh37, instead use genome build hg19"
        except TypeError:
            pass

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
        self.selected_build = genome_build

    
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
    def __init__(self, variant_description, genome_build, vfo, transcript_model=None, specify_transcripts=None,
                 checkOnly=False, liftover=False):

        self.variant_description = variant_description
        self.vfo = vfo
        # Add warning level
        self.warning_level = None
        gen_error = None
        self.liftover = liftover
        
        if genome_build not in ['GRCh37', 'GRCh38', 'hg19', 'hg38']:
            p_vcf = None
            g_hgvs = None
            hgvs_ref_bases = None
            un_norm_hgvs = None
            gen_error = "genome_build must be one of: 'GRCh37'; 'GRCh38'; 'hg19'; 'hg38'"
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                      variant_description)
            self.genomic_descriptions = gds
            self.warning_level = 'genomic_variant_warning'
            return
        else:
            self.genome_build = genome_build

        if transcript_model is None:
            transcript_model = 'all' 
        if transcript_model not in ['ensembl', 'refseq', 'all']:
            p_vcf = None
            g_hgvs = None
            hgvs_ref_bases = None
            un_norm_hgvs = None
            gen_error = "transcript_model must be one of: 'ensembl'; 'refseq'; 'all'"
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                      variant_description)
            self.genomic_descriptions = gds
            self.warning_level = 'genomic_variant_warning'
            return
        else:
            self.transcript_model = transcript_model
        
        self.specify_transcripts = specify_transcripts
                
        # hgvs2vcf route
        if re.match('N[CTW]_', self.variant_description):
            try:
                hgvs_genomic = formatter.parse(self.variant_description, self.vfo)
                vfo.vr.validate(hgvs_genomic)
            except Exception as e:
                validation = vfo.validate(self.variant_description, self.genome_build, 'all',
                                          liftover_level=None).format_as_dict(test=True)

                # Can the variant be auto-corrected
                # if "warning" not in validation["flag"]:
                reset_variant = None
                edit_warnings = None
                # Code to auto-recover dud description where possible
                for val_key, val_val in validation.items():

                    try:
                        if "primary_assembly_loci" in val_val.keys():
                            reset_variant = val_val["primary_assembly_loci"][genome_build.lower()
                                                                             ]["hgvs_genomic_description"]
                            validation_warned = val_val["validation_warnings"]
                            edit_warnings = []
                            for wrn in validation_warned:
                                if "automapped to" in wrn:
                                    edit_warnings.append(wrn)

                            break
                        else:
                            continue
                    except AttributeError:
                        continue
                    except KeyError:
                        validation_warned = val_val["validation_warnings"]
                        edit_warnings = ", ".join(validation_warned)
                        p_vcf = None
                        g_hgvs = None
                        hgvs_ref_bases = None
                        un_norm_hgvs = None
                        gen_error = edit_warnings
                        gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                                  variant_description)
                        self.genomic_descriptions = gds
                        self.warning_level = 'genomic_variant_warning'
                        return

                hgvs_genomic = formatter.parse(reset_variant, self.vfo)
                try:
                    recovery_error = edit_warnings[0]
                except IndexError:
                    recovery_error = None
                self.warning_level = 'genomic_variant_warning'

            # Continuation - No exception
            try:
                vcf_dictionary = formatter.hgvs_genomic2vcf(hgvs_genomic, self.genome_build, self.vfo)
                if vcf_dictionary['grc_chr'] == "NC_001807.4" and genome_build == "hg19":
                    chr_num = vcf_dictionary['ucsc_chr']
                else:
                    chr_num = vcf_dictionary['grc_chr']
                vcf_list = [chr_num, vcf_dictionary['pos'], vcf_dictionary['ref'],
                            vcf_dictionary['alt']]

                p_vcf = ':'.join(vcf_list)
            except Exception as e:
                if "Variant span is outside sequence bounds" in str(e):
                    e = "The specified coordinate is outside the boundaries of reference sequence %s" % self.\
                        variant_description.split(':')[0]
                p_vcf = None
                g_hgvs = None
                hgvs_ref_bases = None
                un_norm_hgvs = None
                gen_error = str(e)
                gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                          variant_description)
                self.genomic_descriptions = gds
                self.warning_level = 'genomic_variant_warning'
                return

            try:
                genomic_level = formatter.vcf2hgvs_genomic(p_vcf, self.genome_build, self.vfo)
            except Exception as e:
                if "Variant span is outside sequence bounds" in str(e):
                    e = "The specified coordinate is outside the boundaries of reference sequence %s" % self.\
                        variant_description.split(':')[0]
                p_vcf = None
                g_hgvs = None
                hgvs_ref_bases = None
                un_norm_hgvs = None
                gen_error = str(e)
                gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                          variant_description)
                self.genomic_descriptions = gds
                self.warning_level = 'genomic_variant_warning'
                return
            else:
                if genomic_level['error'] != '':
                    p_vcf = None
                    g_hgvs = None
                    hgvs_ref_bases = None
                    un_norm_hgvs = None
                    gen_error = genomic_level['error']
                    gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                              variant_description)
                    self.genomic_descriptions = gds
                    self.warning_level = 'genomic_variant_warning'
                    return
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']

                # Check for auto-update of variant description
                if str(g_hgvs.posedit.pos) not in str(self.variant_description):
                    self.warning_level = 'genomic_variant_warning'
                    gen_error = self.variant_description + " updated to " + str(formatter.remove_reference(g_hgvs))
                else:
                    gen_error = None
            
        # vcf2hgvs route
        elif re.match('chr[\w\d]+\-', self.variant_description) or re.match(
                'chr[\w\d]+:', self.variant_description) or re.match('[\w\d]+\-', self.variant_description) \
                or re.match('[\w\d]+:', self.variant_description):
            try:
                genomic_level = formatter.vcf2hgvs_genomic(self.variant_description, self.genome_build, self.vfo)
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
                    gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                              variant_description)
                    self.genomic_descriptions = gds
                    self.warning_level = 'genomic_variant_warning'
                    return
                p_vcf = self.variant_description
                g_hgvs = genomic_level['hgvs_genomic']
                un_norm_hgvs = genomic_level['un_normalized_hgvs_genomic']
                hgvs_ref_bases = genomic_level['ref_bases']
        
        # Not recognised
        else:
            p_vcf = None
            g_hgvs = None
            hgvs_ref_bases = None
            un_norm_hgvs = None
            gen_error = 'Variant description ' + self.variant_description + ' is not in a supported format'
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error, genome_build,
                                      variant_description)
            self.genomic_descriptions = gds
            self.warning_level = 'submission_warning'
            return

        # Create genomic_descriptions object
        try:
            if recovery_error is not None:
                gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, recovery_error,
                                          genome_build, variant_description)
            else:
                gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error,
                                          genome_build, variant_description)
        except UnboundLocalError:
            gds = GenomicDescriptions(p_vcf, g_hgvs, un_norm_hgvs, hgvs_ref_bases, gen_error,
                                      genome_build, variant_description)
        self.genomic_descriptions = gds

        # Return on checkonly
        if checkOnly is True:
            return

        # Add transcript and protein data
        prelim_transcript_descriptions = {}
        transcript_list = []

        # Requires hgvs mapping so update genome build if hg
        if self.genome_build == 'hg19':
            self.genome_build = 'GRCh37'
        if self.genome_build == 'hg38':
            self.genome_build = 'GRCh38'

        # Transcripts specified
        if self.specify_transcripts is not None:
            trans_list = str(self.specify_transcripts).split('|')
            for tx in trans_list:
                transcript_list.append([tx,''])

        # No transcripts specified
        else:
            transcript_list = formatter.fetch_aligned_transcripts(g_hgvs, self.transcript_model, self.vfo)

        # Remove malform IDs
        cp_transcript_list = copy.copy(transcript_list)
        transcript_list = []
        for tx in cp_transcript_list:
            # Known UTA ID malforms
            if re.search('\/', tx[0]):
                continue
            else:
                transcript_list.append(tx)

        # Create a variable to trap direct g_g liftover
        g_to_g_lift = {}

        # Create transcript level descriptions
        for tx_alignment_data in transcript_list:
            tx_id = tx_alignment_data[0]
            hgvs_transcript_dict = formatter.hgvs_genomic2hgvs_transcript(g_hgvs, tx_id, self.vfo)
            
            # Gap checking
            try:
                am_i_gapped = formatter.gap_checker(hgvs_transcript_dict['hgvs_transcript'], g_hgvs,
                                                    self.genome_build, self.vfo)
            except Exception:
                self.warning_level = 'processing_error'
                if hgvs_transcript_dict['error'] == '':
                    hgvs_transcript_dict['error'] = None

                am_i_gapped = {'hgvs_transcript': None, 'position_lock': False, 'gapped_alignment_warning': None,
                               'corrective_action': None, 'gap_position': None, 'transcript_accession': tx_id,
                               'error': hgvs_transcript_dict['error'], 'hgvs_protein_tlc': None,
                               'hgvs_protein_slc': None}

            # add to dictionary
            else:
                if hgvs_transcript_dict['error'] == '':
                    hgvs_transcript_dict['error'] = None

                # Adds LOVD requested to transcript level only validation to checkOnly. Will not produce p.
                if "tx" in str(checkOnly):
                    hgvs_protein_tlc = None
                    hgvs_protein_slc = None

                # Standard checkOnly will produce protein
                else:
                    # map to Protein
                    if am_i_gapped['hgvs_transcript'].type == 'c':
                        try:
                            hgvs_protein_tlc = formatter.hgvs_transcript2hgvs_protein(am_i_gapped['hgvs_transcript'],
                                                                                      self.genome_build,
                                                                                      self.vfo)
                            # Handle edits that have been stringified
                            try:
                                hgvs_protein_tlc.posedit.edit.ref
                            except AttributeError:
                                hgvs_protein_tlc = formatter.parse(str(hgvs_protein_tlc), self.vfo)

                            hgvs_protein_slc = formatter.single_letter_protein(hgvs_protein_tlc)
                        except NotImplementedError as e:
                            hgvs_protein_tlc = None
                            hgvs_protein_slc = None
                            hgvs_transcript_dict['error'] = str(e)
                        except vvhgvs.exceptions.HGVSDataNotAvailableError as e:
                            hgvs_protein_tlc = None
                            hgvs_protein_slc = None
                            hgvs_transcript_dict['error'] = str(e)

                    if am_i_gapped['hgvs_transcript'].type == 'n':
                        hgvs_protein_tlc = None
                        hgvs_protein_slc = None
            
                # add to dictionary
                if hgvs_protein_tlc is not None:
                    am_i_gapped['hgvs_protein_tlc'] = str(hgvs_protein_tlc)
                    am_i_gapped['hgvs_protein_slc'] = str(hgvs_protein_slc)
                else:
                    am_i_gapped['hgvs_protein_tlc'] = hgvs_protein_tlc
                    am_i_gapped['hgvs_protein_slc'] = hgvs_protein_slc
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

            # Add liftover information if requested
            if self.liftover is not False:
                if self.genomic_descriptions.selected_build == 'hg19' or self.genomic_descriptions.selected_build \
                        == 'GRCh37':
                    build_to = 'GRCh38'
                elif self.genomic_descriptions.selected_build == 'hg38' or self.genomic_descriptions.selected_build \
                        == 'GRCh38':
                    build_to = 'GRCh37'

                # Look for previous g_to_g lift
                if (order_my_tp['transcript_variant_error'] is not None and g_to_g_lift == {}) or (
                        order_my_tp['transcript_variant_error'] is None):

                    if order_my_tp['t_hgvs'] is not None:
                        specified_tx_variant = formatter.parse(order_my_tp['t_hgvs'], self.vfo)
                    else:
                        specified_tx_variant = None

                    current_lift = lo.liftover(self.genomic_descriptions.g_hgvs,
                                               self.genomic_descriptions.selected_build,
                                               build_to,
                                               vfo.splign_normalizer,
                                               vfo.reverse_splign_normalizer,
                                               None,
                                               vfo,
                                               specify_tx=tx_id,
                                               liftover_level=self.liftover,
                                               gap_map=formatter.gap_checker,
                                               vfo=self.vfo,
                                               specified_tx_variant=specified_tx_variant
                                               )

                    if "am_i_gapped" in current_lift.keys():
                        if order_my_tp['gapped_alignment_warning'] is "":
                            order_my_tp['gapped_alignment_warning'] = current_lift['am_i_gapped'][
                                'gapped_alignment_warning']
                        if order_my_tp['gap_statement'] is "":
                            order_my_tp['gap_statement'] = current_lift['am_i_gapped']['gap_position']
                        current_lift.pop("am_i_gapped")

                    if g_to_g_lift == {}:
                        g_to_g_lift = current_lift

                elif order_my_tp['transcript_variant_error'] is not None and g_to_g_lift != {}:
                    current_lift = g_to_g_lift

                # Copy the liftover and split into primary and alt
                cp_current_lift = copy.deepcopy(current_lift)
                scaff_lift = copy.deepcopy(current_lift)
                alt_list = []
                for key, val in current_lift.items():
                    for chr_type in val.keys():
                        if 'NC_' not in chr_type:
                            del cp_current_lift[key][chr_type]
                        if '24.' in chr_type:
                            del cp_current_lift[key][chr_type]
                for key, val in scaff_lift.items():
                    for chr_type in val.keys():
                        if 'NC_' in chr_type and "24." not in chr_type:
                            continue
                        else:
                            for compile_list_key, compile_list_val in val.items():
                                if "NC" in compile_list_key:
                                    if "24." in compile_list_key:
                                        alt_list.append({key: compile_list_val})
                                    else:
                                        continue
                                else:
                                    alt_list.append({key: compile_list_val})

                order_my_tp['primary_assembly_loci'] = cp_current_lift
                order_my_tp['alt_genomic_loci'] = alt_list

            # add to output dictionary keyed by tx_ac
            prelim_transcript_descriptions[tx_id] = order_my_tp
            
        self.t_and_p_descriptions = prelim_transcript_descriptions
    
    # Create ordered output
    def stucture_data(self):        
        bring_order = collections.OrderedDict()
        
        # Add the data to the ordered dictionary structure
        bring_order['p_vcf'] = self.genomic_descriptions.p_vcf
        if self.genomic_descriptions.g_hgvs is not None:
            # Handle mitochondrial genomic variants
            if 'NC_012920.1' in self.genomic_descriptions.g_hgvs or 'NC_001807.4' in self.genomic_descriptions.g_hgvs:
                self.genomic_descriptions.g_hgvs = self.genomic_descriptions.g_hgvs.replace(':g.', ':m.')
                try:
                    self.genomic_descriptions.gen_error = self.genomic_descriptions.gen_error.replace(":g.", ":m.")
                except AttributeError:
                    pass
        bring_order['g_hgvs'] = self.genomic_descriptions.g_hgvs # Is the removed ref version!
        bring_order['selected_build'] = self.genomic_descriptions.selected_build
        bring_order['genomic_variant_error'] = self.genomic_descriptions.gen_error
        try:
            if self.t_and_p_descriptions == {}:
                bring_order['hgvs_t_and_p'] = {'intergenic': {'primary_assembly_loci': None}}
                bring_order['hgvs_t_and_p'] = {'intergenic': {'alt_genomic_loci': None}}
                if self.liftover is not False:
                    if self.genomic_descriptions.selected_build == 'hg19' or self.genomic_descriptions.selected_build \
                            == 'GRCh37':
                        build_to = 'GRCh38'
                    elif self.genomic_descriptions.selected_build == 'hg38' or self.genomic_descriptions.selected_build \
                            == 'GRCh38':
                        build_to = 'GRCh37'
                    current_lift = lo.liftover(self.genomic_descriptions.g_hgvs,
                                               self.genomic_descriptions.selected_build,
                                               build_to,
                                               self.vfo.splign_normalizer,
                                               self.vfo.reverse_splign_normalizer,
                                               None,
                                               self.vfo,
                                               specify_tx=False,
                                               liftover_level=self.liftover
                                               )

                    # Copy the liftover and split into primary and alt
                    cp_current_lift = copy.deepcopy(current_lift)
                    scaff_lift = copy.deepcopy(current_lift)
                    cp_scaff_lift = copy.deepcopy(current_lift)
                    for key, val in current_lift.items():
                        for chr_type in val.keys():
                            if not 'NC_' in chr_type:
                                del cp_current_lift[key][chr_type]
                    for key, val in scaff_lift.items():
                        for chr_type in val.keys():
                            if 'NC_' in chr_type:
                                del cp_scaff_lift[key][chr_type]

                    bring_order['hgvs_t_and_p']['intergenic']['primary_assembly_loci'] = cp_current_lift
                    bring_order['hgvs_t_and_p']['intergenic']['alt_genomic_loci'] = cp_scaff_lift

            else:
                bring_order['hgvs_t_and_p'] = self.t_and_p_descriptions
        except AttributeError:
            bring_order['hgvs_t_and_p'] = None
        brought_order = {str(self.variant_description): bring_order}
        return brought_order
        
    def collect_metadata(self):
        meta = collections.OrderedDict()
        meta['api_version'] = self.vfo.version
        meta['hgvs_version'] = self.vfo.hgvsVersion
        meta['uta_schema'] = self.vfo.utaVersion
        meta['seqrepo_db'] = self.vfo.seqrepoVersion
        return meta
                

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

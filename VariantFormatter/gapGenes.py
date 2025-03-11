# -*- coding: utf-8 -*-

# Import python modules
import re
import vvhgvs.exceptions
import vvhgvs.assemblymapper
import vvhgvs.variantmapper
import VariantValidator
import VariantValidator.modules.seq_data as seq_data
import VariantValidator.modules.gapped_mapping
from VariantValidator.modules.variant import Variant


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


def compensate_g_to_t(hgvs_tx,
                      hgvs_genomic,
                      vm,
                      hn,
                      reverse_normalizer,
                      primary_assembly,
                      hdp,
                      vfo,
                      transcript_model="refseq",
                      mixed_transcrript_model=True):

    # Set Variable
    re_hash_hgvs_genomic = hgvs_genomic

    # Not required in these instances
    gene_symbol = hdp.get_tx_identity_info(hgvs_tx.ac)[6]

    # Check the blacklist
    gap_compensation = gap_black_list(gene_symbol)

    # Set alt_aln_method
    if mixed_transcrript_model is not True:
        if transcript_model == "ensembl":
            alt_aln_method = 'genebuild'
        else:
            alt_aln_method = 'splign'
    else:
        if "ENST" in hgvs_tx.ac:
            alt_aln_method = 'genebuild'
        else:
            alt_aln_method = 'splign'

    if gap_compensation is False:
        normalized_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn,
                                        reverse_normalizer, vm, vfo)
        hgvs_tx_returns = [normalized_tx, False, None, None, None]
    else:

        """
        Code now modified to use native VV gap mapper
        """
        # VV functions require evm instances
        no_norm_evm = vvhgvs.assemblymapper.AssemblyMapper(hdp,
                                                           assembly_name=primary_assembly,
                                                           alt_aln_method=alt_aln_method,
                                                           normalize=False,
                                                           replace_reference=True
                                                           )
        evm = vvhgvs.assemblymapper.AssemblyMapper(hdp,
                                                   assembly_name=primary_assembly,
                                                   alt_aln_method=alt_aln_method,
                                                   normalize=True,
                                                   replace_reference=True
                                                   )

        # Set validator attributes
        vfo.select_transcripts = 'all'
        vfo.alt_aln_method = alt_aln_method

        # Set variant specific attributes
        variant = Variant(str(hgvs_genomic))
        variant.hgvs_genomic = hgvs_genomic
        variant.reverse_normalizer = reverse_normalizer
        variant.hn = hn
        vfo.reverse_hn = reverse_normalizer
        vfo.hn = hn
        variant.evm = evm
        variant.no_norm_evm = no_norm_evm
        variant.vm = vm
        variant.primary_assembly = primary_assembly
        variant.post_format_conversion = hgvs_genomic

        # Map the gap
        gap_mapper = VariantValidator.modules.gapped_mapping.GapMapper(variant, vfo)
        data, nw_rel_var = gap_mapper.gapped_g_to_c([str(hgvs_tx)], select_transcripts_dict={})
        ori = vfo.tx_exons(tx_ac=hgvs_tx.ac, alt_ac=hgvs_genomic.ac, alt_aln_method=vfo.alt_aln_method)
        re_hash_hgvs_genomic, suppress_c_normalization, hgvs_coding = gap_mapper.g_to_t_compensation(ori,
                                                                                                     nw_rel_var[0],
                                                                                                     "")

        # Populate output list
        gap_compensated_tx_2 = [nw_rel_var[0]]
        if "does not represent a true variant" in data["gapped_alignment_warning"] \
                or "may be an artefact of" in data["gapped_alignment_warning"]:
            gap_compensated_tx_2.append(True)
        else:
            gap_compensated_tx_2.append(False)
        try:
            if data["disparity_deletion_in"][0] == 'transcript':
                corrective_action_taken = 'Automap has deleted ' + str(
                    data["disparity_deletion_in"][1][0]) + ' bp from chromosomal reference sequence ' + str(
                    hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence'
            if data["disparity_deletion_in"][0] == 'chromosome':
                corrective_action_taken = 'Automap has added ' + str(
                    data["disparity_deletion_in"][1][0]) + ' bp to chromosomal reference sequence ' + str(
                    hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence '
            gap_compensated_tx_2.append(corrective_action_taken)
        except Exception:
            gap_compensated_tx_2.append(None)
        gap_compensated_tx_2.append(data["gapped_alignment_warning"].replace("the transcripts listed below",
                                                                             nw_rel_var[0].ac))
        gap_compensated_tx_2.append(data["auto_info"].replace("\n", ""))

        if gap_compensated_tx_2[1] is False:
            refresh_hgvs_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn,
                                              reverse_normalizer, vm, vfo)
            gap_compensated_tx_2[0] = refresh_hgvs_tx
        hgvs_tx_returns = gap_compensated_tx_2

    # Create response dictionary for gap mapping output
    hgvs_tx_dict = {'hgvs_transcript': hgvs_tx_returns[0],
                    'gapped_alignment_warning': hgvs_tx_returns[3],
                    'corrective_action': hgvs_tx_returns[2],
                    'gap_position': hgvs_tx_returns[-1],
                    'transcript_accession': hgvs_tx_returns[0].ac,
                    'hgvs_genomic': re_hash_hgvs_genomic
                    }

    return hgvs_tx_dict


"""
Fully normalizes the hgvs_tx variant from the hgvs_genomic variant
Is only activated if the g_to_t_compensation_code IS NOT USED!
"""


def fully_normalize(hgvs_tx, hgvs_genomic, hn, reverse_normalizer, vm, vfo):
    
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

    try:
        hgvs_tx = vm.g_to_t(hgvs_genomic, hgvs_tx.ac)
    except vvhgvs.exceptions.HGVSError:
        pass
    try:
        hgvs_tx = hn.normalize(hgvs_tx)
    except vvhgvs.exceptions.HGVSError:
        pass
        
    return hgvs_tx

    
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

# -*- coding: utf-8 -*-
"""
supportFunctions.py is an extension module containing additional random functions that support the main function set.

The majority of these functions require hgvs Python package top-level functions or sub-functions contained in uta.py and
seqfetcher.py

"""

# IMPORT REQUIRED PYTHON MODULES
import re
# BioPython modules
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

"""
Function which predicts the protein effect of c. inversions
"""


def pro_inv_info(prot_ref_seq, prot_var_seq):
    info = {
        'variant': 'true',
        'prot_del_seq': '',
        'prot_ins_seq': '',
        'edit_start': 0,
        'edit_end': 0,
        'terminate': 'false',
        'ter_pos': 0,
        'error': 'false'
    }

    # Is there actually any variation?
    if prot_ref_seq == prot_var_seq:
        info['variant'] = 'false'
    else:
        # Deal with terminations
        term = re.compile("\*")
        if term.search(prot_var_seq):
            # Set the termination reporter to true
            info['terminate'] = 'true'
            # The termination position will be equal to the length of the variant sequence because it's a TERMINATOR!!!
            info['ter_pos'] = len(prot_var_seq)
            # cut the ref sequence to == size
            prot_ref_seq = prot_ref_seq[0:info['ter_pos']]
            prot_var_seq = prot_var_seq[0:info['ter_pos']]

            # Whether terminated or not, the sequences should now be the same length
            # Unless the termination codon has been disrupted
            if len(prot_var_seq) < len(prot_ref_seq):
                info['error'] = 'true'
                return info
            else:
                # Set the counter
                aa_counter = 0

                # Make list copies of the sequences to gather the required info
                ref = list(prot_ref_seq)
                var = list(prot_var_seq)

                # Loop through ref list to find the first missmatch position
                for aa in ref:
                    if ref[aa_counter] == var[aa_counter]:
                        aa_counter = aa_counter + 1
                    else:
                        break

                # Enter the start position
                info['edit_start'] = aa_counter + 1
                # Remove those elements form the list
                del ref[0:aa_counter]
                del var[0:aa_counter]

                # the sequences should now be the same length
                # Except if the termination codon was removed
                if len(ref) > len(var):
                    info['error'] = 'true'
                    return info
                else:
                    # Reset the aa_counter but to go backwards
                    aa_counter = 0
                    # reverse the lists
                    ref = ref[::-1]
                    var = var[::-1]
                    # Reverse loop through ref list to find the first missmatch position
                    for aa in ref:
                        if var[aa_counter] == '\*':
                            break
                        if aa == var[aa_counter]:
                            aa_counter = aa_counter + 1
                        else:
                            break
                    # Remove those elements form the list
                    del ref[0:aa_counter]
                    del var[0:aa_counter]
                    # re-reverse the lists
                    ref = ref[::-1]
                    var = var[::-1]

                    # If the var is > ref, the ter has been removed, need to re-add ter to each
                    if len(ref) < len(var):
                        ref.append('*')
                        if prot_var_seq[-1] == '*':
                            var.append('*')
                    # the sequences should now be the same length
                    # Except if the ter was removed
                    if len(ref) > len(var):
                        info['error'] = 'true'
                        return info
                    else:
                        # Enter the sequences
                        info['prot_del_seq'] = ''.join(ref)
                        info['prot_ins_seq'] = ''.join(var)
                        info['edit_end'] = info['edit_start'] + len(ref) - 1
                        return info


"""
Function which predicts the protein effect of c. inversions
"""

def pro_delins_info(prot_ref_seq, prot_var_seq):
    info = {
        'variant': 'true',
        'prot_del_seq': '',
        'prot_ins_seq': '',
        'edit_start': 0,
        'edit_end': 0,
        'terminate': 'false',
        'ter_pos': 0,
        'error': 'false'
    }

    # Is there actually any variation?
    if prot_ref_seq == prot_var_seq:
        info['variant'] = 'false'
    else:
        # Deal with terminations
        term = re.compile("\*")
        if term.search(prot_var_seq):
            # Set the termination reporter to true
            info['terminate'] = 'true'
            # The termination position will be equal to the length of the variant sequence because it's a TERMINATOR!!!

            # Set the te pos dependant on the shortest sequence
            if len(prot_var_seq) <= len(prot_ref_seq):
                info['ter_pos'] = len(prot_ref_seq)
            else:
                info['ter_pos'] = len(prot_var_seq)

            # cut the ref sequence to == size
            prot_ref_seq = prot_ref_seq[0:info['ter_pos']]
            prot_var_seq = prot_var_seq[0:info['ter_pos']]

            # Set the counter
            aa_counter = 0

            # Make list copies of the sequences to gather the required info
            ref = list(prot_ref_seq)
            var = list(prot_var_seq)

            # Loop through ref list to find the first missmatch position
            for aa in ref:
                if ref[aa_counter] == var[aa_counter]:
                    aa_counter = aa_counter + 1
                else:
                    break

            # Enter the start position
            info['edit_start'] = aa_counter + 1
            # Remove those elements form the list
            del ref[0:aa_counter]
            del var[0:aa_counter]

            # Reset the aa_counter but to go backwards
            aa_counter = 0
            # reverse the lists
            ref = ref[::-1]
            var = var[::-1]
            # Reverse loop through ref list to find the first missmatch position
            for aa in ref:
                try:
                    if var[aa_counter] == '\*':
                        break
                except IndexError:
                    break
                if aa == var[aa_counter]:
                    aa_counter = aa_counter + 1
                else:
                    break
            # Remove those elements form the list
            del ref[0:aa_counter]
            del var[0:aa_counter]
            # re-reverse the lists
            ref = ref[::-1]
            var = var[::-1]

            # Enter the sequences
            info['prot_del_seq'] = ''.join(ref)
            info['prot_ins_seq'] = ''.join(var)
            info['edit_end'] = info['edit_start'] + len(ref) - 1
            return info


"""
Translate c. reference sequences, including those that have been modified 
must have the CDS in the specified position
"""


def translate(ed_seq, cds_start):
    # ed_seq = ed_seq.replace('\n', '')
    ed_seq = ed_seq.strip()
    # Ensure the starting codon is in the correct position
    met = ed_seq[cds_start:cds_start + 3]
    if (met == 'ATG') or (met == 'atg'):
        # Remove the 5 prime UTR
        sequence = ed_seq[cds_start:]
        coding_dna = Seq(str(sequence), IUPAC.unambiguous_dna)
        # Translate
        trans = coding_dna.translate()
        aain = list(trans)
        aaout = []
        count = 0
        while aain:
            if aain[count] != '*':
                aaout.append(aain[count])
                count = count + 1
            else:
                aaout.append(aain[count])
                break
        translation = ''.join(aaout)
        # Apply a width of 60 characters to the string output
        # translation = textwrap.fill(translation, width=60)
        return translation
    else:
        translation = 'error'
        return translation


"""
Convert single letter amino acid code to 3 letter code
"""


def one_to_three(seq):
    aacode = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
        'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
        'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Ter'}

    oned = list(seq)
    out = []
    for aa in oned:
        get_value = aacode.get(aa)
        out.append(get_value)

    threed_up = ''.join(out)

    return threed_up


""" 
Takes a reference sequence and inverts the specified position
"""


# n. Inversions - This comes from VariantValidator, not validation!!!!
def n_inversion(ref_seq, del_seq, inv_seq, interval_start, interval_end):
    # Open a list to store the fasta file
    sequence = ''

    # Use string indexing to check whether the sequences are the same
    test = ref_seq[interval_start - 1:interval_end]

    if test == del_seq:
        sequence = ref_seq[0:interval_start - 1] + inv_seq + ref_seq[interval_end:]
        return sequence
    else:
        sequence = 'error'
        return sequence


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
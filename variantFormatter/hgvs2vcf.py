# -*- coding: utf-8 -*-
"""
hgvs2vcf.py
report_hgvs2vcf
Used to report the Most true representation of the VCF i.e. 5 prime normalized but no
additional bases added. NOTE: no gap handling capabilities

hgvs2vcf

Used to report a loose representation of the VCF i.e. 5 prime normalized but no
additional bases added. NOTE: no gap handling capabilities
"""

# Import  modules
import re
import supportedChromosomeBuilds

# Import Biopython modules
from Bio.Seq import Seq



def report_hgvs2vcf(hgvs_genomic_variant, primary_assembly, hdp, reverse_normalize, sf):
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)

    # UCSC Chr
    ucsc_chr = supportedChromosomeBuilds.to_chr_num_ucsc(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if ucsc_chr is not None:
        pass
    else:
        ucsc_chr = reverse_normalized_hgvs_genomic.ac

    # GRC Chr
    grc_chr = supportedChromosomeBuilds.to_chr_num_refseq(reverse_normalized_hgvs_genomic.ac, primary_assembly)
    if grc_chr is not None:
        pass
    else:
        grc_chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        alt_start = start - 1  #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), alt_start, end - 1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq

    # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, start)
        # Assemble
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base


    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        bs = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start - 1, adj_start)
        # Assemble
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())

    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(
            reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base - 1)
        adj_start = start - 1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), start, end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        # pos = str(start)
        # ref = vcf_del_seq
        # alt = vcf_del_seq[:1] + ins_seq
        pos = str(start + 1)
        ref = vcf_del_seq[1:]
        alt = ins_seq

    # Duplications
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)  #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2  #
        start = start - 1  #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac), adj_start, end)
        # Assemble
        pos = str(start+1)
        ref = vcf_ref_seq[1:]
        alt = vcf_ref_seq[1:] + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''

    # Dictionary the VCF
    vcf_dict = {'pos': str(pos), 'ref': ref, 'alt': alt, 'ucsc_chr': ucsc_chr, 'grc_chr': grc_chr,
                'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict


def hgvs2vcf(hgvs_genomic):
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalize.normalize(hgvs_genomic_variant)
    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Chr
    chr = supported_chromosome_builds.to_chr_num(reverse_normalized_hgvs_genomic.ac)
    if chr is not None:
        pass
    else:
        chr = reverse_normalized_hgvs_genomic.ac

    if re.search('[GATC]+\=', str(reverse_normalized_hgvs_genomic.posedit)):
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos.start)
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.ref

    # Insertions    
    elif (re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('del', str(reverse_normalized_hgvs_genomic.posedit))):                       
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base) 
        alt_start = start - 1 #
        # Recover sequences
        ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),alt_start,end-1)
        ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        # Assemble  
        pos = start
        ref = ref_seq
        alt = ref_seq + ins_seq                     

    # Substitutions
    elif re.search('>', str(reverse_normalized_hgvs_genomic.posedit)):
        ref = reverse_normalized_hgvs_genomic.posedit.edit.ref
        alt = reverse_normalized_hgvs_genomic.posedit.edit.alt
        pos = str(reverse_normalized_hgvs_genomic.posedit.pos)

    # Deletions
    elif re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and not re.search('ins', str(reverse_normalized_hgvs_genomic.posedit)):                     
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2
        start = start - 1
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''        
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
        pre_base = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,start)
        # Assemble  
        pos = str(start)
        ref = pre_base + hgvs_del_seq
        alt = pre_base

    
    # inv
    elif re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):                        
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start -1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''        
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
        bs  = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start-1, adj_start)
        # Assemble  
        pos = str(start)
        # pos = str(start-1)
        # ref = bs + vcf_del_seq
        ref = vcf_del_seq
        alt = ins_seq
        if re.search('inv', str(reverse_normalized_hgvs_genomic.posedit)):
            my_seq = Seq(vcf_del_seq)
            # alt = bs + str(my_seq.reverse_complement())
            alt = str(my_seq.reverse_complement())
             
    
    # Delins
    elif (re.search('del', str(reverse_normalized_hgvs_genomic.posedit)) and re.search('ins', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base)
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base -1)
        adj_start = start -1
        start = start
        try:
            ins_seq = reverse_normalized_hgvs_genomic.posedit.edit.alt
        except:
            ins_seq = ''
        else:
            if str(ins_seq) == 'None':
                ins_seq = ''        
        # Recover sequences
        hgvs_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),start,end)
        vcf_del_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
        # Assemble  
        pos = str(start)
        ref = vcf_del_seq
        alt = vcf_del_seq[:1] + ins_seq         
    
    
    # Duplications                              
    elif (re.search('dup', str(reverse_normalized_hgvs_genomic.posedit))):
        end = int(reverse_normalized_hgvs_genomic.posedit.pos.end.base) #
        start = int(reverse_normalized_hgvs_genomic.posedit.pos.start.base)
        adj_start = start - 2 #
        start = start - 1 #
        # Recover sequences
        dup_seq = reverse_normalized_hgvs_genomic.posedit.edit.ref
        vcf_ref_seq = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),adj_start,end)
        # Assemble  
        pos = str(start)
        ref = vcf_ref_seq
        alt = vcf_ref_seq + dup_seq
    else:
        chr = ''
        ref = ''
        alt = ''
        pos = ''
        
    
    # ensure as 5' as possible
    if chr != '' and pos != '' and ref != '' and alt != '':
        if len(ref) > 1:
            rsb = list(str(ref))
            if reverse_normalized_hgvs_genomic.posedit.edit.type == 'identity':
                pos = int(pos) - 1
                prev = sf.fetch_seq(str(reverse_normalized_hgvs_genomic.ac),pos-1,pos)
                pos = str(pos)              
                ref = prev + ref
                alt = prev + alt
    

    # Dictionary the VCF
    vcf_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'normalized_hgvs': reverse_normalized_hgvs_genomic}
    return vcf_dict




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
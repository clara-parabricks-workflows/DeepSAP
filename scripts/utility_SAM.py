# Copyright 2025 NVIDIA CORPORATION
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import subprocess
from datetime import datetime

from utility_GTF import get_region_by_locus

class CIGAR_Eelment:
    def __init__(self, character, ref_bases, read_bases):
        self.character  = character
        self.ref_bases  = ref_bases
        self.read_bases = read_bases

class Read_Info:
    def __init__(self, 
                is_unmapped:bool,
                is_bad_read:bool,
                read_score:int, 
                junctions_score:int, 
                regions:set, 
                junctions_ids:list, 
                n_junctions:int, 
                is_novel:bool,
                is_clipped:bool, 
                new_CIGAR:str, 
                start_pos_ref:int, 
                end_pos_ref:int, 
                n_bases_read:int,  
                nj_tag:str, 
                jt_tag:str, 
                js_tag:str):
        """
            Read information class 
        Args:
            read_score (int): _description_
            junctions_score (int): _description_
            junctions_ids (list): _description_
            n_junctions (int): _description_
            is_clipped (bool): _description_
            new_CIGAR (str): _description_
            start_pos_ref (int): _description_
            end_pos_ref (int): _description_
            n_bases_read (int): _description_
            nj_tag (str): _description_
            jt_tag (str): _description_
            js_tag (str): _description_
        """
        self.is_unmapped      = is_unmapped
        self.is_bad_read      = is_bad_read
        self.read_score       = read_score
        self.junctions_score  = junctions_score
        self.regions          = regions
        self.junctions_ids    = junctions_ids
        self.n_junctions      = n_junctions
        self.is_novel         = is_novel
        
        self.is_clipped       = is_clipped
        self.new_CIGAR        = new_CIGAR
        self.start_pos_ref    = start_pos_ref
        self.end_pos_ref      = end_pos_ref
        self.n_bases_read     = n_bases_read
        
        self.nj_tag = nj_tag
        self.jt_tag = jt_tag
        self.js_tag = js_tag
    
    def __str__(self):
        return f"Read Info: is_unmapped={self.is_unmapped}\nis_bad_read={self.is_bad_read}\nread_score={self.read_score}\njunctions_score={self.junctions_score}\nis_clipped={self.is_clipped}\n)"

def get_CIGAR(CIGAR_elements:list):
    CIGAR = ""
    
    for element in CIGAR_elements:
        character = element.character
        
        if character in ['M', 'X', '=', 'I', 'S', 'H']:
            CIGAR += str(element.read_bases)
        elif character in ['D', 'N']:
            CIGAR += str(element.ref_bases)
        else:
            assert True == False

        CIGAR += character
    
    return CIGAR


def get_strand(FLAG:int, sam_record_fields:list):
    if "XS:A:-" in sam_record_fields:
        return '-'
    elif "XS:A:+" in sam_record_fields:
        return '+'
    else:
        if ('N' in sam_record_fields[5]):
            print("[{}]\t[Warning]\tProblematic record due to spliced-read without XS tag: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), sam_record_fields[0]))
            
        if is_reverse_strand(FLAG):
            return '-'
        else:
            return '+'

    # Old method
    # if is_reverse_strand(FLAG):
    #     if "XS:A:-" in sam_record_fields:
    #         return '-'
    #     else:
    #         return '+'
    # else:
    #     if "XS:A:-" in sam_record_fields:
    #         return '-'
    #     else:
    #         return '+'
                                        
def is_mate_1(FLAG:int):
    return FLAG & 0x40

def is_mate_2(FLAG:int):
    return FLAG & 0x80

def is_secondary_alignment(FLAG:int):
    return (int(FLAG) & 0x100) != 0 # Check if the 0x100 bit (512th bit) is set

def is_proper_pair(FLAG:int):
    return FLAG & 0x1 and FLAG & 0x2 # Read is paired and mapped in proper pair

def is_reverse_strand(FLAG:int):
    return (int(FLAG) & 0x10) != 0  # Check if the 5th bit is set

def are_mates_from_same_transcript(found_junctions:dict, mate1_junctions:list, mate2_junctions:list):
    return True

def get_first_read(path_sam:str, samtools:str):
    read_id = ""
    
    if path_sam[-4:] == ".bam":
        records = subprocess.Popen([samtools, "view", "-h", path_sam], stdout=subprocess.PIPE)
    else:
        records = subprocess.Popen(["cat", path_sam], stdout=subprocess.PIPE)
    
    while True:  
        line = records.stdout.readline().decode()
        
        if len(line) == 0:
            break
        elif(line.startswith('@')):
            continue
        else:
            read_id  = line.split('\t')[0]
            break
            
    return read_id

def is_bad_novel_deletion(deletion_size:int, n_overhang_bases:int, match_reward:int, gap_open_penalty:int, gap_ext_penalty:int):
    penalty = gap_open_penalty + (gap_ext_penalty * deletion_size)
    reward  = match_reward * n_overhang_bases

    if penalty >= reward:
        return True 
    else:
        return False

def adjust_min_splice_junction_penalty(min_novel_score, n_overhang_bases):
    rate = 0.33
    min_penalty = max(0, min_novel_score - 20)
    adjusted_penalty = min_novel_score - rate * n_overhang_bases
    return max(min_penalty, adjusted_penalty)

def adjust_min_qual(min_averge_qual, n_overhang_bases):
    rate = 0.1
    min_qual = max(0, min_averge_qual - 5)
    adjusted_qual = min_averge_qual - rate * n_overhang_bases
    return max(min_qual, adjusted_qual)

def is_bad_novel_junction(junction_score, min_novel_score, 
                          n_overhang_bases, min_overhang_bases_novel, 
                          bases_qualities, min_averge_qual, ASCII_base):
    """
    ex1:
        junction_score           = 20
        novel_junction_penalty   = 30
        n_overhang_bases         = 10
        min_overhang_bases_novel = 8
        
        n_overhang_bases < min_overhang_bases_novel = False
        
        score = (100 - 20 ) / 100      = 0.8
        score * novel_junction_penalty = 24
        n_overhang_bases <  score * novel_junction_penalty = 
        10 < 24 = False --> Soft-clip The Junction
        
    ex2:
        junction_score           = 90
        novel_junction_penalty   = 30
        n_overhang_bases         = 10
        min_overhang_bases_novel = 8
        
        n_overhang_bases < min_overhang_bases_novel = False
        
        score = (100 - 90 ) / 100      = 0.1
        score * novel_junction_penalty = 3
        n_overhang_bases <  score * novel_junction_penalty = 
        10 < 3 = False --> Keep The Junction
    """
        
    assert(n_overhang_bases == len(bases_qualities))
    
    is_bad_junctions = False
    
    if n_overhang_bases < min_overhang_bases_novel:
        is_bad_junctions = True

    elif junction_score < min_novel_score:
        is_bad_junctions = True

    else:
        averge_qual = sum([ (ord(Q) - ASCII_base) for Q in bases_qualities])/len(bases_qualities)
        
        if averge_qual < min_averge_qual :
            is_bad_junctions = True
    
    return is_bad_junctions

def is_bad_annotated_junction(junction_score, min_annotated_score,
                              n_overhang_bases,
                              min_overhang_bases_annotated_solo,  
                              min_overhang_bases_annotated_multi,
                              bases_qualities, min_averge_qual, ASCII_base,
                              n_hits):
    
    assert(n_overhang_bases == len(bases_qualities))
    
    is_bad_junctions   = False
    min_overhang_bases = min_overhang_bases_annotated_multi
    
    if n_hits <= 2:
        min_overhang_bases = min_overhang_bases_annotated_solo
        
    if n_overhang_bases < min_overhang_bases:
        is_bad_junctions = True
    
    elif junction_score < min_annotated_score:
        is_bad_junctions = True

    else:
        averge_qual = sum([ (ord(Q) - ASCII_base) for Q in bases_qualities])/len(bases_qualities)
        
        if averge_qual < min_averge_qual :
            is_bad_junctions = True
    
    return is_bad_junctions

def get_n_soft_clipped_bases(CIGAR):
    n_bases = ""
    
    for character in CIGAR:
        if character.isdigit():
            n_bases += character
        
        elif character == 'S':
            return int(n_bases)
        else:
            n_bases = ""
            continue
    return 0

def get_read_score_junctions_regions(found_junctions:dict, 
                       junctions_scores:dict,
                       regions_array_table:dict,
                       min_overhang_bases_novel:int,
                       min_overhang_bases_annotated_solo:int,  
                       min_overhang_bases_annotated_multi:int,
                       min_novel_intron_size:int,
                       annotated_junction_reward:int,
                       novel_junction_penalty:int,
                       min_novel_score:int,
                       min_annotated_score:int,
                       gap_open_penalty:int,
                       gap_ext_penalty:int,
                       mismatch_penalty:int,
                       match_reward:int,
                       soft_clipping,
                       min_averge_qual, ASCII_base,
                       n_hits:int,
                       sam_record_fields:list):
    
    def add_regions(locus_ref:int, n_bases:int, regions:set):
        first_region = get_region_by_locus(chrom, locus_ref,               regions_array_table) # all regions are 0-based
        last_region  = get_region_by_locus(chrom, locus_ref + n_bases - 1, regions_array_table)
        
        # print("first_region: {}", first_region)
        # print("last_region:  {}", last_region)
        
        if first_region != None and last_region != None:
            if first_region == -1 or last_region == -1:
                regions.add(first_region)
                regions.add(last_region)
            else:
                for region in range(first_region, last_region + 1):
                    regions.add(region)
            
    FLAG    = int(sam_record_fields[1])
    chrom   = sam_record_fields[2]
    start   = int(sam_record_fields[3])
    MAPQ    = int(sam_record_fields[4])
    CIGAR   = str(sam_record_fields[5])
    strand  = get_strand(FLAG, sam_record_fields)
    n_bases = ""
    
    n_read_bases  = len(sam_record_fields[9])
    is_clipped    = False
    start_pos_ref = start - 1
    end_pos_ref   = start - 1
    n_bases_read  = 0
    
    locus_ref      = start - 1 # 0 based
    locus_read     = 0
    
    regions        = set()
    junctions_ids  = []
    CIGAR_elements = []
    read_score     = 0
    n_edits        = 0
    n_bases_matched = 0  # bases aligned using M or = 

    is_right_anchor_clipped = False         # To allow calculating the read score even if soft-clipping was done
                
    for character_index, character in enumerate(CIGAR):
        if character.isdigit():
            n_bases += character
            
        elif character == 'M':
            
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('M', n_bases, n_bases))
            
            add_regions(locus_ref, n_bases, regions)
            
            n_bases_matched += n_bases
            
            locus_ref    += n_bases
            locus_read   += n_bases
            n_bases_read += n_bases
            read_score   += (n_bases * match_reward)
            n_bases       = ""
        
        elif character == '=':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('=', n_bases, n_bases))
            
            add_regions(locus_ref, n_bases, regions)
            
            n_bases_matched += n_bases
            
            locus_ref    += n_bases
            n_bases_read += n_bases
            locus_read   += n_bases
            read_score   += (n_bases * match_reward)
            n_bases       = ""
            
            is_mordern_sam = True
            
        elif character == 'X':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('X', n_bases, n_bases))
            
            add_regions(locus_ref, n_bases, regions)
            
            locus_ref    += n_bases
            n_bases_read += n_bases
            locus_read   += n_bases
            read_score   -= (n_bases * mismatch_penalty)
            n_bases       = ""
            
        elif character == 'D':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('D', n_bases, 0))
            
            # if n_bases > min_novel_intron_size:
                # print("[{}]\t[Warning]\tProblematic record due to deletion size is greater than min_novel_intron_size: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), sam_record_fields[0]))
                
            locus_ref   += n_bases
            read_score  -= (gap_open_penalty + (n_bases * gap_ext_penalty))
            n_bases      = ""
            n_edits     += 1
            
        elif character == 'I':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('I', 0, n_bases))
            
            n_bases_read += n_bases
            locus_read   += n_bases
            read_score   -= (gap_open_penalty + (n_bases * gap_ext_penalty))
            n_bases       = ""
            n_edits      += 1
            
        elif character == 'S':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('S', 0, n_bases))
            
            locus_read += n_bases
            n_bases     = ""
            
        elif character == 'H':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('H', 0, n_bases))
            n_bases     = ""
            
        elif character == 'N':
            n_bases     = int(n_bases)
            CIGAR_elements.append(CIGAR_Eelment('N', n_bases, 0))
            intron_size   = n_bases
            donor_site    = locus_ref
            locus_ref    += n_bases
            
            acceptor_site = locus_ref
            n_bases       = ""
            n_edits      += 1
            
            splice_id     = chrom + "__" + strand + "__" + str(donor_site) + "__" + str(acceptor_site)
            junction_type = found_junctions[splice_id].transcript_type
            assert splice_id in found_junctions.keys()
            
            # if junction_type == 'Novel' and intron_size < min_novel_intron_size:
            #     CIGAR_elements[-1].character = 'D'
            #     read_score  -= (gap_open_penalty + (intron_size * gap_ext_penalty))
            #     is_clipped   = True
            #     continue
            
            junctions_ids.append(splice_id)
            
            if soft_clipping and is_right_anchor_clipped == False:
                junction_score  = junctions_scores[splice_id]
                is_bad_left_anchor  = False
                is_bad_right_anchor = False
                
                n_soft_clipped = get_n_soft_clipped_bases(CIGAR[character_index + 1:]) # get the number of softclipped bases to avoid short-anchor and soft-clipped bases 85M14755N6M9S
                
                bases_qualities = ""
                
                if len(sam_record_fields[9]) ==  len(sam_record_fields[10]):
                    bases_qualities = sam_record_fields[10]
                else:
                    # If SAM file doesnt have qualities 
                    good_quality = ASCII_base + 25
                    bases_qualities = len(sam_record_fields[9]) * chr(good_quality)

                left_flank  = bases_qualities[:n_bases_matched]
                right_flank = bases_qualities[locus_read:n_read_bases - n_soft_clipped]

                if len(left_flank) == 0 or len(right_flank) == 0:
                        print("[{}]\t[Error]\tBad overhanged bases of a junction, Left: {}, Right: {} Read: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),left_flank , right_flank, sam_record_fields[0]))

                elif junction_type == 'Novel' and intron_size < min_novel_intron_size:
                    del junctions_ids[-1]
                    del CIGAR_elements[-1]
                    CIGAR_elements.append(CIGAR_Eelment('D', intron_size, 0))
                    is_clipped     = True

                    is_bad_left_anchor  = is_bad_novel_deletion(intron_size, len(left_flank),  match_reward, gap_open_penalty, gap_ext_penalty)
                    is_bad_right_anchor = is_bad_novel_deletion(intron_size, len(right_flank), match_reward, gap_open_penalty, gap_ext_penalty)

                elif junction_type == 'Novel':
                    is_bad_left_anchor = is_bad_novel_junction(junction_score, min_novel_score, 
                                                               n_bases_matched, 
                                                               min_overhang_bases_novel,
                                                               left_flank, min_averge_qual, ASCII_base)
                    
                    is_bad_right_anchor = is_bad_novel_junction(junction_score, min_novel_score,  
                                                                n_read_bases - locus_read - n_soft_clipped,
                                                                min_overhang_bases_novel,
                                                                right_flank, min_averge_qual, ASCII_base)
                else:
                    is_bad_left_anchor = is_bad_annotated_junction(junction_score, min_annotated_score,
                                                                   n_bases_matched, 
                                                                   min_overhang_bases_annotated_solo,  
                                                                   min_overhang_bases_annotated_multi,
                                                                   left_flank, min_averge_qual, ASCII_base,
                                                                   n_hits)
                    
                    is_bad_right_anchor = is_bad_annotated_junction(junction_score, min_annotated_score,
                                                                    n_read_bases - locus_read - n_soft_clipped, 
                                                                    min_overhang_bases_annotated_solo,  
                                                                    min_overhang_bases_annotated_multi,
                                                                    right_flank, min_averge_qual, ASCII_base,
                                                                    n_hits)
                if is_bad_left_anchor:
                    is_clipped     = True
                    start_pos_ref  = locus_ref
                    n_bases_read   = 0
                    CIGAR_elements = [CIGAR_Eelment('S', 0, locus_read)]
                    junctions_ids  = list()
                    regions        = set()
                    
                elif is_bad_right_anchor:
                    is_clipped     = True
                    is_right_anchor_clipped = True
                    
                    junctions_ids  = junctions_ids[:-1]
                    CIGAR_elements = CIGAR_elements[:-1]
                    CIGAR_elements.append(CIGAR_Eelment('S', 0, n_read_bases - locus_read))
                    
                    tmp_junctions_ids  = list(junctions_ids)      # To allow calculating the read score even if soft-clipping was done
                    tmp_CIGAR_elements = list(CIGAR_elements)
                    tmp_regions        = set(regions)
                    
                    end_pos_ref = locus_ref - intron_size
                    
    
    if is_right_anchor_clipped == True:
        junctions_ids  = tmp_junctions_ids
        CIGAR_elements = tmp_CIGAR_elements
        regions        = tmp_regions
        
    else:
        end_pos_ref = locus_ref 
        
    if len(junctions_ids) != 0:
        
        nj_tag = "nj:i:{}".format(len(junctions_ids))
        jt_tag = "jt:Z:{}".format(",".join([found_junctions[junction_id].transcript_type for junction_id in junctions_ids]))
        js_tag = "js:Z:{}".format(",".join([str(junctions_scores[junction_id]) for junction_id in junctions_ids]))
        
        junctions_score = 0
        is_novel        = False
        
        for junction_id in junctions_ids:
            score = 0
            
            if found_junctions[junction_id].transcript_type == "Novel": 
                score += junctions_scores[junction_id]
                score -= novel_junction_penalty
                
                is_novel = True
            else:
                score += junctions_scores[junction_id] 
                score += annotated_junction_reward     # I am not sure! no matter how good the model is, there are annotated junctions that will have low score
            
            if score < 0:
                score = 0
            elif score > 99:
                score = 99
                
            junctions_score += score
            
        junctions_score = int(junctions_score / len(junctions_ids))
            
        read_score += junctions_score
        
        return Read_Info(False, False, read_score, junctions_score, regions, junctions_ids, len(junctions_ids), is_novel,
                         is_clipped, get_CIGAR(CIGAR_elements), start_pos_ref + 1, end_pos_ref, n_bases_read,
                         nj_tag, jt_tag, js_tag)
    else:
        junctions_score = 99
        read_score     += junctions_score
        
        return Read_Info(False, False, read_score, junctions_score, regions, junctions_ids, 0, False,
                         is_clipped, get_CIGAR(CIGAR_elements), start_pos_ref + 1, end_pos_ref, n_bases_read,
                         "", "", "")


def get_MAPQ(mate_hits:set):
    if len(mate_hits) == 1:
        return 255
    elif len(mate_hits) == 2:
        return 3
    elif len(mate_hits) == 3:
        return 2
    elif len(mate_hits) == 4:
        return 1
    else:
        return 0

def is_primary_alignment(flag):
    return (flag & 0x100) == 0  # Check if the 0x100 bit is set

def set_primary_alignment_flag(flag):
    if is_primary_alignment(flag):
        return flag
    else:
        return flag & ~0x100 # Clear the 0x100 bit to mark as primary alignment
    
def set_secondary_alignment_flag(flag):
    if is_primary_alignment(flag):
        return flag | 0x100  # Set the 0x100 bit to mark as secondary alignment
    else:
        return flag

def is_pair(mate1:list, mate2:list):
    if mate1[2] == mate2[2] and mate1[3] == mate2[7] and mate1[7] == mate2[3] :
        return True
    else:
        return False
    
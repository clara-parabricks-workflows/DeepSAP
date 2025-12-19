
#!/usr/bin/env python3

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


import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import math
import subprocess
import argparse
import os
import sys
import json
from statistics import median, mean, mode
from sys import exit
import pandas as pd
from math import log10
import random
from datetime import datetime
import time
import concurrent.futures
import shutil
import gc
import re 
from cryptography.fernet import Fernet
import tarfile
import re

from predict import predict
from utility_SAM import Read_Info
from utility_SAM import get_strand, get_read_score_junctions_regions, get_first_read, is_proper_pair, is_pair, is_mate_1, is_mate_2, get_MAPQ, set_primary_alignment_flag, set_secondary_alignment_flag
from utility_FASTA import parse_fasta, predict_splice_junctions_from_FASTA, parse_motif_pairs, parse_motifs_arg
from utility_GTF import parse_gtf, get_transcripts_info, Junction, generate_regions, are_mates_from_same_transcript,JunctionFromRead, collect_junctions_sequences_from_GTF

os.environ["HUGGINGFACE_HOME"] = "/root/.cache/huggingface"


def is_bad_read(read_hit, min_averge_qual_per_read, ASCII_base):
    is_bad = False

    if len(read_hit[9]) == len(read_hit[10]):
        bases_qualities = read_hit[10]
        averge_qual = sum([ (ord(Q) - ASCII_base) for Q in bases_qualities])/len(bases_qualities)

        if averge_qual < min_averge_qual_per_read:
            read_hit[4] = 0
            is_bad = True
    return is_bad

def reverse_complement(sequence):
    reversed_bases = {'N':'N', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    reversed_sequence = "".join([reversed_bases[base] for base in sequence[::-1]])
        
    return reversed_sequence

def filter_and_add_junctions_info_to_SAM(found_junctions, 
                                         junctions_scores,
                                         config_file,
                                         path_sam_in:str, 
                                         path_out:str,
                                         output_prefix:str,
                                         samtools:str,
                                         regions_loci_table , regions_array_table, regions_relations_table, regions_transcripts_table, jumps_transcripts_table,
                                         version, command,
                                         is_verbose=False,
                                         wanted_read=""):
    
    pre_read_id      = get_first_read(path_sam_in, samtools)
    read_hits = []
    
    with open(config_file, "r") as f:
        parameters_config = json.load(f)
    
    soft_clipping = parameters_config['soft_clipping']
    min_overhang_bases_novel           = parameters_config['min_overhang_bases_novel']
    min_overhang_bases_annotated_solo  = parameters_config['min_overhang_bases_annotated_solo']
    min_overhang_bases_annotated_multi = parameters_config['min_overhang_bases_annotated_multi']
    min_novel_intron_size    = parameters_config['min_novel_intron_size']
    novel_transcript_penalty  = parameters_config['novel_transcript_penalty']
    novel_junction_penalty    = parameters_config['novel_junction_penalty']
    min_novel_score           = parameters_config['min_novel_score']
    min_annotated_score       = parameters_config['min_annotated_score']
    
    annotated_junction_reward = parameters_config['annotated_junction_reward']
    max_score_margin          = parameters_config['max_score_margin']
    score_margin              = parameters_config['score_margin']
    
    min_averge_qual_per_read = parameters_config['min_averge_qual_per_read']
    min_averge_qual = parameters_config['min_averge_qual']
    ASCII_base      = parameters_config['ASCII_base']
    
    gap_open_penalty = parameters_config['gap_open_penalty']
    gap_ext_penalty  = parameters_config['gap_ext_penalty']
    mismatch_penalty = parameters_config['mismatch_penalty']
    match_reward     = parameters_config['match_reward']
    
    # Reading input BAM 
    if path_sam_in[-4:] == ".bam":
        records = subprocess.Popen([samtools, "view", "-h", path_sam_in], stdout=subprocess.PIPE)
    else:
        records = subprocess.Popen(["cat", path_sam_in], stdout=subprocess.PIPE)
    
    # Writing output BAM
    output_bam = os.path.join(path_out , output_prefix + ".bam")

    samtools_process   = subprocess.Popen([samtools, "view", "-bS", "-o", output_bam , "-"], stdin=subprocess.PIPE)
    command_entry      = "@PG\tID:DeepSAP\tPN:DeepSAP\tVN:{}\tCL:{}\n".format(version, command)
    is_command_written = False
    
    n_sam_records     = 0 
    n_reads_ids       = 0
    n_processed_reads_ids = 0
    is_end_of_file = False

    while True:  
        line = records.stdout.readline().decode()
        is_spliced = False
        
        if(line.startswith('@')):
            if line.startswith("@PG") and is_command_written == False:
                samtools_process.stdin.write(command_entry.encode())
                is_command_written = True
                
            samtools_process.stdin.write(line.encode())
            continue
        else:
            if len(line) == 0 :
                is_end_of_file = True
                read_id = ""
            else:
                n_sam_records += 1
                
                fields  = line[:-1].split('\t')      # remove the new line \n
                read_id = str(fields[0])
        
            if read_id == pre_read_id:
                read_hits.append(fields)
            else:
                n_reads_ids     += 1
                n_hits           = len(read_hits)
                is_spliced       = False 
                are_proper_pairs = True
                
                if is_verbose:
                    if wanted_read in pre_read_id:
                        print(pre_read_id)
                        pass
                    else:
                        pre_read_id = read_id
                        read_hits   = [fields]
                        continue
                
                # Check if a hit is spliced, and regardless get the hit score & junctions
                for i in range(n_hits):
                    if 'N' in read_hits[i][5]:
                        is_spliced = True  
                        
                    if not is_proper_pair(int(read_hits[i][1])): 
                        are_proper_pairs = False
                
                if is_spliced and (are_proper_pairs or n_hits <=2):
                    n_processed_reads_ids += 1
                    
                    read_infos = []
                    
                    for i in range(n_hits):
                        if read_hits[i][2] == '*':
                            read_info = Read_Info(True, False, 0, 0, set(), [], 0, False, False, "", 0, 0, 0, '', '', '')

                        else:
                            read_info = get_read_score_junctions_regions(found_junctions, junctions_scores, regions_array_table,
                                                                        min_overhang_bases_novel, 
                                                                        min_overhang_bases_annotated_solo,
                                                                        min_overhang_bases_annotated_multi,
                                                                        min_novel_intron_size,
                                                                        annotated_junction_reward, novel_junction_penalty, 
                                                                        min_novel_score,
                                                                        min_annotated_score,
                                                                    gap_open_penalty, gap_ext_penalty, mismatch_penalty, match_reward,
                                                                    soft_clipping,
                                                                    min_averge_qual, ASCII_base,
                                                                    n_hits, read_hits[i])
                            
                            if is_bad_read(read_hits[i], min_averge_qual_per_read, ASCII_base):
                                read_info.is_bad_read = True
                                read_hits[i][4] = '0'

                        if is_verbose:
                            print("\nRead:      {}".format(read_hits[i]))
                            print("Junctions:")
                            for jct in read_info.junctions_ids:
                                print("\t{} : {}, {}".format(jct, found_junctions[jct].transcript_type , junctions_scores[jct]))
                        
                        read_infos.append(read_info)
                        
                    if n_hits == 1:
                        if read_infos[0].is_clipped:
                            read_hits[0][3] = str(read_infos[0].start_pos_ref)
                            read_hits[0][5] = read_infos[0].new_CIGAR
                    
                    elif n_hits == 2:
                        if read_infos[0].is_clipped or read_infos[1].is_clipped: 
                            
                            """
                            1- Forward
                            472:1,12061,+,100M;1,12171,-,56M418N44M:AAAAAAAC        99      chr1    12062   40      100M    =       12175   211     GGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC    I=DDIIEFGEEFDD8GGGGGGGGH
                            HHHIIIIIIIIIIIIIIIIIIIIIIIIIGGBGEEIIIIHIHHFDFEHIE<AGIGFFIHFHHFIHIIGHFIIHIIGC    XM:Z:53M418N44M XD:Z:53C43      XN:i:1  MD:Z:100        NH:i:1  HI:i:1  NM:i:0  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XX:Z:GENE.4076:208..307:S:.1.,GENE.289622:194..293:S:.1.
                            ,GENE.505049:189..288:S:.1.,GENE.505050:189..288:S:.1.  XS:A:+
                            472:1,12061,+,100M;1,12171,-,56M418N44M:AAAAAAAC        147     chr1    12175   40      53M418N44M      =       12062   -211    AAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCATATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATC       =BBDBEGHIIIGDDHI
                            HIIIIHIIHIIIIFGBDGGG@(GGDEG8>GFGGGG@@;;0GIHHHHFIGHIIIIIEHHIIHDGGGIIHIIIIIIGGIIII#       XM:Z:100M       XD:Z:100        XN:i:0  MD:Z:53C43      NH:i:1  HI:i:1  NM:i:1  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XX:Z:GENE.505049:302..398:S:.1s|s2.     XY:Z:GEN
                            E.4076:321..450:S:.1s|x2.,GENE.289622:307..436:S:.1s|x2.,GENE.505050:302..449:S:.1s|x2. XS:A:+

                            12062      12162
                                                12175    12272
                            j_start + j_n_bases - i_start + 1 
                            12272 - 12062 + 1= 211

                            2- Forward with Insertion
                            3352:1,790069,+,100M;1,790242,-,41M5I54M        99      chr1    790070  40      100M    =       790245  274     TGATGGCACACTCTCCTTTCTCAGCATTCTTGAAAATGGCGTGTGAGCCAAGCCATCTTACTGAGAAATTACCTTAAACCATTTCATGAACCAAAAATAT    GGGDBBDGHHHIHHHGEBGGIEBBED>BGEHIIIHHFA;45:(DIHIHIFEGGIIIII@@ACDDEIIFHIGIGGCGIHD@E:88@DFHFCGBCEACDA@@    XM:Z:39M5I54M   XD:Z:93 XN:i:5  MD:Z:84G15      NH:i:1  HI:i:1  NM:i:1  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XX:Z:GENE.724:1873..1972:B:.5.,GENE.4212:318..417:B:.1.
                            3352:1,790069,+,100M;1,790242,-,41M5I54M        147     chr1    790245  40      39M5I54M        =       790070  -274    GAAAAATCACAGCTTGGGTGCATGTGAAGCCAGAGGAGCCTCAAACTCCTGGGGCCTGCGGTGCCATCACTCAGCTCCCCTGGGCTTCACATGGCCAT      BHEIE8?AAA@GECCG@IIIHHHFHHHHHGGDB@GDIIIIEEFIHIIIIIFIIHHIIHG@GIIHFIIDEIIIII>BGGGDGGGGBB@GIIHIIIIIIH      XM:Z:100M       XD:Z:84G15      XN:i:1  MD:Z:93 NH:i:1  HI:i:1  NM:i:5  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XX:Z:GENE.724:2048..2140:B:.5.,GENE.4212:493..585:B:.1.

                            790070         790170
                                                    790245     790338+5 
                            j_start + j_n_bases - i_start
                            790343 - 790070 + 1 = 274

                            3- Reverse
                            206:1,13392,+,100M;1,13508,-,100M       83      chr15   102517678       40      100M    =       102517562       -216    GTGCAGCTGCCTGTCAGGAAGAGGCCTACTTCTGGTGAGACTGGGCCGACAAAAGGCAGTGAGAAATGTGATCTCGGGGTGGTGGAGGCTCTAGGGAAAG    CIFDDIIIIGFDDEIEIIICDIIIIICIIIIIIIIICCDDDDVEIIHCCC?HCEICVHHFCCIIVVCCC<EVVHCCIIIIIIFFICIIIIIDDDDDIIIH    XM:Z:100M       XD:Z:57C42      XN:i:1  MD:Z:0C37A7A53  NH:i:1  HI:i:1  NM:i:3  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XY:Z:GENE.345549:1877..1877:U:i16i,GENE.345551:2062..2062:U:i16i
                            206:1,13392,+,100M;1,13508,-,100M       163     chr15   102517562       40      100M    =       102517678       216     NTGGACTCCACACTCTCCTGGGTTTCACCTTTGTAGCAGGATCCCTGCAGACCAGGCTCATGACAAACACCGTCTCCAGCGGGCAGAGCAAAGGAAGGGC    #CVVDDDDDIIIIIIIIIIDCIIDDIIIIIIIICH>HECCDIIGGECCEEIIIDDIIIFEVEV'<5>9VV?50>T<<VGC??GFHHIIIECHHVEEHVHI    XM:Z:100M       XD:Z:0C37A7A53  XN:i:3  MD:Z:57C42      NH:i:1  HI:i:1  NM:i:1  SM:i:40 XQ:i:40 X2:i:0  XO:Z:CU XY:Z:GENE.345549:1877..1877:U:i16i,GENE.345551:2062..2062:U:i16i

                                                102517678  102517778 
                            102517562 102517662
                            i start + n_bases - j strt
                            102517778 - 102517562
                            
                            -4 With Intron
                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      99      chr1    145304633       2       20M1292N52M3903N27M     =       145309982       255     TGCTGTACACATTATACCAGAAAATGAAAGTGATGATGAGGAAGAGGAAGAAAAAGGGCCAGTGTCTCCCAGGAATCTGCAGGAGTCTGAAGAGGAGGA    GIIHHIIGDHHHIIIIIIHDEB?BIIIHIIIIFFHIIIHHHEGGHH>GGFIIHIEEBGIHGDIIIGGIIIHIDDGEEEGIIIFHFHEDD<@CEBIHIHH     XM:Z:82M1039N18M        XD:Z:100        XN:i:0  MD:Z:15T83      NH:i:3HI:i:1   NM:i:1  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.2615:2731..2829:S:.15s|s16s|s17.,GENE.288677:943..1041:S:.6s|s7s|s8.,GENE.468819:1008..1106:S:.7s|s8s|s9.     XS:A:+
                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      147     chr1    145309982       2       82M1039N18M     =       145304633       -255    ATGTTGGCCTCGTACCAGTCTTACAGCAGCACATTTCACTCATTAGAGGAACAGCAAGTCTGCATGGCTGTTGACATAGGCAGACATCGGTGGGATCAAG   <?:==DB?<6>BCG:=9<BGDEIIEEEDBHHGHEEBEHGEIGGIIIIIIIIIIIHHHHHGA.;4EIGBEGHIIIHHG<CCCC@=@?:==2/4=7;38DII    XM:Z:20M1292N52M3903N27M        XD:Z:15T83      XN:i:1  MD:Z:100        NH:i:3HI:i:1   NM:i:0  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.2615:2885..2984:S:.17s|s18.,GENE.288677:1097..1196:S:.8s|s9.,GENE.468819:1162..1261:S:.9s|s10.        XS:A:+

                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      355     chr1    144151705       2       20M1288N52M3907N27M     =       144157054       255     TGCTGTACACATTATACCAGAAAATGAAAGTGATGATGAGGAAGAGGAAGAAAAAGGGCCAGTGTCTCCCAGGAATCTGCAGGAGTCTGAAGAGGAGGA    GIIHHIIGDHHHIIIIIIHDEB?BIIIHIIIIFFHIIIHHHEGGHH>GGFIIHIEEBGIHGDIIIGGIIIHIDDGEEEGIIIFHFHEDD<@CEBIHIHH     XM:Z:82M1042N18M        XD:Z:100        XN:i:0  MD:Z:15T83      NH:i:3HI:i:2   NM:i:1  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.288633:186..284:S:.1s|s2s|s3. XS:A:+
                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      403     chr1    144157054       2       82M1042N18M     =       144151705       -255    ATGTTGGCCTCGTACCAGTCTTACAGCAGCACATTTCACTCATTAGAGGAACAGCAAGTCTGCATGGCTGTTGACATAGGCAGACATCGGTGGGATCAAG   <?:==DB?<6>BCG:=9<BGDEIIEEEDBHHGHEEBEHGEIGGIIIIIIIIIIIHHHHHGA.;4EIGBEGHIIIHHG<CCCC@=@?:==2/4=7;38DII    XM:Z:20M1288N52M3907N27M        XD:Z:15T83      XN:i:1  MD:Z:100        NH:i:3HI:i:2   NM:i:0  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.288633:340..439:S:.3s|s4.     XS:A:+

                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      339     chr1    148330261       2       27M3907N52M1288N20M     =       148329064       -254    TCCTCCTCTTCAGACTCCTGCAGATTCCTGGGAGACACTGGCCCTTTTTCTTCCTCTTCCTCATCATCACTTTCATTTTCTGGTATAATGTGTACAGCA    DDIDIVEG?>HHEDFDFIIICEEECHHIDIIICCIIIHCDICVEEIDIIFCC<DDCCEDDDIIIDFFIIIIDIIIV?VEHDIIIIIIDDDHCIIDDIIC     XM:Z:18M1042N82M        XD:Z:100        XN:i:0  MD:Z:83A15      NH:i:3HI:i:3   NM:i:1  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.16365:588..686:S:.5s|s6s|s7.  XS:A:-
                            398198:1,145304632,+,20M1292N52M3903N28M;1,145309981,-,82M1039N18M      419     chr1    148329064       2       18M1042N82M     =       148330261       254     CTTGATCCCACCGATGTCTGCCTATGTCAACAGCCATGCAGACTTGCTGTTCCTCTAATGAGTGAAATGTGCTGCTGTAAGACTGGTACGAGGCCAACAT   IIH83;7=4/2==:??=?GGGG>CDDIIIDCEVCIE4;.TCDDDDDIIIIIIIIIIICCIECDEVEEDCDDVHEEEIIEHCV>9=:CGV<6>?VH==:?>    XM:Z:27M3907N52M1288N20M        XD:Z:83A15      XN:i:1  MD:Z:100        NH:i:3HI:i:3   NM:i:0  SM:i:2  XQ:i:40 X2:i:40 XO:Z:CM XX:Z:GENE.16365:742..841:S:.7s|s8.      XS:A:-

                            145,304,633       145,304,732
                                                            145,309,982   145,310,064
                            """
                            
                            if read_infos[0].is_unmapped == False and read_infos[0].is_bad_read == False:
                                read_hits[0][3] = str(read_infos[0].start_pos_ref)
                                read_hits[0][5] = read_infos[0].new_CIGAR

                            if read_infos[1].is_unmapped == False and read_infos[1].is_bad_read == False:
                                read_hits[1][3] = str(read_infos[1].start_pos_ref)
                                read_hits[1][5] = read_infos[1].new_CIGAR
                            
                            if read_infos[0].is_unmapped == False and read_infos[1].is_unmapped == False:
                                read_hits[0][7] = str(read_hits[1][3])
                                read_hits[1][7] = str(read_hits[0][3])

                                if int(read_hits[0][8]) > 0:
                                    
                                    distance = read_infos[1].start_pos_ref + read_infos[1].n_bases_read  - read_infos[0].start_pos_ref + 1
                                    
                                    read_hits[0][8] = str( distance)
                                    read_hits[1][8] = str(-distance)
                                else:
                                    
                                    distance = read_infos[0].start_pos_ref + read_infos[0].n_bases_read  - read_infos[1].start_pos_ref
                                    
                                    read_hits[0][8] = str(-distance)
                                    read_hits[1][8] = str( distance)
                                
                    else:
                        mate1_hits = []
                        mate2_hits = []
                        
                        for i in range(n_hits):
                            FLAG = int(read_hits[i][1])
                            
                            if is_mate_1(FLAG)   and i not in mate1_hits:
                                mate1_hits.append(i)

                            elif is_mate_2(FLAG) and i not in mate2_hits:
                                mate2_hits.append(i)
                        
                        max_score    = 0
                        pairs_scores = dict()
                        
                        for i in mate1_hits: 
                            for j in mate2_hits:
                                if is_pair(read_hits[i], read_hits[j]):
                                    
                                    score = read_infos[i].read_score + read_infos[j].read_score
                                    chrom = read_hits[i][2]
                                    # print(read_hits[i])
                                    # print(read_hits[j])
                                    # pairs_lengths[len(read_hits[i][9]) + len(read_hits[j][9])]  += 1
                                    
                                    if not are_mates_from_same_transcript(chrom, 
                                                                            read_infos[i].is_novel, read_infos[j].is_novel, 
                                                                            read_infos[i].regions, read_infos[j].regions, 
                                                                            regions_relations_table, regions_transcripts_table, jumps_transcripts_table):
                                        score -= novel_transcript_penalty
                                    
                                    if score > max_score:
                                        max_score = score 
                                        
                                    pairs_scores[str(i) + '__' + str(j)] = score
                                    
                        mate_1_good_hits = set()
                        mate_2_good_hits = set()
                        
                        for pair, score in pairs_scores.items():
                            if score >= max_score - max_score_margin:
                                i, j = pair.split("__")
                                mate_1_good_hits.add(i)
                                mate_2_good_hits.add(j)
                        
                        MAPQ_i = get_MAPQ(mate_1_good_hits)
                        MAPQ_j = get_MAPQ(mate_2_good_hits)
                        
                        primary_i = -1
                        primary_j = -1
                        is_yet_to_set_primary = True
                        
                        for pair, score in pairs_scores.items():
                            i, j = pair.split("__")
                            i = int(i)
                            j = int(j)
                            
                            if read_infos[i].is_clipped or read_infos[j].is_clipped:
                                
                                read_hits[i][3] = str(read_infos[i].start_pos_ref)
                                read_hits[i][5] = read_infos[i].new_CIGAR

                                read_hits[j][3] = str(read_infos[j].start_pos_ref)
                                read_hits[j][5] = read_infos[j].new_CIGAR
                            
                                read_hits[i][7] = str(read_hits[j][3])
                                read_hits[j][7] = str(read_hits[i][3])
                                
                                if int(read_hits[i][8]) > 0:
                                    
                                    distance = read_infos[j].start_pos_ref + read_infos[j].n_bases_read  - read_infos[i].start_pos_ref + 1
                                    
                                    read_hits[i][8] = str( distance)
                                    read_hits[j][8] = str(-distance)
                                else:
                                    
                                    distance = read_infos[i].start_pos_ref + read_infos[i].n_bases_read  - read_infos[j].start_pos_ref
                                    read_hits[i][8] = str(-distance)
                                    read_hits[j][8] = str( distance)
                            
                            if score >= max_score - max_score_margin:
                                if read_infos[i].is_bad_read == False:
                                    read_hits[i][4] = str(MAPQ_i)
                                if read_infos[j].is_bad_read == False:
                                    read_hits[j][4] = str(MAPQ_j)
                                
                                if is_yet_to_set_primary:
                                    read_hits[i][1] = str(set_primary_alignment_flag(int(read_hits[i][1])))
                                    read_hits[j][1] = str(set_primary_alignment_flag(int(read_hits[j][1])))
                                    primary_i = i
                                    primary_j = j
                                    is_yet_to_set_primary = False
                                else:
                                    if i != primary_i:
                                        read_hits[i][1] = str(set_secondary_alignment_flag(int(read_hits[i][1])))
                                    if j != primary_j:
                                        read_hits[j][1] = str(set_secondary_alignment_flag(int(read_hits[j][1])))
                            else :
                                if i not in mate_1_good_hits:
                                    read_hits[i][4] = "0"
                                    read_hits[i][1] = str(set_secondary_alignment_flag(int(read_hits[i][1])))
                                if j not in mate_2_good_hits:
                                    read_hits[j][4] = "0"
                                    read_hits[j][1] = str(set_secondary_alignment_flag(int(read_hits[j][1])))
                                    
                            
                        # if MAPQ_i == 40 and MAPQ_j == 40:
                        if is_verbose:
                            print("-----------------------------------------------------")
                            print("max_score {} ".format(max_score))
                            
                            print("mate1_hits:", mate1_hits)
                            print("mate2_hits:", mate2_hits)

                            print("Mate1 new MAPQ:", mate_1_good_hits,  MAPQ_i)
                            print("Mate2 new MAPQ:", mate_2_good_hits,  MAPQ_j)

                            print("\nWinners:")
                            for pair, score in pairs_scores.items():
                                if score >= max_score - max_score_margin:
                                    i, j = pair.split("__")
                                    i = int(i)
                                    j = int(j)
                                    print(pair, score)
                                    print(read_hits[i])
                                    print(read_infos[i].jt_tag, read_infos[i].js_tag, read_infos[i].new_CIGAR)
                                    print()
                                    print(read_hits[j])
                                    print(read_infos[j].jt_tag, read_infos[j].js_tag, read_infos[j].new_CIGAR)
                                    print()
                            
                            print("\nLosers:")
                            for pair, score in pairs_scores.items():
                                if score < max_score - max_score_margin:
                                    i, j = pair.split("__")
                                    i = int(i)
                                    j = int(j)
                                    print(pair, score)
                                    print(read_hits[i])
                                    print(read_infos[i].jt_tag, read_infos[i].js_tag, read_infos[i].new_CIGAR)
                                    print()
                                    print(read_hits[j])
                                    print(read_infos[j].jt_tag, read_infos[j].js_tag, read_infos[j].new_CIGAR)
                                    print()
                    
                    for i in range(n_hits):
                        if read_hits[i][5] == '':
                            # print('Reads with empty CIGAR')
                            # print(read_hits[i])
                            # print(read_infos[i])
                            read_hits[i][5] = '*'

                        if read_infos[i].n_junctions != 0:
                            read_hits[i].append(read_infos[i].nj_tag)
                            read_hits[i].append(read_infos[i].jt_tag)
                            read_hits[i].append(read_infos[i].js_tag)

                        sam_output_record = "\t".join([str(field) for field in read_hits[i] ]) + "\n"
                        
                        samtools_process.stdin.write(sam_output_record.encode())
                else:
                    # we are not interested in these reads for now ! 
                    for i in range(n_hits):
                        sam_output_record = "\t".join([str(field) for field in read_hits[i] ]) + "\n"
                        samtools_process.stdin.write(sam_output_record.encode())
                        
                # add next read 
                pre_read_id = read_id
                read_hits   = [fields]

            if is_end_of_file:
                break

    # Close the subprocess stdin, wait for completion, and handle errors
    samtools_process.stdin.close()
    samtools_process.wait()

    message = ""

    if samtools_process.returncode != 0:
        message += "Error occurred while using samtools.\n"
    else:
        message += "[{}]\t[LOG]\tFinished writing BAM successfully into '{}'\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), output_bam)
    
        message += "[{}]\t[LOG]\tNumber of SAM records: {} \n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n_sam_records)
        message += "[{}]\t[LOG]\tNumber of reads IDs:   {} \n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n_reads_ids)
        message += "[{}]\t[LOG]\tNumber of processed reads IDs: {}  {:.2f}% \n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n_processed_reads_ids, n_processed_reads_ids * 100 / n_reads_ids)
    
    return message 

def generate_made_up_junctions_sequences(junctions_table, sequence_table, n_junctions, window, canonical_percent=99, random_seed=777, key_prefix="", add_junction_seq=False):
    
    print("[{}]\t[LOG]\tStarted fabricating made-up splice junctions from the unmasked regions of the genome:".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print("[{}]\t[LOG]\tStarted masking the genome".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    random.seed(random_seed)
    
    w = window // 2
    sequence_np_table = dict()

    for key in sequence_table.keys():
        sequence_np_table[key] = np.array(list(sequence_table[key]), dtype=str)
    
    for key in junctions_table.keys():
        chrom = junctions_table[key].chrom
        donor_end = junctions_table[key].donor_end
        acceptor_start = junctions_table[key].acceptor_start
        
        donor_sequnce_start = donor_end - window
        donor_sequnce_end = donor_end + window
        
        acceptor_sequnce_start = acceptor_start - window
        acceptor_sequnce_end = acceptor_start + window
        
        if donor_sequnce_start < 0:
            donor_sequnce_start = 0
        
        if acceptor_sequnce_end >= len(sequence_table[chrom]):
            acceptor_sequnce_end = len(sequence_table[chrom]) - 1
            
        sequence_np_table[chrom][donor_sequnce_start : donor_sequnce_end] = 'N'
        sequence_np_table[chrom][acceptor_sequnce_start : acceptor_sequnce_end] = 'N'
    
    
    chroms      = list(sequence_table.keys())
    n_chroms    = len(chroms)
    n           = 0

    fake_junctions_table = dict()
    
    print("[{}]\t[LOG]\tStarted fabricating {} made-up junctions".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), int(n_junctions)))
    
    while n < n_junctions:
        random_chrom = random.randint(0, n_chroms - 1)
        chrom        = chroms[random_chrom]
        
        if len(sequence_table[chrom]) < 8 * w:
            continue
        
        donor_pos = random.randint(2 * w, len(sequence_table[chrom]) - (2 * w))
        donor_seq = "".join(list(sequence_np_table[chrom][donor_pos - w : donor_pos + w]))
        
        dna_alphabet = [base for base in "ACTGactg"]
        
        signal = "GTAG"
        is_canonical_signal = True 
        signal_random = random.randint(0, 100)
        
        if signal_random > canonical_percent:
            is_canonical_signal = False
            
        if (is_canonical_signal == False or donor_seq[w : w+2] == signal[:2]) and all(base in dna_alphabet for base in donor_seq):
            n_attempts = 0
            while n_attempts < 10:
                n_attempts += 1  
                
                try:
                    acceptor_pos = random.randint(donor_pos + (2 * w), min(donor_pos + 5000, len(sequence_table[chrom]) - (2 * w)))
                except ValueError:
                    break
                
                acceptor_seq = "".join(list(sequence_np_table[chrom][acceptor_pos - w : acceptor_pos + w]))

                if (is_canonical_signal == False or acceptor_seq[w -2:w] == signal[2:]) and all(base in dna_alphabet for base in acceptor_seq):  
                    
                    junction_id =  key_prefix + chrom + "__+__" + str(donor_pos) + "__" + str(acceptor_pos)
                    
                    if junction_id not in fake_junctions_table.keys():
                        if add_junction_seq == True:
                            junction_sequence = donor_seq [:w] + acceptor_seq[w:]
                        else:
                            junction_sequence = ''

                        fake_junctions_table[junction_id] = Junction(chrom, 
                                                                    "", "", 
                                                                    donor_pos, acceptor_pos, 
                                                                    "+", "GTAG", "FAKE", 
                                                                    donor_seq.upper(), acceptor_seq.upper(), junction_sequence.upper())
                        
                        n += 1
                        if n % 100000 == 0:
                            print("[{}]\t[LOG]\tDone making {} junctions".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n))
                        break
                    
    return fake_junctions_table

def write_junctions_info(found_junctions:dict, 
                        junctions_scores:dict,
                        path_out:str, 
                        output_prefix:str):

    with open(path_out + output_prefix + "_junctions.tsv", "w") as file:
        
        file.write("ID\tSignla\tType\tScore\tDonor\tAcceptor\n")

        for splice_id in found_junctions.keys():
            ID=splice_id
            # Chr=found_junctions[splice_id].chrom
            # Strand=found_junctions[splice_id].strand
            # DonorEnd=found_junctions[splice_id].donor_end
            # AcceptorStart=found_junctions[splice_id].acceptor_start
            Signal=found_junctions[splice_id].signal
            Donor=found_junctions[splice_id].donor_kmers
            Acceptor=found_junctions[splice_id].acceptor_kmers
            # Junction=found_junctions[splice_id].junction_kmers
            Type=found_junctions[splice_id].transcript_type
            Score=str(junctions_scores[splice_id])

            file.write(f"{ID}\t{Signal}\t{Type}\t{Score}\t{Donor}\t{Acceptor}\n")

def collect_junctions_from_TP_FP(sequences_table: dict, 
                                 path: str, 
                                 window: int, 
                                 is_tp=True):
    
    junctions_table = dict()
    w = window // 2

    with open(path, "r") as file:
        n_reads = 0
        n_spliced_reads = 0

        while True:
            line = file.readline()

            if len(line) == 0 :
                break
            elif(line[0] == '@'):
                continue
            else:
                fields  = line[:-1].split('\t')

                if is_tp:
                    fields = fields[1:]
                
                FLAG    = int(fields[1])
                chrom   = fields[2]
                start   = int(fields[3])
                MAPQ    = int(fields[4])
                CIGAR   = fields[5]
                strand  = get_strand(FLAG, fields)
                n_bases = ""
                locus   = start - 1 # 0-based
                n_reads += 1

                nNs = CIGAR.count('N')

                if nNs != 1 or 'js:' not in line:
                    continue
            
                n_spliced_reads += 1
                
                for character in CIGAR:
                    if character.isdigit():
                        n_bases += character
                    elif character in ['M', '=', 'X']:
                        locus += int(n_bases)
                        n_bases = ""
                    elif character == 'D':
                        locus += int(n_bases)
                        n_bases = ""
                    elif character in ['I', 'S', 'H']:
                        n_bases = ""
                    elif character == 'N':
                        donor_site = locus
                        locus     += int(n_bases)
                        n_bases    = ""
                        acceptor_site = locus
                        
                        splice_id       = chrom + "__" + strand + "__" + str(donor_site) + "__" + str(acceptor_site)
                        transcript_type = ""

                        if "Annotated" in line:
                            transcript_type = "Annotated"
                        else:
                            transcript_type = "Novel"
                        
                        pattern_js = r'js:Z:\d+'
                        pattern_rs = r'rs:Z:\d+'

                        matches_js = re.findall(pattern_js, line)
                        matches_rs = re.findall(pattern_rs, line)

                        score_js = int(matches_js[0].split(':')[2])
                        score_rs = 0

                        if len(matches_rs) == 1:
                            score_rs = int(matches_rs[0].split(':')[2])

                        if splice_id not in junctions_table.keys():

                            if acceptor_site - donor_site > window:
                                new_w = w
                            elif acceptor_site - donor_site > 6:
                                new_w = (acceptor_site - donor_site - 1) // 2
                            else:
                                new_w = 3

                            donor_sequence_start = donor_site - new_w
                            donor_sequence_end   = donor_site + new_w
                            
                            acceptor_sequence_start = acceptor_site - new_w
                            acceptor_sequence_end   = acceptor_site + new_w
                            
                            if donor_sequence_start < 0:
                                donor_sequence_start = 0
                            
                            if donor_sequence_end >= len(sequences_table[chrom]):
                                donor_sequence_end = len(sequences_table[chrom]) - 1
                            
                            if acceptor_sequence_start < 0:
                                acceptor_sequence_start = 0
                            
                            if acceptor_sequence_end >= len(sequences_table[chrom]):
                                acceptor_sequence_end = len(sequences_table[chrom]) - 1
                            
                            if donor_site >= len(sequences_table[chrom]) or acceptor_site >= len(sequences_table[chrom]):
                                print("[Error]\tProblematic record, locus outside of chromosom\nDonor site: {}\nAcceptor site: {}\nRecord line: {}".format(donor_site, acceptor_site, line))
                                exit()
                            
                            signal            = sequences_table[chrom][donor_site: donor_site + 2] + sequences_table[chrom][acceptor_site - 2: acceptor_site]
                            junction_sequence = sequences_table[chrom][donor_sequence_start    : donor_site] + sequences_table[chrom][acceptor_site: acceptor_sequence_end]
                            donor_sequence    = sequences_table[chrom][donor_sequence_start    : donor_sequence_end]
                            acceptor_sequence = sequences_table[chrom][acceptor_sequence_start : acceptor_sequence_end]
                            
                            signal            = signal.upper()
                            junction_sequence = junction_sequence.upper()
                            donor_sequence    = donor_sequence.upper()
                            acceptor_sequence = acceptor_sequence.upper()
                            
                            if "N" in signal or  "N" in junction_sequence or "N" in donor_sequence or "N" in acceptor_sequence:
                                print("[Warning]\tProblematic record due to N:\nN replaced with A:\nDonor site: {}\nAcceptor site: {}\nRecord line: {}".format(donor_site, acceptor_site, line))
                                print("signal: {}".format(signal))
                                print("junction_sequence: {}".format(junction_sequence))
                                print("donor_sequence: {}".format(donor_sequence))
                                print("acceptor_sequence: {}\n".format(acceptor_sequence))
                                
                                signal            = signal.replace('N', 'A')
                                junction_sequence = junction_sequence.replace('N', 'A')
                                donor_sequence    = donor_sequence.replace('N', 'A')
                                acceptor_sequence = acceptor_sequence.replace('N', 'A')
                                
                            if strand == "-":
                                signal            = reverse_complement(signal)
                                junction_sequence = reverse_complement(junction_sequence)
                                temp              = donor_sequence
                                donor_sequence    = reverse_complement(acceptor_sequence)
                                acceptor_sequence = reverse_complement(temp)
                            
                            junctions_table[splice_id] = JunctionFromRead( 
                                                        signal, 
                                                        transcript_type,
                                                        donor_sequence, acceptor_sequence, junction_sequence, 
                                                        fields[9], score_js, score_rs)

    return junctions_table

def collect_junctions_from_SAM(GTF_junctions:dict, sequences_table: dict, 
                               path: str, 
                            #    k: int, 
                               window: int, 
                               samtools:str): # junctions_path: str, kmers_path: str, w):
    
    print("[{}]\t[LOG]\tCollecting splice junctions from SAM/BAM file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), path))
    if path[-4:] == ".bam":
        records = subprocess.Popen([samtools, "view", path], stdout=subprocess.PIPE)
    else:
        records = subprocess.Popen(["cat", path], stdout=subprocess.PIPE)
    
    n_reads         = 0
    n_spliced_reads = 0
    junctions_table = dict()
    min_sequence    = 8
    w       = window // 2
    XS_A_re = 0
    XS_A_fr = 0
    
    while True:  
        line = records.stdout.readline().decode()
        if len(line) == 0 :
            break
        elif(line[0] == '@'):
            continue
        else:
            """6:1,13960,+,100M;1,14094,-,100M:AAAAAAAAAAAN/877P       99      chr1    13961   40      79M     =       14132   234     CCACCTTCTTCCTGAGTCATTCCTGCAGCCTTGCTCCCTAACCTGCCCCACAGCCTTGCCTGGATTTCTATCTCCCTGG IBB82=@BBBFIIIIIIIIIIIIIIIFHFDDIIIIIIHIHIIHIIIIIIIEGDGGHGGGGG@CCE38<C8AGHIIGG@@
            """
            fields  = line[:-1].split('\t')
            FLAG    = int(fields[1])
            chrom   = fields[2]
            start   = int(fields[3])
            MAPQ    = int(fields[4])
            CIGAR   = fields[5]
            strand  = get_strand(FLAG, fields)
            n_bases = ""
            locus   = start - 1 # 0-based
            n_reads += 1
            
            if "N" not in CIGAR:
                continue
            
            n_spliced_reads += 1
            
            for character in CIGAR:
                if character.isdigit():
                    n_bases += character
                elif character in ['M', '=', 'X']:
                    locus += int(n_bases)
                    n_bases = ""
                elif character == 'D':
                    locus += int(n_bases)
                    n_bases = ""
                elif character in ['I', 'S', 'H']:
                    n_bases = ""
                elif character == 'N':
                    donor_site = locus
                    locus     += int(n_bases)
                    n_bases    = ""
                    acceptor_site = locus
                    
                    splice_id       = chrom + "__" + strand + "__" + str(donor_site) + "__" + str(acceptor_site)
                    transcript_type = ""
                    
                    if splice_id in GTF_junctions.keys():
                        transcript_type = GTF_junctions[splice_id].transcript_type
                    else:
                        transcript_type = "Novel"
                    
                    if splice_id not in junctions_table.keys():

                        if acceptor_site - donor_site > window:
                            new_w = w
                        elif acceptor_site - donor_site > 6:
                            new_w = (acceptor_site - donor_site - 1) // 2
                        else:
                            new_w = 3
                        
                        donor_sequence_start = donor_site - new_w
                        donor_sequence_end   = donor_site + new_w
                        
                        acceptor_sequence_start = acceptor_site - new_w
                        acceptor_sequence_end   = acceptor_site + new_w
                        
                        if donor_sequence_start < 0:
                            donor_sequence_start = 0

                        if chrom not in sequences_table.keys():
                            print("[{}]\t[Error]\tRead chromosome name {} not found in the reference chromosomes {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), chrom, sequences_table.keys()))
                            exit()

                        if donor_sequence_end >= len(sequences_table[chrom]):
                            donor_sequence_end = len(sequences_table[chrom]) - 1
                        
                        if acceptor_sequence_start < 0:
                            acceptor_sequence_start = 0
                        
                        if acceptor_sequence_end >= len(sequences_table[chrom]):
                            acceptor_sequence_end = len(sequences_table[chrom]) - 1
                        
                        if donor_site >= len(sequences_table[chrom]) or acceptor_site >= len(sequences_table[chrom]):
                            print("[Error]\tProblematic record, locus outside of chromosom\nDonor site: {}\nAcceptor site: {}\nRecord line: {}".format(donor_site, acceptor_site, line))
                            exit()
                        
                        signal            = sequences_table[chrom][donor_site: donor_site + 2] + sequences_table[chrom][acceptor_site - 2: acceptor_site]
                        junction_sequence = sequences_table[chrom][donor_sequence_start    : donor_site] + sequences_table[chrom][acceptor_site: acceptor_sequence_end]
                        donor_sequence    = sequences_table[chrom][donor_sequence_start    : donor_sequence_end]
                        acceptor_sequence = sequences_table[chrom][acceptor_sequence_start : acceptor_sequence_end]
                        
                        signal            = signal.upper()
                        junction_sequence = junction_sequence.upper()
                        donor_sequence    = donor_sequence.upper()
                        acceptor_sequence = acceptor_sequence.upper()
                        
                        if "N" in signal or  "N" in junction_sequence or "N" in donor_sequence or "N" in acceptor_sequence:
                            print("[Warning]\tProblematic record due to N:\nN replaced with A:\nDonor site: {}\nAcceptor site: {}\nRecord line: {}".format(donor_site, acceptor_site, line))
                            print("signal: {}".format(signal))
                            print("junction_sequence: {}".format(junction_sequence))
                            print("donor_sequence: {}".format(donor_sequence))
                            print("acceptor_sequence: {}\n".format(acceptor_sequence))
                            
                            signal            = signal.replace('N', 'A')
                            junction_sequence = junction_sequence.replace('N', 'A')
                            donor_sequence    = donor_sequence.replace('N', 'A')
                            acceptor_sequence = acceptor_sequence.replace('N', 'A')
                            
                        if strand == "-":
                            signal            = reverse_complement(signal)
                            junction_sequence = reverse_complement(junction_sequence)
                            temp              = donor_sequence
                            donor_sequence    = reverse_complement(acceptor_sequence)
                            acceptor_sequence = reverse_complement(temp)
                        
                        if strand   == '-':
                            XS_A_re += 1 
                        elif strand == '+':
                            XS_A_fr += 1
                        
                        junctions_table[splice_id] = Junction( chrom, " ", " ",  
                                                    int(donor_site), 
                                                    int(acceptor_site), 
                                                    strand,
                                                    signal, 
                                                    transcript_type,
                                                    donor_sequence, acceptor_sequence, junction_sequence)

                        if len(signal) != 4 or len(junctions_table[splice_id].acceptor_kmers) == 0 or len(junctions_table[splice_id].donor_kmers) == 0:
                            print(acceptor_sequence_end, len(sequences_table[chrom]))
                            print(line)
                            print(junction_sequence)
                            print(donor_sequence)
                            print(acceptor_sequence)
                            
                            print("Signal:\t{}".format(junctions_table[splice_id].signal))
                            print("Donor:\t{}".format(junctions_table[splice_id].donor_kmers))
                            print("Acceptor:\t{}".format(junctions_table[splice_id].acceptor_kmers))
                            
                            print(junctions_table[splice_id]())
                            exit()

    print("[{}]\t[INFO]\tSense junctions {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), XS_A_fr))
    print("[{}]\t[INFO]\tAntisense junctions {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), XS_A_re))
    print("[{}]\t[INFO]\tTotal number of reads {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n_reads))
    print("[{}]\t[INFO]\tTotal number of spliced reads {} {}%".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), n_spliced_reads, n_spliced_reads/n_reads * 100))
    
    return junctions_table

def addition_scaling_rounding(donor_prob, acceptor_prob):
    return round(50 * (donor_prob + acceptor_prob))

def multiplication_scaling_rounding(donor_prob, acceptor_prob):
    return round(100 * donor_prob * acceptor_prob)

def negative_log_10(donor_prob, acceptor_prob):
    return round(50 * (-1) * math.log10((1 - donor_prob) * (1 - acceptor_prob)))

def write_prediction_datapoints(junctions_table:dict, output_path:str):
    
    print("[{}]\t[LOG]\tWritting dev.csv file for predicting into '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), output_path))
    
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    os.makedirs(output_path)
    
    eval_file = open(output_path + "/dev.csv",   "w")
    eval_file.write("sequence,label\n")

    ## Adding datapoints from junctions_table:
    for junction_id in junctions_table.keys():
        donor_seq    = junctions_table[junction_id].donor_kmers 
        acceptor_seq = junctions_table[junction_id].acceptor_kmers
            
        eval_file.write("{},0\n{},1\n".format(donor_seq, acceptor_seq)) 

    eval_file.write("NNNNNNNNNNNN,0\n") # since all the models are tuned on three labels
    eval_file.write("NNNNNNNNNNNN,1\n") # since all the models are tuned on three labels
    eval_file.write("NNNNNNNNNNNN,2\n") # since all the models are tuned on three labels
    
    eval_file.close()
    
    print("[{}]\t[LOG]\tdev.csv file contains:   0: {}, 1: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),   len(junctions_table.keys()), len(junctions_table.keys())))

def write_prediction_datapoints_from_reads(path_sam_in:str, 
                                           path_out:str,
                                           output_prefix:str,
                                           samtools:str,
                                           predictions_output_dir:str, 
                                           n_reads:int):
    print("[{}]\t[LOG]\tWritting dev.csv file for predicting reads into {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), predictions_output_dir))
    
    if os.path.isdir(predictions_output_dir):
        shutil.rmtree(predictions_output_dir)
        
    os.makedirs(predictions_output_dir)
    
    eval_file = open(predictions_output_dir + "/dev.csv",   "w")
    eval_file.write("sequence,label\n")

    path_sam_in = "{}/{}.bam".format(path_out, output_prefix)

    records     = subprocess.Popen([samtools, "view", "-h", path_sam_in], stdout=subprocess.PIPE)
    n = 0

    while True:
        line = records.stdout.readline().decode()

        if len(line) == 0 or n >= n_reads:
            break
        elif line[0] == '@':
            continue
        
        n += 1
        read_seq = line.split('\t')[9]
        eval_file.write("{},2\n".format(read_seq)) 

    eval_file.write("NNNNNNNNNNNN,0\n") # since all the models are tuned on three labels
    eval_file.write("NNNNNNNNNNNN,1\n") # since all the models are tuned on three labels
    eval_file.write("NNNNNNNNNNNN,2\n")
    eval_file.close()
    
    print("[{}]\t[LOG]\tdev.csv file contains:   0: {}, 1: {}, 2:{}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 0, 0, n))

    return n + 3

def add_reads_scores_to_sam(path_sam_in:str, 
                            path_out:str,
                            output_prefix,
                            samtools:str,
                            output_dir:str, 
                            n_reads:int, 
                            method="Add"):
    
    print("[{}]\t[LOG]\tWritting reads scores to SAM".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    probs = np.load(output_dir + '/probs.npy')
    
    path_sam_in  = "{}/{}.bam".format(path_out, output_prefix)
    path_sam_out = "{}/{}_with_reads.bam".format(path_out, output_prefix)

    records            = subprocess.Popen([samtools, "view", "-h",        path_sam_in], stdout=subprocess.PIPE)
    samtools_process   = subprocess.Popen([samtools, "view", "-bS", "-o", path_sam_out , "-"], stdin=subprocess.PIPE)

    n = 0

    while True:
        line = records.stdout.readline().decode()

        if len(line) == 0:
            break
        elif line[0] == '@' or n >= n_reads:
            samtools_process.stdin.write(line.encode())
            continue
        
        prob_donor    = probs[n][0]
        prob_acceptor = probs[n][1]
        
        if method=="Add":
            score_read = addition_scaling_rounding(prob_donor, prob_acceptor)
        else:
            score_read = multiplication_scaling_rounding(prob_donor, prob_acceptor)

        # score_read  = int(100 * (1 - probs[n][2]))
        n         += 1

        line = line[0:-1] + "\trs:i:{}\n".format(score_read) 
        samtools_process.stdin.write(line.encode())

    # Close the subprocess stdin, wait for completion, and handle errors
    samtools_process.stdin.close()
    samtools_process.wait()

    os.remove(path_sam_in) 

def write_tuning_datapoints(junctions_table:dict, fake_junctions_table:dict, output_path:str, max_n_junctions:int,
                           train_percent=70, eval_percent=15, factor=1.0, shuffle_csv=True, random_seed=777):
    
    print("[{}]\t[LOG]\tWritting train.csv, dev.csv and test.csv files for tuning into {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), output_path))
    
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    
    train_file = open(output_path + "/train.csv", "w")
    eval_file   = open(output_path + "/dev.csv",   "w")
    test_file  = open(output_path + "/test.csv",  "w")
    
    train_file.write("sequence,label\n")
    eval_file.write("sequence,label\n")
    test_file.write("sequence,label\n")
    
    random.seed(random_seed) # to be able to reproduce the experiment
    
    donors_set   = set() # to avoid any shared datapoints between the three datasets
    acceptors_set = set() 
    
    counts = [0, 0, 0, 0, 0, 0, 0, 0, 0]   # [train_0, train_1, train_2, eval_0, eval_1, eval_2, test_0, test_1, test_2]
    
    ## Adding datapoints from junctions_table:
    n_junctions = 0 # to count how many junctions are written
    junctions_ids_list = list(junctions_table.keys())
    random.shuffle(junctions_ids_list)
    
    for junction_id in junctions_ids_list:
        if n_junctions >= max_n_junctions:
            break 
        
        donor_seq    = junctions_table[junction_id].donor_kmers 
        acceptor_seq = junctions_table[junction_id].acceptor_kmers

        if donor_seq in donors_set or acceptor_seq in acceptors_set:
            continue 
        
        n_junctions += 1
        
        donors_set.add(donor_seq)
        acceptors_set.add(acceptor_seq)
        
        random_num = random.randint(0, 100)
        
        if (random_num >= train_percent + eval_percent):
            file = test_file 
            counts[6] += 1
            counts[7] += 1
            
        elif (random_num >= train_percent):
            file = eval_file
            counts[3] += 1
            counts[4] += 1
        else:
            file = train_file
            counts[0] += 1
            counts[1] += 1
            
        file.write("{},0\n{},1\n".format(donor_seq, acceptor_seq)) 
    
    fake_donors_set    = set() # to avoid any shared datapoints between the three datasets and to avoid fake junctions that is reported in the junctions_table
    fake_acceptors_set = set() 
    
    ## Adding datapoints from fake_junctions_table:
    
    negative_max_n_junctions = n_junctions * factor
    n_junctions     = 0
    
    target_is_reached = False
    
    for junction_id in fake_junctions_table.keys():
        if n_junctions >= negative_max_n_junctions:
            target_is_reached = True
            break
        
        fake_donors_seq   = fake_junctions_table[junction_id].donor_kmers 
        fake_acceptos_seq = fake_junctions_table[junction_id].acceptor_kmers
        
        if fake_donors_seq in fake_donors_set or fake_donors_seq in fake_acceptors_set or fake_donors_seq in donors_set  or fake_donors_seq in acceptors_set:
            continue
        
        if fake_acceptos_seq in fake_acceptors_set or fake_acceptos_seq in fake_donors_set or fake_acceptos_seq in donors_set or fake_acceptos_seq in acceptors_set:
            continue
        
        n_junctions += 1
        
        fake_donors_set.add(fake_donors_seq)
        fake_acceptors_set.add(fake_acceptos_seq)
        
        random_num = random.randint(0, 100)
        
        if (random_num >= train_percent + eval_percent):
            file = test_file 
            counts[8] += 2
            
        elif (random_num >= train_percent):
            file = eval_file
            counts[5] += 2
        else:
            file = train_file
            counts[2] += 2
        
        file.write("{},2\n{},2\n".format(fake_donors_seq, fake_acceptos_seq)) 
        
    train_file.close()
    eval_file.close()
    test_file.close()
    
    if target_is_reached:
        print("[{}]\t[LOG]\tThe target number of made-up junctions for the negative dataset is reached successfully, negative_max_n_junctions={}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), negative_max_n_junctions))
    else:
        print("[{}]\t[Error]\tThe target number of made-up junctions for the negative dataset was not reached, negative_max_n_junctions={}, please increase n_made_up_junctions".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), negative_max_n_junctions))
        exit()
        
    if shuffle_csv == True:
        # Shuffle train
        df = pd.read_csv(output_path + "/train.csv")
        shuffled_df = df.sample(frac=1, random_state=random_seed)
        shuffled_df.to_csv(output_path + "/train.csv", index=False)

        # Shuffle dev
        df = pd.read_csv(output_path + "/dev.csv")
        shuffled_df = df.sample(frac=1, random_state=random_seed)
        shuffled_df.to_csv(output_path + "/dev.csv", index=False)
        
        # Shuffle test
        df = pd.read_csv(output_path + "/test.csv")
        shuffled_df = df.sample(frac=1, random_state=random_seed)
        shuffled_df.to_csv(output_path + "/test.csv", index=False)
    
    print("[{}]\t[LOG]\ttrain.csv file contains: 0: {}, 1: {}, 2:{}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), counts[0], counts[1], counts[2]))
    print("[{}]\t[LOG]\tdev.csv file contains:   0: {}, 1: {}, 2:{}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),   counts[3], counts[4], counts[5]))
    print("[{}]\t[LOG]\ttest.csv file contains:  0: {}, 1: {}, 2:{}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),  counts[6], counts[7], counts[8]))

def get_junctions_scores(path, junctions_table, method="Add"):
    probs = np.load(path + '/probs.npy')
    junctions_scores = dict()

    i = 0
    for splice_id in junctions_table.keys():
        prob_donor    = probs[i    ][0]
        prob_acceptor = probs[i + 1][1]
        i += 2
        
        if method=="Add":
            junctions_scores[splice_id] = addition_scaling_rounding(prob_donor, prob_acceptor)
        else:
            junctions_scores[splice_id] = multiplication_scaling_rounding(prob_donor, prob_acceptor)
    
    return junctions_scores


class DictBatchIterator:
    def __init__(self, dictionary, batch_size):
        self.dictionary = dictionary
        self.batch_size = batch_size
        self.keys = list(dictionary.keys())
        self.values = list(dictionary.values())
        self.num_batches = (len(self.keys) + batch_size - 1) // batch_size
        self.current_batch = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_batch < self.num_batches:
            start_idx = self.current_batch * self.batch_size
            end_idx = min((self.current_batch + 1) * self.batch_size, len(self.keys))
            batch_keys = self.keys[start_idx:end_idx]
            batch_values = self.values[start_idx:end_idx]
            batch_dict = {k: v for k, v in zip(batch_keys, batch_values)}
            self.current_batch += 1
            return batch_dict
        else:
            raise StopIteration


import subprocess



def run_gsnap_alignment(gsnap_aln_flags, threads, gsnap_index, gsnap_idx_flags, mate_1, mate_2, sm_view_flags, sam_sort_flags, output_bam, output_prefix, log_file):
    """
    Runs GSNAP alignment with the specified parameters and saves output to a BAM file.
    
    Parameters:
    - gsnap_aln_flags (str): Additional flags for GSNAP alignment.
    - threads (int): Number of threads to use.
    - gsnap_index (str): Directory of the GSNAP index.
    - mate_1 (str): Path to the first read file.
    - mate_2 (str): Path to the second read file.
    - sm_view_flags (str): Flags for `samtools view`.
    - sam_sort_flags (str): Flags for `samtools sort`.
    - output_bam (str): Path to the output BAM file.
    - output_prefix (str): Prefix for temporary sorting files.
    - log_file (str): Path to the log file for capturing output and error messages.
    """
    
    cmd_gsnap = (
        f"gsnap {gsnap_aln_flags} -t {threads} -D {gsnap_index} {gsnap_idx_flags} "
        f"{mate_1} {mate_2} | "
        f"samtools view -u -b -h {sm_view_flags} - | "
        f"samtools sort {sam_sort_flags} -o {output_bam} -T {output_prefix} -@ {threads}"
    )
    
    # with open(log_file, "a") as log:
    #     try:
    #         subprocess.run(cmd_gsnap, shell=True, stdout=log, stderr=log, check=True)
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error occurred during GSNAP alignment: {e}")
    #         exit(1)
    with open(log_file, "a") as log:
        try:
            # Run the command and capture the output
            result = subprocess.run(
                cmd_gsnap,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            # Write stdout and stderr to the log file
            log.write(result.stdout)
            log.write(result.stderr)

            # Print stdout and stderr to the console
            # print(result.stdout)
            # print(result.stderr)

        except subprocess.CalledProcessError as e:
            print(f"Error occurred during GSNAP alignment: {e}")
            if e.stdout:
                print("STDOUT:")
                print(e.stdout)
                log.write("STDOUT:\n")
                log.write(e.stdout)
            if e.stderr:
                print("STDERR:")
                print(e.stderr)
                log.write("STDERR:\n")
                log.write(e.stderr)
            exit(1)

def load_folder_from_encrypted_file(encrypted_file, encryption_key, output_directory):
    """
    Decrypt an encrypted file and extract its contents as a folder.

    Args:
        encrypted_file (str): Path to the encrypted file.
        encryption_key (str): Encryption key used to encrypt the file.
        output_directory (str): Path to extract the contents of the decrypted file.

    Returns:
        str: Path to the extracted folder.
    """
    # Initialize the Fernet cipher with the provided key
    cipher_suite = Fernet(encryption_key.encode())

    # Decrypt the encrypted file
    decrypted_tar_path = "/tmp/decrypted_model.tar.gz"  # Temporary path for the decrypted .tar.gz
    with open(encrypted_file, "rb") as file:
        encrypted_data = file.read()
    
    decrypted_data = cipher_suite.decrypt(encrypted_data)

    # Write the decrypted data to a temporary .tar.gz file
    with open(decrypted_tar_path, "wb") as file:
        file.write(decrypted_data)
    
    # Extract the .tar.gz file into the output directory
    with tarfile.open(decrypted_tar_path, "r:gz") as tar:
        tar.extractall(path=output_directory)
    
    # Cleanup: Remove the temporary .tar.gz file
    import os
    os.remove(decrypted_tar_path)

    # print(f"Folder extracted to: {output_directory}")
    return output_directory



#### DeepSAP v0.1.0 

import os
import subprocess
from datetime import datetime
import shutil
import pysam
from pyfaidx import Fasta


class SpliceJunction:
    """
    A class to store detailed information about a unique splice junction,
    including its annotation status and read support.
    """
    def __init__(self, chrom, start, end, strand, motif, ref_donor, ref_acceptor, biotype="novel"):
        self.chrom = chrom
        self.start = start          # 0-based coordinate of the last base of the donor exon
        self.end = end              # 0-based coordinate of the first base of the acceptor exon
        self.strand = strand
        self.motif = motif
        self.biotype = biotype      # "novel" or the transcript_type from GTF
        self.ref_donor = ref_donor
        self.ref_acceptor = ref_acceptor
        self.read_support = {}      # {read_id: (donor_diffs, acceptor_diffs)}

    def add_read_support(self, read_id, read_donor_diffs, read_acceptor_diffs):
        """Adds evidence from a single read alignment."""
        self.read_support[read_id] = (read_donor_diffs, read_acceptor_diffs)

def parse_md_tag(md_string):
    """Parses an MD tag string into a list of operations."""
    return re.findall(r'(\d+)|(\^[A-Z]+)|([A-Z])', md_string)

def match_motif_with_mismatch(sequence, motif):
    """Checks if a 4-base sequence matches a 4-base motif with at most one mismatch."""
    if len(sequence) != 4 or len(motif) != 4:
        return False
    mismatches = sum(1 for i in range(4) if sequence[i] != motif[i])
    return mismatches <= 1

def run_gsnap_and_filter_stream(gsnap_flags, threads, gsnap_index, gsnap_index_opts, 
                                mate_1, mate_2, sm_view_flags, sm_sort_flags, 
                                output_prefix, log_file, genome_fasta, annotated_junctions_table, window_size):
    
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tStarting single-pass GSNAP streaming workflow.")

    clean_bam_path = f"{output_prefix}.bam"
    junctions_tsv_path = f"{output_prefix}_junctions.tsv"

    gsnap_cmd = f"gsnap {gsnap_flags} -t {threads} -D {gsnap_index} {gsnap_index_opts} {mate_1} {mate_2}"
    clean_samtools_cmd = f"samtools view -u -b -h {sm_view_flags} - | samtools sort {sm_sort_flags} -o {clean_bam_path} -T {output_prefix}_clean_temp -@ {threads}"
    
    junction_table = {}

    with open(log_file, "a") as log:
        try:
            gsnap_process = subprocess.Popen(gsnap_cmd, shell=True, stdout=subprocess.PIPE, stderr=log, text=True)
            clean_process = subprocess.Popen(clean_samtools_cmd, shell=True, stdin=subprocess.PIPE, stderr=log, text=True)

            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tProcesses started. Filtering stream...")

            current_read_id, read_buffer, has_n_cigar, header_lines = None, [], False, []

            for line in gsnap_process.stdout:
                if line.startswith('@'):
                    header_lines.append(line)
                    clean_process.stdin.write(line)
                    continue

                fields = line.strip().split('\t')
                read_id = fields[0]

                if read_id != current_read_id and current_read_id is not None:
                    if has_n_cigar:
                        update_junction_table(junction_table, read_buffer, "".join(header_lines), genome_fasta, annotated_junctions_table, window_size)
                    else:
                        for record in read_buffer:
                            clean_process.stdin.write(record)
                    
                    read_buffer, has_n_cigar = [], False
                
                current_read_id = read_id
                read_buffer.append(line)
                if 'N' in fields[5]:
                    has_n_cigar = True

            if current_read_id is not None:
                if has_n_cigar:
                    update_junction_table(junction_table, read_buffer, "".join(header_lines), genome_fasta, annotated_junctions_table, window_size)
                else:
                    for record in read_buffer:
                        clean_process.stdin.write(record)

            clean_process.stdin.close()
            gsnap_retcode = gsnap_process.wait()
            clean_retcode = clean_process.wait()

            if gsnap_retcode != 0 or clean_retcode != 0:
                raise subprocess.CalledProcessError(gsnap_retcode or clean_retcode, gsnap_cmd)

            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tIndexing final clean BAM file...")
            pysam.index(clean_bam_path)
            
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tWriting {len(junction_table)} collected junctions to {junctions_tsv_path}...")
            with open(junctions_tsv_path, 'w') as f_out:
                f_out.write("junction_id\tbiotype\tmotif\tref_donor\tref_acceptor\tsupporting_reads\n")
                for junc_id, junc_obj in junction_table.items():
                    reads_str = ",".join(junc_obj.read_support.keys())
                    junc_id_1based = f"{junc_obj.chrom}__{junc_obj.strand}__{junc_obj.start + 1}__{junc_obj.end + 1}"
                    f_out.write(f"{junc_id_1based}\t{junc_obj.biotype}\t{junc_obj.motif}\t{junc_obj.ref_donor}\t{junc_obj.ref_acceptor}\t{reads_str}\n")
            
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tCalculating motif statistics...")
            forward_signals, reverse_signals = {}, {}
            for junc in junction_table.values():
                canonical_motif = junc.motif if junc.strand == '+' else get_reverse_complement(junc.motif)
                if junc.strand == '+':
                    forward_signals[canonical_motif] = forward_signals.get(canonical_motif, 0) + 1
                else:
                    reverse_signals[canonical_motif] = reverse_signals.get(canonical_motif, 0) + 1
            
            all_signals = sorted(list(set(list(forward_signals.keys()) + list(reverse_signals.keys()))))
            df = pd.DataFrame(columns=['Signal', 'Forward', 'Reverse', 'Percentage'])
            total = sum(forward_signals.values()) + sum(reverse_signals.values())
            
            if total > 0:
                for i, signal in enumerate(all_signals):
                    fwd_count = forward_signals.get(signal, 0)
                    rev_count = reverse_signals.get(signal, 0)
                    percentage = (fwd_count + rev_count) * 100 / total
                    df.loc[i] = [signal, fwd_count, rev_count, round(percentage, 2)]
                
                print("\n--- Top 10 Splicing Signals Types ---")
                print(df.sort_values('Percentage', ascending=[False]).head(10).to_string(index=False))
                print("-------------------------------------\n")

            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[SUCCESS]\tWorkflow complete.")

        except Exception as e:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[ERROR]\tAn error occurred during the streaming workflow: {e}")
            exit(1)

def update_junction_table(junction_table, sam_records, sam_header, genome, annotated_junctions_table, window_size):
    """
    Parses SAM records, extracts junction info, determines strand from XS tag or motif,
    and annotates as 'novel' or with a GTF biotype.
    """
    header = pysam.AlignmentHeader.from_text(sam_header)
    w = window_size // 2

    for sam_string in sam_records:
        try:
            record = pysam.AlignedSegment.fromstring(sam_string, header)
            if record.is_unmapped or 'N' not in record.cigarstring:
                continue

            ref_pos = record.reference_start # 0-based
            for op, length in record.cigartuples:
                if op == 3: # Splice (N)
                    donor_coord_0based = ref_pos - 1
                    acceptor_coord_0based = ref_pos + length
                    
                    chrom = record.reference_name
                    
                    # --- NEW STRAND AND MOTIF LOGIC ---
                    raw_motif = str(genome[chrom][donor_coord_0based + 1 : donor_coord_0based + 3]).upper() + \
                                str(genome[chrom][acceptor_coord_0based - 2 : acceptor_coord_0based]).upper()

                    try:
                        strand = record.get_tag('XS')
                    except KeyError:
                        if match_motif_with_mismatch(raw_motif, "GTAG"):
                            strand = '+'
                        elif match_motif_with_mismatch(raw_motif, "CTAC"):
                            strand = '-'
                        else:
                            strand = '-' if record.is_reverse else '+'
                    
                    final_motif = raw_motif if strand == '+' else get_reverse_complement(raw_motif)
                    junction_id_1based = f"{chrom}__{strand}__{donor_coord_0based + 1}__{acceptor_coord_0based + 1}"

                    if junction_id_1based not in junction_table:
                        donor_seq_start = donor_coord_0based - w + 1
                        donor_seq_end = donor_coord_0based + w + 1
                        acceptor_seq_start = acceptor_coord_0based - w + 1
                        acceptor_seq_end = acceptor_coord_0based + w + 1

                        ref_donor_seq = str(genome[chrom][donor_seq_start:donor_seq_end]).upper()
                        ref_acceptor_seq = str(genome[chrom][acceptor_seq_start:acceptor_seq_end]).upper()
                        
                        biotype = annotated_junctions_table.get(junction_id_1based, {}).get('type', 'novel')
                        
                        junction_table[junction_id_1based] = SpliceJunction(
                            chrom, donor_coord_0based, acceptor_coord_0based, strand, final_motif, 
                            ref_donor_seq, ref_acceptor_seq, biotype
                        )
                    
                    junction_table[junction_id_1based].add_read_support(record.query_name, {}, {})

                if op in [0, 2, 3]: # M, D, N advance the reference pointer
                    ref_pos += length
        except (ValueError, KeyError, IndexError):
            continue

### old ###

def run_gsnap_index(anno_gtf, genome_fasta, gsnap_index, gsnap_idx_flags, log_file):
    """
    Creates GSNAP index files using specified parameters.
    
    Parameters:
    - anno_gtf (str): Path to the annotation GTF file.
    - genome_fasta (str): Path to the genome FASTA file.
    - gsnap_index (str): Directory where the GSNAP index files will be created.
    - log_file (str): Path to the log file for command output and error messages.
    """
    
    cmd_genes = f"cat {anno_gtf} | gtf_genes > {gsnap_index}/chr.genes.txt"
    cmd_gmap_build_genome = f"gmap_build -D {gsnap_index} -d index {genome_fasta}"
    cmd_gmap_build_transcriptome = f"gmap_build -D {gsnap_index} {gsnap_idx_flags} --genes={gsnap_index}/chr.genes.txt"

     # Open the log file in append mode
    with open(log_file, "a") as log:
        try:
            subprocess.run(cmd_genes, shell=True, stdout=log, stderr=log, check=True)
            subprocess.run(cmd_gmap_build_genome, shell=True, stdout=log, stderr=log, check=True)
            subprocess.run(cmd_gmap_build_transcriptome, shell=True, stdout=log, stderr=log, check=True)
            
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during GSNAP indexing: {e}")
            exit()

def run_gsnap_and_filter_stream_old(gsnap_aln_flags, gsnap_index, gsnap_idx_flags, threads, 
                                mate_1, mate_2, sam_view_flags, sam_sort_flags, 
                                output_prefix, log_file):
    """
    Runs a full GSNAP alignment and filtering workflow in a single, memory-efficient stream.

    This function pipes GSNAP output to a Python filter, which then routes reads
    to one of two parallel samtools pipelines. The "clean" reads are sorted and saved
    as a BAM file, while the "filtered" reads are saved directly as an unsorted BAM file.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tStarting single-pass GSNAP streaming workflow.")

    # Define final output paths
    clean_bam_path = f"{output_prefix}.bam"
    filtered_bam_path = f"{output_prefix}_filtered.bam"

    # Construct the pipeline commands
    gsnap_cmd = f"gsnap {gsnap_aln_flags} -t {threads} -D {gsnap_index} {gsnap_idx_flags} {mate_1} {mate_2}"
    
    # Command for the "clean" BAM file (sorted)
    clean_samtools_cmd = f"samtools view -u -b -h {sam_view_flags} - | samtools sort {sam_sort_flags} -o {clean_bam_path} -T {output_prefix}_clean_temp -@ {threads//2 or 1}"
    
    # MODIFIED: Command for the "filtered" BAM file (unsorted)
    filtered_samtools_cmd = f"samtools view -b -h {sam_view_flags} -o {filtered_bam_path} -"

    with open(log_file, "a") as log:
        try:
            # Start the GSNAP process
            gsnap_process = subprocess.Popen(gsnap_cmd, shell=True, stdout=subprocess.PIPE, stderr=log, text=True)
            
            # Start the two parallel samtools processes
            clean_process = subprocess.Popen(clean_samtools_cmd, shell=True, stdin=subprocess.PIPE, stderr=log, text=True)
            filtered_process = subprocess.Popen(filtered_samtools_cmd, shell=True, stdin=subprocess.PIPE, stderr=log, text=True)

            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tProcesses started. Filtering stream...")

            current_read_id = None
            read_buffer = []
            has_n_cigar = False

            # Process the stream from GSNAP line by line
            for line in gsnap_process.stdout:
                if line.startswith('@'):
                    # Write header lines to both downstream processes
                    clean_process.stdin.write(line)
                    filtered_process.stdin.write(line)
                    continue

                fields = line.strip().split('\t')
                read_id = fields[0]

                if read_id != current_read_id:
                    # We've encountered a new read, so process the previous one
                    if current_read_id is not None:
                        if has_n_cigar:
                            for record in read_buffer:
                                filtered_process.stdin.write(record)
                        else:
                            for record in read_buffer:
                                clean_process.stdin.write(record)
                    
                    # Reset for the new read
                    current_read_id = read_id
                    read_buffer = []
                    has_n_cigar = False

                # Buffer the current line and check its CIGAR
                read_buffer.append(line)
                cigar = fields[5]
                if 'N' in cigar:
                    has_n_cigar = True

            # Process the very last read in the file
            if current_read_id is not None:
                if has_n_cigar:
                    for record in read_buffer:
                        filtered_process.stdin.write(record)
                else:
                    for record in read_buffer:
                        clean_process.stdin.write(record)

            # Close the stdin pipes to signal that we're done sending data
            clean_process.stdin.close()
            filtered_process.stdin.close()

            # Wait for all processes to finish
            gsnap_retcode = gsnap_process.wait()
            clean_retcode = clean_process.wait()
            filtered_retcode = filtered_process.wait()

            if gsnap_retcode != 0 or clean_retcode != 0 or filtered_retcode != 0:
                raise subprocess.CalledProcessError(
                    gsnap_retcode or clean_retcode or filtered_retcode, gsnap_cmd
                )

            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tIndexing final clean BAM file...")
            pysam.index(clean_bam_path)
            
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tClean alignments saved to: {clean_bam_path}")
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tFiltered alignments saved to: {filtered_bam_path}")

        except Exception as e:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[ERROR]\tAn error occurred during the streaming workflow: {e}")
            exit(1)



####
if __name__ == '__main__':
    version = "v0.1.0"
    command = "DeepSAP "
    command += " ".join([arg for arg in sys.argv[1:]])

    parser = argparse.ArgumentParser(description = 'DeepSAP')

    parser.add_argument('-o','--out',       help='Path to output folder',                                      required=True, type=str)
    parser.add_argument('--prefix',         help="Prefix to ouput files",                                   dest='prefix', required=True, type=str, default ="")
    parser.add_argument('-s','--sam',       help='Path to the SAM/BAM file or director of files',              required=False, type=str, default ="")
    parser.add_argument('-f','--fasta',     help='Path to FASTA genome file compatible with the BAM file',     required=True, type=str)
    parser.add_argument('-g','--gtf',       help='Path to GTF annotation file compatible with the BAM file',   required=False, type=str, default ="")
    parser.add_argument('-c','--config',    help='Config .json file',                                          required=True, type=str)
    # parser.add_argument('--samtools',       help='Path to samtools',                                           required=False, type=str, default="samtools")
    parser.add_argument('--model_name',     help='Tranformer pre-trained model name',                       required=True, type=str)
    parser.add_argument('--model_path',     help='Tranformer model path',                                   required=True, type=str)
    parser.add_argument('--score_method',   help='Splice-junction scoring method',                          required=False, default="Add")
    parser.add_argument('--batch',          help='Batch size when performing prediction ',                    required=False, type=int, default=2048)
    parser.add_argument('--set_size',       help='Datapoints set is splitted into sets of this size to avoid performing prediction  on the whole datasets at once', required=False, type=int, default= 1024 * 100)
    parser.add_argument('--max_seq',        help='Maximum sequence length',                                 required=True, type=int)
    parser.add_argument('--seed',           help='Randomness seed',                                         required=False, type=int, default=42)
    parser.add_argument('-w','--window',    help='Window size around splice-junction used by the Transformer Model', required=True, type=int)
    parser.add_argument('-k','--kmer',      help='Kmer size used by the Transformer Model',                          required=True, type=int)
    parser.add_argument('-t','--threads',   help='Number of threads',                          required=False, type=int, default=os.cpu_count())
    parser.add_argument('--fp16',           help='Use fp16', action='store_true')
    parser.add_argument('--no-fp16',        help="Don't use fp16", dest='fp16', action='store_false')
    parser.add_argument('--score_reads',    help='Classify reads using Transformer model and ddd reads scores to SAM', action='store_true')
    parser.add_argument('--n_reads',        help='Number of reads to be classified using Transformer model if --score_reads is used ', required=False, type=int, default=100000)
    parser.add_argument('--tokens',         help='Get first 1k tokens', action='store_true')
    parser.add_argument('--test',           help="Don't perform prediction", dest='test', action='store_true')
    parser.add_argument('--mate_1',         help="Fastq file of mate 1 (RNA-seq short reads)", dest='mate_1',    required=False, type=str, default ="")
    parser.add_argument('--mate_2',         help="Fastq file of mate 2 (RNA-seq short reads)", dest='mate_2',    required=False, type=str, default ="")
    
    parser.add_argument('--gsnap_idx',      help="GSNAP index path ",                      dest='gsnap_idx',        required=False, type=str, default ="")
    parser.add_argument('--gsnap_idx_flags',help="GSNAP indexing flags ",                  dest='gsnap_idx_flags',  required=False, type=str, default ="-d index -c transcriptome")
    parser.add_argument('--gsnap_aln_flags',help="GSNAP alignment flags ",                 dest='gsnap_aln_flags',  required=False, type=str, default ="--gunzip -A sam --novelsplicing 1")

    parser.add_argument('--samtools',       help='Samtools binary path',                   dest='samtools',         required=False, type=str, default="samtools")
    parser.add_argument('--sam_sort_flags', help='Samtools sort flags',                    dest='sam_sort_flags',   required=False, type=str, default="-m 500M")
    parser.add_argument('--chr_coverage',   help='Max bases per chromosome to scan',                default=10**12, type=int)
    parser.add_argument('--motifs',         help='4-mers, space/comma-separated (e.g. "GTAG ATAC GCAG")', default="GTAG ATAC GCAG", type=str)
    parser.add_argument('--min_prob',       required=False, type=float, default=0.10, help='Minimum probability to report a site')

    parser.set_defaults(fp16=True)
    
    print("\n[{}]\t[INFO]\tRunning DeepSAP {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), version))

    args = parser.parse_args()

    ## load model 
    clean_up = False
    if str(args.model_path).endswith(".tar.gz"):
        model_tag = os.path.basename(args.model_path).removesuffix(".tar.gz")
        clean_up = True 

        if model_tag == 'ADNRBE_T11MS50':
            encryption_key = 'tV5xffLxBbOJxZI5laxLBicMlvGd8svmYehI9Hg0cOs='
            model_path_tmp = "/tmp/76/12/09/bin/dir/contains/include/"
            os.makedirs(model_path_tmp, exist_ok=True)

            load_folder_from_encrypted_file(args.model_path, encryption_key, model_path_tmp)
            args.model_path = model_path_tmp + "/DNABERT1_MS150"

    # Setting up output folder
    if not args.out.endswith("/"):
        args.out += "/"


    os.makedirs(args.out, exist_ok=True)
    
    if args.fasta != "" and args.gtf == "" and args.sam == "" and args.mate_1 == "" and args.mate_2 == "":
        """
        Motif	Spliceosome Type	            Frequency in Humans	        Reason to Include
        GT-AG	Major (U2-type)	                Canonical (>98%)	        This is the most common splicing signal.
        GC-AG	Major (U2-type)	                Common non-canonical (~1%)	The most frequent alternative to GT-AG and is still processed by the main splicing machinery.
        AT-AC	Both Major (U2) & Minor (U12)	Rare (~0.1%)	            Capture rare events. This motif is recognized by both splicing systems and is crucial for identifying a subset of rare introns.
        AT-AA	Minor (U12-type) only	        Very rare (~0.05%)	        Capture a different pathway. The primary signal for the distinct U12 splicing pathway, which often regulates fundamentally important genes.

        D+: i, D-: i+1, A+: i+1, A-: i

        | Biological site (on gene's strand)      | What we scan on **forward reference** (2-mer at `seq[i:i+2]`)  | Class | Locus to report (`pos`)                                                   |
        | --------------------------------------- | -------------------------------------------------------------- | ----- | ------------------------------------------------------------------------- |
        | **Donor, + strand (5′ splice site)**    | `GT`, `GC`, `AT`                                               | `D`   | `pos = i` (the **first** base of the 2-mer; e.g., the **G** in `GT`)      |
        | **Donor, - strand (5′ splice site)**    | `AC` (rc of `GT`), `GC` (pal), `AT` (pal)                      | `D`   | `pos = i + 1` (the **second** base of the 2-mer; e.g., the **C** in `AC`) |
        | **Acceptor, + strand (3′ splice site)** | `AG`, `AC`                                                     | `A`   | `pos = i + 1` (the **second** base of the 2-mer; e.g., the **G** in `AG`) |
        | **Acceptor, - strand (3′ splice site)** | `CT` (rc of `AG`), `GT` (rc of `AC`)                           | `A`   | `pos = i` (the **first** base of the 2-mer; e.g., the **C** in `CT`)      |
        """

        # Motifs
        motif_pairs = parse_motifs_arg(args.motifs) or ["GTAG","ATAC","GCAG"]
        donors_fwd, acceptors_fwd = parse_motif_pairs(motif_pairs)

        # Flanks: center window on boundary (acceptor left is 2 less to keep total window)
        donor_left  = int(args.window // 2)
        accept_left = int(args.window // 2) - 2

        # Build predict_args (clean, non-duplicated)
        predict_args = {
            'model_name': "/DNA_bert_6/",
            'model_path': args.model_path,
            'kmer': args.kmer,
            'batch': args.batch,
            'seed': args.seed,
            'trust_remote_code': False,
            'local_files_only': True,
            'num_labels': 3,

            'donor_left': donor_left,
            'accept_left': accept_left,
            'donor_label_idx': 0,
            'acceptor_label_idx': 1,
            'donors_fwd': donors_fwd,
            'acceptors_fwd': acceptors_fwd,

            # precision & speed (target stack)
            'dtype': 'bf16',            # weights in bfloat16
            'amp': True,                # autocast
            'allow_tf32': True,         # TF32 matmul
            'compile': False,           # keep False unless Triton toolchain is ready
            'bettertransformer': True,  # fused fastpath
            'low_cpu_mem_usage': True,
            'device_map': 'cuda',

            # tokenizer behavior (fixed windows -> no padding)
            'padding': 'none',
            'add_special_tokens': True,
            'max_seq': 0,               # ignored with padding='none'

            # scan coverage
            'chr_coverage': int(args.chr_coverage),

            # reporting filter
            'min_prob': float(args.min_prob),
        }

        # Only do FASTA→CSV path in this entrypoint
        if args.fasta and not any([args.gtf, args.sam, args.mate_1, args.mate_2]):
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tStarting splice-site prediction")

            t0 = time.time()
            predict_splice_junctions_from_FASTA(
                fasta_path=args.fasta,
                window_size=args.window,
                n_datapoints_per_run=args.set_size,
                chr_coverage=args.chr_coverage,
                output_path=args.out,
                predict_args=predict_args
            )
            t1 = time.time()
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tDone in {(t1 - t0)/60:.2f} min")
        exit()

    elif args.gsnap_idx == "" and args.fasta != "" and args.gtf != "" and args.sam == "" and args.mate_1 == "" and args.mate_2 == "":
        
        start_time = time.time()
        args.gsnap_idx = os.path.join(args.out , args.prefix)
        os.makedirs(args.gsnap_idx, exist_ok=True)

        gsnap_indexing_log = os.path.join(args.gsnap_idx , "gsnap_indexing_log.txt")

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tBuilding GSNAP TGGA index '{args.gsnap_idx}'")

        run_gsnap_index(
            anno_gtf=args.gtf,
            genome_fasta=args.fasta,
            gsnap_index=args.gsnap_idx,
            gsnap_idx_flags=args.gsnap_idx_flags,
            log_file=gsnap_indexing_log)
        
        end_time = time.time()
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tTotal execution time: {(end_time - start_time)/60:.2f} minutes.")

        exit()
    
    elif args.gsnap_idx != "" and args.fasta != "" and args.gtf != "" and args.mate_1 != "" and args.mate_2 != "" and args.sam == "":
        start_time = time.time()
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tProcessing paired-end short RNA-seq reads ")

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tParsing GTF file '{args.gtf}'")

        genome_fasta = Fasta(args.fasta)

        genes_table, transcripts_table, exons_table = parse_gtf(args.gtf)
        
        annotated_junctions_table = collect_junctions_sequences_from_GTF(
            genome_fasta,
            transcripts_table, 
            exons_table, 
            args.window, 
            "", # transcript_types
            is_verbose=False,
            is_training=True,
            is_strict=True,
            add_junction_seq=True
        )
        
        # print("\n--- Example Junctions ---")
        # count = 0
        # for junc_id, junc_obj in annotated_junctions_table.items():
        #     print(f"ID: {junc_id}, Motif: {junc_obj.signal}, Donor Seq Length: {len(junc_obj.donor_kmers)}")
        #     count += 1
        #     if count >= 100:
        #         break

        gsnap_aligning_log = os.path.join(args.out, args.prefix + "_gsnap_aligning_log.txt")
        final_output_prefix = os.path.join(args.out, args.prefix)
        
        

        os.makedirs(args.out, exist_ok=True)
        
        # Call the new single-pass streaming function
        run_gsnap_and_filter_stream(
            gsnap_flags=args.gsnap_aln_flags,
            threads=args.threads, 
            gsnap_index=args.gsnap_idx,
            gsnap_index_opts=args.gsnap_idx_flags,
            mate_1=args.mate_1, 
            mate_2=args.mate_2, 
            sm_view_flags="", 
            sm_sort_flags=args.sam_sort_flags, 
            output_prefix=final_output_prefix,
            log_file=gsnap_aligning_log,
            genome_fasta=genome_fasta,
            annotated_junctions_table=annotated_junctions_table,
            window_size=args.window
        )
        
        end_time = time.time()
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[INFO]\tTotal execution time: {(end_time - start_time)/60:.2f} minutes.")
        exit()



        # gsnap_aligning_log = os.path.join(args.out , args.prefix + "_gsnap_aligning_log.txt")
        # args.sam = os.path.join(args.out , args.prefix + "_gsnap.bam")

        # gsnap_tmp = os.path.join(args.out , args.prefix + "_gsnap_tmp/")
        # os.makedirs(gsnap_tmp, exist_ok=True)

        # print("[{}]\t[LOG]\tRunning GSNAP TGGA '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.sam))

        # run_gsnap_alignment(gsnap_aln_flags="--gunzip -A sam --novelsplicing 1", 
        #                     threads=args.threads, 
        #                     gsnap_index=args.gsnap_idx,
        #                     gsnap_idx_flags=args.gsnap_idx_flags,
        #                     mate_1 = args.mate_1, 
        #                     mate_2 = args.mate_2, 
        #                     sm_view_flags = "", 
        #                     sam_sort_flags = "", 
        #                     output_bam = args.sam, 
        #                     output_prefix = gsnap_tmp,
        #                     log_file = gsnap_aligning_log)
        

        # // another function that takes the ouput from run_gsnap_alignment and do two things:
        # 1- Filter out read Id if any alignment record has N cigar and write to {args.prefix}_filtered.bam  
        # 2- otherwise, write to bam file. {args.prefix}.bam
        
        exit()
        

    if args.mate_1 != "" and args.mate_2 != "":
        print("[{}]\t[LOG]\tRunning GSNAP".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

        if args.gsnap_idx == "":
            args.gsnap_idx = os.path.join(args.out , "gsnap_idx")
            os.makedirs(args.gsnap_idx, exist_ok=True)

            gsnap_indexing_log = os.path.join(args.out , "gsnap_indexing_log.txt")
            print("[{}]\t[LOG]\tBuilding GSNAP TGGA index '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.gsnap_idx))
            
            run_gsnap_index(
                anno_gtf=args.gtf,
                genome_fasta=args.fasta,
                gsnap_index=args.gsnap_idx,
                log_file=gsnap_indexing_log
            )

        gsnap_aligning_log = os.path.join(args.out , args.prefix + "_gsnap_aligning_log.txt")
        args.sam = os.path.join(args.out , args.prefix + "_gsnap.bam")

        gsnap_tmp = os.path.join(args.out , args.prefix + "_gsnap_tmp/")
        os.makedirs(gsnap_tmp, exist_ok=True)

        print("[{}]\t[LOG]\tRunning GSNAP TGGA '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.sam))

        run_gsnap_alignment(gsnap_aln_flags="--gunzip -A sam --novelsplicing 1", 
                            threads=args.threads, 
                            gsnap_index=args.gsnap_idx, 
                            mate_1 = args.mate_1, 
                            mate_2 = args.mate_2, 
                            sm_view_flags = "", 
                            sam_sort_flags = "", 
                            output_bam = args.sam, 
                            output_prefix = gsnap_tmp,
                            log_file = gsnap_aligning_log)
    
    # Collecting input
    if os.path.isfile(args.sam):
        sams_files = [args.sam]
    else:
        dir_content = os.listdir(args.sam)
        sams_files = ["{}/{}".format(args.sam, i) for i in dir_content if i.endswith('.sam') or i.endswith('.bam')]
    
    # Collecting FASTA sequences
    # print("[{}]\t[LOG]\tParsing the FASTA file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.fasta))
    genome_sequences_table = parse_fasta(args.fasta)

    # Collect Gene annotation
    # print("[{}]\t[LOG]\tParsing the GTF file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.gtf))
    genes_table, transcripts_table, exons_table = parse_gtf(args.gtf)
    
    # Collect Transcript information
    # print("[{}]\t[LOG]\tTranscript information".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    get_transcripts_info(transcripts_table, exons_table)
    
    # Generate kmers from splice junctions
    print("[{}]\t[LOG]\tCollecting splice junctions from GTF".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    annotated_junctions_table = collect_junctions_sequences_from_GTF(genome_sequences_table,
                                                                     transcripts_table, 
                                                                     exons_table, 
                                                                     args.window, 
                                                                     "",
                                                                     is_verbose  = False,
                                                                     is_training = False,
                                                                     is_strict   = False)
    
    found_junctions_table = dict()
    # Submit tasks to the executor and gather the results
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        tasks = [executor.submit(collect_junctions_from_SAM,
                                 annotated_junctions_table, genome_sequences_table, 
                                 sam_file, args.window, args.samtools) 
                                 for sam_file in sams_files]
        
        # for future in concurrent.futures.as_completed(tasks):
        #         found_junctions_table.update(future.result())
        #         print("[{}]\t[LOG]\tFinished parsing a SAM file, len(found_junctions_table)= {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), len(found_junctions_table)))
        
        for future in concurrent.futures.as_completed(tasks):
            try:
                result = future.result()  # If an exception occurred in the worker, it's raised here
                found_junctions_table.update(result)
                print("[{}]\t[LOG]\tFinished parsing a SAM file, len(found_junctions_table)= {}".format(
                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    len(found_junctions_table)
                ))
            except Exception as e:
                # Print or log the stack trace
                import traceback
                traceback.print_exc()
                print(f"[ERROR] Future raised an exception: {e}")
                # Optionally re-raise or sys.exit(1) if you need to bail out
                raise

    # Deleting un-necrssory objects
    del annotated_junctions_table, genome_sequences_table
    gc.collect()
    
    n_junctions_per_call = args.set_size
    iterator   = DictBatchIterator(found_junctions_table, n_junctions_per_call)

    if args.test:
        for batch_num, junctions_batch in enumerate(iterator, start=1):
            dataset_path = os.path.join(args.out , f"{args.prefix}_prediction_batch_{batch_num}/")

            if os.path.isdir(dataset_path):
                shutil.rmtree(dataset_path)
                        
            os.makedirs(dataset_path)

            zeros_array  = np.zeros((2 * len(junctions_batch.keys()), 3))
            np.save(dataset_path + "probs.npy", zeros_array)

    else:
        # Iterate through junctions batches to perform prediction (Batching to ovoid high memory usage): 
        for batch_num, junctions_batch in enumerate(iterator, start=1):
            dataset_path = os.path.join(args.out , f"{args.prefix}_prediction_batch_{batch_num}/")

            if os.path.isdir(dataset_path):
                shutil.rmtree(dataset_path)

            os.makedirs(dataset_path)

            print("[{}]\t[LOG]\tGenerating splice-junction prediction dataset batch: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), batch_num))
            write_prediction_datapoints(junctions_batch, dataset_path)

            print("[{}]\t[LOG]\tPredicting found splice junctions using DNABERT MS150".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            predict(args.model_name, args.model_path, dataset_path , args.kmer, args.max_seq, dataset_path, args.batch, args.seed, is_fp16=bool(args.fp16))

    # Generating prediction datapoints
    # print("[{}]\t[LOG]\tGenerating splice-junction prediction dataset".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # dataset_path = args.out + "/prediction/"
    # write_prediction_datapoints(found_junctions_table, dataset_path)

    # # Prediction using transfomer model
    # print("[{}]\t[LOG]\tPredicting found splice junctions using {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args.model_name))
    # predict(args.model_name, args.model_path, dataset_path , 
    #         args.kmer, args.max_seq, dataset_path, args.batch, args.seed, initate_process=bool(args.process), is_fp16=bool(args.fp16))
    
    # Generate Regions
    print("[{}]\t[LOG]\tGenerating genome regions ".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    genome_sequences_table = parse_fasta(args.fasta)
    regions_loci_table , regions_array_table, regions_relations_table, regions_transcripts_table, jumps_transcripts_table = generate_regions(genome_sequences_table, transcripts_table, exons_table, args.threads)
    
    # Deleting un-necrssory objects
    del genome_sequences_table, genes_table, transcripts_table, exons_table
    gc.collect()

    # Generate junctions score
    # junctions_scores = get_junctions_scores(dataset_path, found_junctions_table, method="Add")

    # Generate junctions score
    iterator   = DictBatchIterator(found_junctions_table, n_junctions_per_call)
    junctions_scores = dict()
    
    for batch_num, junctions_batch in enumerate(iterator, start=1):
        dataset_path = args.out + f"/{args.prefix}_prediction_batch_{batch_num}/"
        junctions_scores.update(get_junctions_scores(dataset_path, junctions_batch, method="Add"))

    # Write Junctions ifo 
    write_junctions_info(found_junctions_table, junctions_scores, args.out, args.prefix)

    # print("[{}]\t[LOG]\tCreating new sam file(s)".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:

        tasks = [executor.submit(filter_and_add_junctions_info_to_SAM,
                                 found_junctions_table, 
                                junctions_scores,
                                args.config,
                                sam_file, 
                                args.out,
                                args.prefix,
                                args.samtools,
                                regions_loci_table , 
                                regions_array_table, 
                                regions_relations_table, 
                                regions_transcripts_table, 
                                jumps_transcripts_table,
                                version, command) for sam_file in sams_files]
        
        for future in concurrent.futures.as_completed(tasks):
                print(future.result())

    # Deleting un-necrssory objects
    del found_junctions_table
    del regions_loci_table , regions_array_table, regions_relations_table, regions_transcripts_table, jumps_transcripts_table
    gc.collect()

    if (bool(args.score_reads)):
        for sam_file in sams_files:
            split_reads = True
            dataset_path = args.out + "/prediction_all_reads/"

            n_reads = write_prediction_datapoints_from_reads(sam_file, args.out, args.prefix, args.samtools, dataset_path, n_reads=args.n_reads)
        
            if split_reads == True:
                chunk_size = 65536
                data = pd.read_csv(args.out + "/prediction_all_reads/dev.csv")
                npy_files = [] 

                for i, chunk_start in enumerate(range(0, n_reads, chunk_size)):
                    
                    chunk_end = min(chunk_start + chunk_size, n_reads)
                    chunk_data = data.iloc[chunk_start:chunk_end].copy()

                    rows = [['NNNNNNNNNNNNNNN',0], ['NNNNNNNNNNNNNNN',1], ['NNNNNNNNNNNNNNN',2]]
                    needed_df = pd.DataFrame(rows, columns=['sequence', 'label'])

                    final_df = pd.concat([chunk_data, needed_df], axis=0)

                    dataset_path = args.out + f"/prediction_all_reads_{i+1}/"
                    npy_file = dataset_path + "/probs.npy"
                    dev_file = dataset_path + "/dev.csv"

                    if os.path.isdir(dataset_path):
                        shutil.rmtree(dataset_path)
                        
                    os.makedirs(dataset_path)

                    final_df.to_csv(dev_file, index=False)
                    npy_files.append(npy_file)

                    if args.test:
                        zeros_array  = np.zeros((chunk_end - chunk_start, 3))
                        np.save(npy_file, zeros_array)
                    else:
                        predict(args.model_name, args.model_path, dataset_path , args.kmer, args.max_seq, dataset_path, args.batch, args.seed, is_fp16=bool(args.fp16))

                probs_arrays = [np.load(npy_file)[:-3] for npy_file in npy_files]

                concatenated_probs_arrays = np.concatenate(probs_arrays)

                np.save(args.out + "/prediction_all_reads/probs", concatenated_probs_arrays)

            else:
                predict(args.model_name, args.model_path, dataset_path , args.kmer, args.max_seq, dataset_path, args.batch, args.seed, is_fp16=bool(args.fp16))
                
            dataset_path = args.out + "/prediction_all_reads/"
            add_reads_scores_to_sam(sam_file, args.out, args.prefix, args.samtools, dataset_path, n_reads=args.n_reads, method="Add")
    

    if clean_up:
        shutil.rmtree(args.model_path)
    print("[{}]\t[LOG]\tFinished successfuly".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
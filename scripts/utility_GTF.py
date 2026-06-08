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


import numpy as np
import gzip
import re
import pandas as pd
from pyfaidx import Fasta
from statistics import median, mean, mode
from datetime import datetime
import concurrent.futures
import threading
import os
from datetime import datetime

class Junction:
  def __init__(self, chrom, donor_ID, acceptor_ID, donor_end, acceptor_start, 
               strand, signal, transcript_type,
               donor_kmers, acceptor_kmers, junction_kmers):
    self.chrom    = chrom
    self.donor_ID = donor_ID
    self.acceptor_ID = acceptor_ID
    self.donor_end = donor_end
    self.acceptor_start = acceptor_start
    self.strand = strand
    self.signal = signal
    self.transcript_type = transcript_type
    self.donor_kmers = donor_kmers
    self.acceptor_kmers = acceptor_kmers
    self.junction_kmers = junction_kmers

    def __str__(self):
        representatin  = "Junction :\nDonor ID: {}\nAcceptor ID: {}\nSignal: {}\nJunction kmers: {}\n".format(self.donor_ID, self.acceptor_ID, self.signal, self.junction_kmers)
        return representatin

class JunctionFromRead:
  def __init__(self, signal, transcript_type,
               donor_seq, acceptor_seq, junction_seq, 
               read, score_js, score_rs):
    self.signal = signal
    self.transcript_type = transcript_type
    self.donor_seq = donor_seq
    self.acceptor_seq = acceptor_seq
    self.junction_seq = junction_seq
    self.read = read
    self.score_js = score_js
    self.score_rs = score_rs
    
class Gene:
    def __init__(self, chrom, start, end, strand, transcripts_set, exons_set):
        self.chrom           = chrom
        self.start           = start
        self.end             = end
        self.strand          = strand
        self.transcripts_set = transcripts_set
        self.exons_set       = exons_set
    
    def __str__(self):
        return "Gene:\nChrom: {}\nStart: {}\nEnd: {}\nStrand: {}\nTranscripts set: {}\nExons set: {}\n".format(self.chrom, self.start, self.end, self.strand, self.transcripts_set, self.exons_set)
         
class Transcript:
    def __init__(self, transcript_id, chrom, start, end, strand, gene_id, exons_list, transcript_type, junctions):
        self.transcript_id = transcript_id
        self.chrom         = chrom
        self.start         = start
        self.end           = end
        self.strand        = strand
        self.gene_id       = gene_id
        self.exons_list    = exons_list
        self.type          = transcript_type
        self.junctions     = junctions
    def __str__(self):
        return "Transcript: {}\nChrom: {}\nStart: {}\nEnd: {}\nStrand: {}\nGene ID {}\nExons list: {}\nType: {}\n".format(self.transcript_id, self.chrom, self.start, self.end, self.strand, self.gene_id, self.exons_list, self.type)
    
class Exon:
    def __init__(self, chrom, start, end, strand, exon_ID, exons_set, genes_set, trascripts_dict): 
        """Exon Object init function

        Args:
            chrom (string): Chromosome name
            start (int): start in chromosome
            end (int): end in chromosome
            strand (char): strand
            exons_set (set): set of exons ID that have the exact chr,start,end,strand values 
            genes_set (set): set of genes that this exon belong to, due to dupications an exon might belong to mulitple genes 
            trascripts_dict (dictionary): dictionary of {transcript_id : exon_number}
        """
        self.chrom      = chrom
        self.start      = start
        self.end        = end
        self.strand     = strand
        self.ID         = exon_ID
        self.exons_set  = exons_set
        self.genes_set  = genes_set
        self.trascripts_dict = trascripts_dict
        
    def __str__(self):
        return "Exon: {}\nChrom: {}\nStart: {}\nEnd: {}\nStrand: {}\nExons_set: {}\nGenes_set: {}\nTranscript dict: {}\n".format(self.ID, self.chrom, self.start, self.end, self.strand, self.exons_set, self.genes_set, self.trascripts_dict)
         
def parse_gtf(path):
    """_summary_
    Return hash tables for genes, transcripts and exons

    gene_id -> list(chrom, start, end, strand, list(trancripts_ids), list(exon_internal_ids))
    transcript_id -> list(chrom, start, end, strand, list(exon_internal_ids))
    exon_internal_id -> list(chrom, start, end, strand, list(exon_ids), list(gene_ids), dict({transcript_id : exon_number})]
    
    Args:
        path (string): Path to the GTF file
    
    """
    print("[{}]\t[LOG]\tParsing GTF file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), path))
    
    genes_table = dict()
    exons_table = dict()
    transcripts_table = dict()
    
    file = open(path, "r")
    zipped = False
    
    if path[-3:] == ".gz":
        file = gzip.open(path, 'rb')
        zipped = True
    
    pattern_gene        = r"([.]*)gene_id \"([^\"]*)\";([.]*)"
    pattern_transcript  = r"([.]*)transcript_id \"([^\"]*)\";([.]*)"
    pattern_type        = r"([.]*)transcript_type \"([^\"]*)\";([.]*)"
    pattern_biotypetype = r"([.]*)transcript_biotype \"([^\"]*)\";([.]*)"
    pattern_exon        = r"([.]*)exon_id \"([^\"]*)\";([.]*)"
    pattern_exon_number = r"([.]*)exon_number (\b[0-9]+);([.]*)"
    
    while True:
        line = file.readline()
        
        if zipped:
            line = line.decode()
            
        if len(line) == 0:
            break
        elif(line.startswith('#')):
            continue
        else:
            # print(line)
            fields = line.split('\t')
            feature = fields[2]
            chrom   = fields[0]
            start   = int(fields[3])
            end     = int(fields[4])
            strand  = fields[6]
            
            if feature == "gene":
                search_gene = re.search(pattern_gene, fields[8])

                try:
                    gene_id = search_gene.groups()[1]
                except:
                    print(line)
                    
                if gene_id == "":
                    print("gene_id is not found, line = '{}'".format(line))
                    exit()
                    
                elif (gene_id in genes_table.keys()):
                    if(genes_table[gene_id].start != 0):
                        print("gene_id is repeated {} in line:\n'{}'".format(gene_id, line))

                        # print(f"\ngenes_table[gene_id].start: {genes_table[gene_id].start}")
                        # print(f"genes_table[gene_id].end: {genes_table[gene_id].end}")

                        assert genes_table[gene_id].chrom  == chrom
                        # assert genes_table[gene_id].strand == strand

                        if genes_table[gene_id].start > start:
                            genes_table[gene_id].start = start

                        if genes_table[gene_id].end < end:
                            genes_table[gene_id].end = end
                        
                        # print(f"genes_table[gene_id].start: {genes_table[gene_id].start}")
                        # print(f"genes_table[gene_id].end: {genes_table[gene_id].end}")
                    else:
                        print("GTF is ill sorted, gene_id is after its exons, line = '{}'".format(line))
                            
                        genes_table[gene_id].chrom  = chrom
                        genes_table[gene_id].start  = start
                        genes_table[gene_id].end    = end
                        genes_table[gene_id].strand = strand
                    
                else:
                    genes_table[gene_id] = Gene(chrom, start, end, strand, set(), set())
                    
            
            elif feature == "transcript":
                search_gene       = re.search(pattern_gene, fields[8])
                search_transcript = re.search(pattern_transcript, fields[8])
                search_type       = re.search(pattern_type, fields[8])

                gene_id         = search_gene.groups()[1]
                transcript_id   = search_transcript.groups()[1]
                transcript_type = ""
                
                if (search_type != None):
                    transcript_type = search_type.groups()[1]
                else:
                    search_type  = re.search(pattern_biotypetype, fields[8])
                    if (search_type != None):
                        transcript_type = search_type.groups()[1]
                    else:
                        transcript_type = "Unknow"
                    
                if gene_id == "":
                    print("gene_id is not found, line = '{}'".format(line))
                    exit()
                    
                elif transcript_id == "":
                    print("transcript_id is not found, line = '{}'".format(line))
                    exit()  
                    
                elif (transcript_id in transcripts_table.keys()):
                    if(transcripts_table[transcript_id].start != 0):
                        print("gene_id is repeated, line = '{}'".format(line))
                        exit()
                        
                    else:
                        print("GTF is ill sorted, transcript is after its exons, line = '{}'".format(line))
                        
                    transcripts_table[transcript_id].chrom  = chrom
                    transcripts_table[transcript_id].start  = start
                    transcripts_table[transcript_id].end    = end
                    transcripts_table[transcript_id].strand = strand
                    transcripts_table[transcript_id].type = transcript_type
                else:
                    transcripts_table[transcript_id] = Transcript(transcript_id, chrom, start, end, strand, gene_id, list(), transcript_type, [])

            elif feature == "exon":
                
                search_gene       = re.search(pattern_gene, fields[8])
                search_transcript = re.search(pattern_transcript, fields[8])
                search_exon       = re.search(pattern_exon, fields[8])
                search_number     = re.search(pattern_exon_number, fields[8])
                
                gene_id       = search_gene.groups()[1]
                transcript_id = search_transcript.groups()[1]
                
                if (search_exon != None):
                    exon_id = search_exon.groups()[1]
                else:
                    exon_id = ""
                    
                if (search_number != None):
                    exon_number = int(search_number.groups()[1])
                else:
                    exon_number = 0
                
                exon_internal_id = chrom + "__" + strand + "__" + str(start) + "__" + str(end)
                
                if gene_id == "":
                    print("gene_id is not found, line = '{}'".format(line))
                    exit()

                if transcript_id == "":
                    print("transcript_id is not found, line = '{}'".format(line))
                    exit()
                    
                if gene_id in genes_table.keys():
                    genes_table[gene_id].transcripts_set.add(transcript_id)
                    genes_table[gene_id].exons_set.add(exon_internal_id)
                else:
                    genes_table[gene_id] = Gene(chrom, 0, 0, strand, set(transcript_id), set(exon_internal_id))
                
                if transcript_id in transcripts_table.keys():
                    assert(exon_internal_id not in transcripts_table[transcript_id].exons_list) # we should not have mulitple copies of the same exon belong to the same transcript
                    transcripts_table[transcript_id].exons_list.append(exon_internal_id)
                else:
                    transcript_type = "Annotated"
                    transcripts_table[transcript_id] = Transcript(transcript_id, chrom, 0, 0, strand, gene_id, [exon_internal_id], transcript_type, [])
                        
                if exon_internal_id in exons_table.keys():
                    if exon_id != "":
                        exons_table[exon_internal_id].exons_set.add(exon_id)
                    exons_table[exon_internal_id].genes_set.add(gene_id)
                    exons_table[exon_internal_id].trascripts_dict[transcript_id] = exon_number
                else:
                    exons_table[exon_internal_id] = Exon(chrom, start, end, strand, exon_internal_id, set(), set(), {transcript_id : exon_number})
                    exons_table[exon_internal_id].exons_set.add(exon_id)
                    exons_table[exon_internal_id].genes_set.add(gene_id)

    file.close()
    
    return genes_table, transcripts_table, exons_table

def get_transcripts_info(transcripts_table, exons_table):
    types       = dict()
    mono_exon   = 0
    mulit_exons = 0
    
    shortest        = 10000000
    longest         = 0
    shortest_id     = ""
    longest_id      = ""
    transcripts_lst = []
    
    longest_intron     = 0
    shortest_intron    = 10000000
    shortest_intron_id = ""
    longest_intron_id  = ""
    introns_lst        = []
    shortest_exon_1_id = ""
    shortest_exon_2_id = ""
    longest_exon_1_id = ""
    longest_exon_2_id = ""
                   
    for transcript_id, transcript in transcripts_table.items():
        if transcript.type in types.keys():
            types[transcript.type] += 1
        else:
            types[transcript.type] = 1
    
        length = transcript.end - transcript.start 
        
        transcripts_lst.append(length)
        
        if length > longest:
            longest = length
            longest_id = transcript_id
            
        if length < shortest:
            shortest = length
            shortest_id = transcript_id
            
            
        if len(transcript.exons_list) > 1:
            mulit_exons += 1
        else:
            mono_exon += 1

        if len(transcript.exons_list) > 1:
            exons = [exons_table[exon] for exon in transcript.exons_list]
            exons = sorted(exons, key = lambda x: (x.start, x.end))             # this is very important
            
            for i in range(0, len(exons) - 1):
                exon_1 = exons[i]
                exon_2 = exons[i + 1]
                
                intron = exon_2.start - exon_1.end - 1
                
                introns_lst.append(intron)
                
                if intron > longest_intron:
                    longest_intron = intron
                    longest_intron_id = transcript_id
                    longest_exon_1_id = exon_1.ID
                    longest_exon_2_id = exon_2.ID
                    
                if intron < shortest_intron:
                    shortest_intron = intron
                    shortest_intron_id = transcript_id
                    shortest_exon_1_id = exon_1.ID
                    shortest_exon_2_id = exon_2.ID
    
    print("[{}]\t[LOG]\tTranscript information: ".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print("Number of transcripts:             {}".format(len(transcripts_table.keys())))
    print("Shortest transcript:               {}\t{}".format(shortest, shortest_id))
    print("Longest transcript:                {}\t{}".format(longest, longest_id))
    print("Transcripts length mean:           {:.2f}".format(mean(transcripts_lst)))
    print("Transcripts length median:         {}".format(median(transcripts_lst)))
    print("Transcripts length mode:           {}".format(mode(transcripts_lst)))
    
    print("Shortest intron:                   {}\t{}: {} -> {}".format(shortest_intron, shortest_intron_id, shortest_exon_1_id, shortest_exon_2_id))
    print("Longest intron:                    {}\t{}: {} -> {}".format(longest_intron, longest_intron_id, longest_exon_1_id, longest_exon_2_id))
    print("Introns length mean:               {:.2f}".format(mean(introns_lst)))
    print("Introns length median:             {}".format(median(introns_lst)))
    print("Introns length mode:               {}".format(mode(introns_lst)))
    
    print("Number of multi exons transcripts: {}\t{:.2f}%".format(mulit_exons, mulit_exons*100/(mulit_exons + mono_exon)))
    print("Number of mono exon transcripts:   {}\t{:.2f}%".format(mono_exon, mono_exon*100/(mulit_exons + mono_exon)))
    print("\nType of transcripts:")
    
    df = pd.DataFrame(columns=['BioType', 'Count', 'Percentage'])
    
    total = 0
    for key, value in types.items():
        total += value
    
    i = 0
    for key, value in types.items():
        df.loc[i] = [key, value, round(value * 100/total, 2)]
        i += 1

    print(df.sort_values('Percentage', ascending=[False]))

def reverse_complement(sequence):
    reversed_bases = {'N':'N', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    reversed_sequence = "".join([reversed_bases[base] for base in sequence[::-1]])
        
    return reversed_sequence

def collect_junctions_sequences_from_GTF(genome_fasta: Fasta, 
                               transcripts_table: dict, 
                               exons_table: dict, 
                               window: int, 
                               transcript_types: str, 
                               is_verbose  = True,
                               is_training = True,
                               is_strict   = True, 
                               add_junction_seq=False) -> dict:
    
    if transcript_types == "":
        transcript_types = []
    else:
        transcript_types = transcript_types.split(",")
    
    print("[{}]\t[LOG]\tCollecting splice junctions in mode={} and window={}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Strict" if is_strict else "NotStrict", window))
    print("[{}]\t[LOG]\tCollecting splice junctions from transcript types: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), transcript_types if len(transcript_types) != 0 else "All"))
    
    junctions_table = dict()
    w = window // 2
    
    forward_signals = dict()
    reverse_signals = dict()
    
    duplicated_junctions   = 0
    removed_junctions_N    = 0
    removed_junction_short_intron = 0
    removed_junction_short_accepotr = 0
    removed_junction_short_donor = 0
    
    for transcript_id, transcript in transcripts_table.items():
        chrom  = transcript.chrom
        chrom_len = len(genome_fasta[chrom])
        strand = transcript.strand
        transcript_type = transcript.type
        
        if len(transcript.exons_list) > 1 and (transcript_type in transcript_types or len(transcript_types) == 0):
            exons = [exons_table[exon] for exon in transcript.exons_list]   # convert exon hash ID to exon object
            exons = sorted(exons, key = lambda x: (x.start, x.end))         # sort exons by genome coordinates
            
            for i in range(0, len(exons) - 1):
                exon_1 = exons[i]
                exon_2 = exons[i + 1]
                
                signal            = ""
                junction_sequence = ""
                donor_sequence    = ""
                acceptor_sequence = ""
                
                donor_start    = exon_1.start - 1 ## All 0-base
                donor_end      = exon_1.end
                acceptor_start = exon_2.start - 1
                acceptor_end   = exon_2.end
                junction_id    = chrom + "__" + strand + "__" + str(donor_end) + "__" + str(acceptor_start)
                
                donor_sequence_start = donor_end - w 
                donor_sequence_end  = donor_end + w 
                
                acceptor_sequence_start = acceptor_start - w
                acceptor_sequence_end   = acceptor_start + w
                
                if donor_sequence_start < 0:
                    donor_sequence_start = 0
                
                if donor_sequence_end >= chrom_len:
                    donor_sequence_end = chrom_len - 1
                
                if acceptor_sequence_start < 0:
                    acceptor_sequence_start = 0
                
                if acceptor_sequence_end >= chrom_len:
                    acceptor_sequence_end = chrom_len - 1
                
                if is_strict and (acceptor_start - donor_end  < w) :
                    removed_junction_short_intron += 1
                    
                    if is_verbose == True: 
                        print("Junction is removed due to short intron:")
                        print("donor_end:      {}".format(donor_end))
                        print("acceptor_start: {}".format(acceptor_start))
                        print("")
                    
                    if is_training:   ## to add even bad junction to the table when training but still be able to report them as annotated when predicting
                        continue
                    
                elif is_strict and (acceptor_sequence_end - acceptor_sequence_start < w or acceptor_end - acceptor_start < w): 
                    removed_junction_short_accepotr += 1
                    
                    if is_verbose == True: 
                        print("Junction is removed due to short acceptor:")
                        print("acceptor_start:          {}".format(acceptor_start))
                        print("acceptor_end:            {}".format(acceptor_end))
                        print("acceptor_sequence_start: {}".format(acceptor_sequence_start))
                        print("acceptor_sequence_end:   {}".format(acceptor_sequence_end))
                        print("")
                        
                    if is_training:   ## to add even bad junction to the table when training but still be able to report them as annotated when predicting
                        continue
                
                elif is_strict and (donor_sequence_end - donor_sequence_start < w  or donor_end - donor_start < w) :
                    removed_junction_short_donor += 1
                    
                    if is_verbose == True: 
                        print("Junction is removed due to short donor:")
                        print("donor_start:          {}".format(donor_start))
                        print("donor_end:            {}".format(donor_end))
                        print("donor_sequence_start: {}".format(donor_sequence_start))
                        print("donor_sequence_end:   {}".format(donor_sequence_end))
                        print("")
                        
                    if is_training:       ## to add even bad junction to the table when training but still be able to report them as annotated when predicting
                        continue
                
                else:
                    if is_training:
                        signal              = str(genome_fasta[chrom][donor_end : donor_end + 2]).upper() + \
                                              str(genome_fasta[chrom][acceptor_start - 2 : acceptor_start]).upper()
                        donor_sequence      = str(genome_fasta[chrom][donor_sequence_start : donor_sequence_end]).upper()
                        acceptor_sequence   = str(genome_fasta[chrom][acceptor_sequence_start : acceptor_sequence_end]).upper()
                        junction_sequence   = str(genome_fasta[chrom][donor_sequence_start : donor_sequence_end]).upper() + \
                                              str(genome_fasta[chrom][acceptor_sequence_start : acceptor_sequence_end]).upper()
                        
                        for letter in junction_sequence + donor_sequence + acceptor_sequence: 
                            if letter not in 'ACTGactg':
                                break
                        
                        dna_alphabet = [base for base in "ACTGactg"]
                
                        if any(base not in dna_alphabet for base in junction_sequence + donor_sequence + acceptor_sequence):
                            removed_junctions_N += 1
                        
                            if is_verbose == True: 
                                print("Junction is removed due to character not ACTG: {}".format(junction_id))
                                print(junction_sequence)
                                print(donor_sequence)
                                print(acceptor_sequence)
                            
                            signal            = ""          ## to add even bad junction to the table when training but still be able to report them as annotated when predicting
                            junction_sequence = ""
                            donor_sequence    = ""
                            acceptor_sequence = ""
                        
                            continue

                    else:
                        signal              = str(genome_fasta[chrom][donor_end : donor_end + 2]).upper() + \
                                              str(genome_fasta[chrom][acceptor_start - 2 : acceptor_start]).upper()
                        junction_sequence = ""
                        donor_sequence    = ""
                        acceptor_sequence = ""
                        
                # Generate kmers from sequences          
                donor_ID    = exon_1.ID
                acceptor_ID = exon_2.ID
                
                transcript.junctions.append(junction_id)
                    
                if junction_id not in junctions_table.keys() or junctions_table[junction_id].signal == "":
                    donor_sequence = donor_sequence.upper()
                    acceptor_sequence = acceptor_sequence.upper()
                    junction_sequence = junction_sequence.upper()
                    signal = signal.upper()
                    
                    if strand == '-':
                        signal            = reverse_complement(signal)

                        if is_training:
                            temp              = donor_sequence
                            
                            junction_sequence = reverse_complement(junction_sequence)
                            donor_sequence    = reverse_complement(acceptor_sequence)
                            acceptor_sequence = reverse_complement(temp)
                        
                        if signal in reverse_signals.keys():
                            reverse_signals[signal] += 1
                        else:
                            reverse_signals[signal] = 1
                    
                    else:
                        if signal in forward_signals.keys():
                            forward_signals[signal] += 1
                        else:
                            forward_signals[signal] = 1
                    
                    if add_junction_seq == False:
                        junction_sequence = ''

                    junctions_table[junction_id] = Junction(chrom, donor_ID, acceptor_ID, donor_end, acceptor_start, strand, signal, transcript_type, donor_sequence, acceptor_sequence, junction_sequence)
                    
                else:
                    duplicated_junctions += 1

    signals = sorted(list(set(list(forward_signals.keys()) + list(reverse_signals.keys()))))

    df = pd.DataFrame(columns=['Signal', 'Forward', 'Reverse', 'Percentage'])
    total = 0
    for signal in signals:
        if signal not in forward_signals.keys():
            forward_signals[signal] = 0
        
        if signal not in reverse_signals.keys():
            reverse_signals[signal] = 0
        
        total +=  forward_signals[signal] + reverse_signals[signal]
    
    i = 0
    for signal in signals:    
        df.loc[i] = [signal, forward_signals[signal], reverse_signals[signal], round((forward_signals[signal] + reverse_signals[signal]) * 100/total, 2)]
        i += 1

    print("Number of duplicated junctions:        {}".format(duplicated_junctions))
    print("Number of short junctions (intron):    {}".format(removed_junction_short_intron))
    print("Number of short junctions (donor):     {}".format(removed_junction_short_donor))
    print("Number of short junctions (acceptor):  {}".format(removed_junction_short_accepotr))
    print("Number of junctions contains N:        {}".format(removed_junctions_N))
    print("Number of accepted junctions:          {}".format(len(junctions_table)))

    if is_verbose == True:
        print("Splicing Signals Types: ")
        print(df.sort_values('Percentage', ascending=[False]).to_string(index=False))
    else:
        print("The First 10 Splicing Signals Types: ")
        print(df.sort_values('Percentage', ascending=[False]).head(10).to_string(index=False))
    return junctions_table 


def resolove_overlapping_exons(temp_exons_loci, regions_loci):
    # print("\nExons:  \t", temp_exons_loci)
    
    start_loci = set()
    end_loci   = set()
    
    for temp_exon in temp_exons_loci:
        start_loci.add(temp_exon[0] - 1)                         ## All regions are 0-based
        end_loci.add(temp_exon[1])
    
    # print("st-Loci: \t", start_loci)
    # print("en-Loci: \t", end_loci)
    
    loci = sorted(list(start_loci) + list(end_loci) )
    
    # print("Loci:   \t", loci)
    for i in range(len(loci) - 1) :
        # print("\tnew_regions: ", [loci[i], loci[i + 1]])
        regions_loci.append( [loci[i], loci[i + 1]] ) # - 1, loci[i + 1]] - 1)
    
    # print("Regions:\t", regions_loci)
    
def get_regions(chrom:str, chrom_len:int,
                exons_table:dict, transcripts_table:dict,
                regions_loci_table:dict, 
                regions_array_table:dict,
                regions_relations_table:dict,
                regions_transcripts_table:dict,
                jumps_transcripts_table:dict,
                lock:threading.Lock):
    """
    chromosome  |---------------------------------------------------------------------------------------------------------------|
    Genes                  |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|         |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|   
    Transcript             |-----------|             |---------|   |--------|         
    Transcript                |--------|           |-----------|   |--------|         
    Transcript                                                                           |----|                    |----|
    Transcript                                                                        |-------|   |-------|     |-------| 
    Regions     |    0     |1 |   2    |     3     |4|    5    | 6 |   7    |    8    |  9 |10| 11|  12   | 13  |14| 15 |   16   |
    """
    
    # 1- Collecting regions from exons boundries
    exons_loci = list()                     
    for exon in exons_table.values():
        if exon.chrom == chrom:
            exon_locus = [int(exon.start), int(exon.end)] # all exons are 1-based
            if exon_locus not in exons_loci:
                exons_loci.append(exon_locus)
    
    if len(exons_loci) == 0:
         return chrom, 0
    
    exons_loci = sorted(exons_loci, reverse=False, key = lambda x: (x[0], x[1]))
    
    temp_exons_loci = [[exons_loci[0][0], exons_loci[0][1]]]
    temp_exon_start = exons_loci[0][0]
    temp_exon_end   = exons_loci[0][1]
    regions_loci    = list()  ## All regions are 0-based

    index = 0
    for exon_start, exon_end in exons_loci[1:]:
        index += 1
        if exon_start > temp_exon_end:
            if len(temp_exons_loci) == 1:
                regions_loci.append([temp_exon_start - 1, temp_exon_end]) ## All regions are 0-based
            else:
                resolove_overlapping_exons(temp_exons_loci, regions_loci)
                
            temp_exons_loci = [[exon_start, exon_end]]
            temp_exon_start = exon_start
            temp_exon_end   = exon_end

        else :
            if exon_end > temp_exon_end:  # very imoportant condition
                temp_exon_end = exon_end
                
            temp_exons_loci.append([exon_start, exon_end])
        
        if index == len(exons_loci) - 1:
            if len(temp_exons_loci) == 1:
                regions_loci.append([temp_exon_start - 1, temp_exon_end]) ## All regions are 0-based
            else:
                resolove_overlapping_exons(temp_exons_loci, regions_loci)

    # print("Final regions_loci: ", regions_loci)
    
    # 2- Filling up the regions_array with regions indices
    regions_array = np.empty(chrom_len, dtype=int)
    regions_array.fill(-1)
    
    i = 0
    for region_start, region_end in regions_loci:
        regions_array[region_start:region_end] = i # very import destinction, regions with value -1 are introns !!
        i += 1
         
    # 3- Filling up regions info tables
    regions_relations   = dict()    # splice-junction relationships between regions
    regions_transcripts = dict()    # for each region, get all the trranscripts that has it
    jumps_transcripts = dict()
    
    for transcript_id, transcript in transcripts_table.items():
        if transcript.chrom == chrom:
            exons  = [exons_table[exon] for exon in transcript.exons_list]       
            exons  = sorted(exons, key = lambda x: (x.start, x.end)) 
            strand = transcript.strand
            
            for i in range(0, len(exons) - 1):
                exon_1 = exons[i]
                exon_2 = exons[i + 1]
                
                regions_1_first = regions_array[exon_1.start - 1]   ## exons are 1-based, while all regions are 0-based
                regions_1_last  = regions_array[exon_1.end - 1]

                regions_2_first = regions_array[exon_2.start - 1]   ## exons are 1-based, while all regions are 0-based
                regions_2_last  = regions_array[exon_2.end - 1]
                
                ## Filling up regions_transcripts
                for region in range(regions_1_first, regions_1_last + 1):
                    
                    assert(region != -1)
                    
                    if region not in regions_transcripts.keys():
                        regions_transcripts[region] = set()
                    regions_transcripts[region].add(transcript_id)
                
                for region in range(regions_2_first, regions_2_last + 1):
                    assert(region != -1)
                    if region not in regions_transcripts.keys():
                        regions_transcripts[region] = set()
                    regions_transcripts[region].add(transcript_id)

                ## Filling up regions_relations
                donor_region    = regions_1_last
                acceptor_region = regions_2_first
                
                if strand == '-':
                    temp = donor_region
                    donor_region = acceptor_region
                    acceptor_region = temp
                    
                if donor_region not in regions_relations.keys():
                    regions_relations[donor_region] = set()
                    
                regions_relations[donor_region].add(acceptor_region)
                
                ## Filling up jumps_transcripts
                jump_id = str(donor_region) + "__" + str(acceptor_region) 
                
                if jump_id not in jumps_transcripts.keys():
                    jumps_transcripts[jump_id] = set()
                        
                jumps_transcripts[jump_id].add(transcript_id)
                
                if strand not in ['+', '-']:
                    print("Problematic transcript with strand issues")
                    print(transcript)
                    exit()
                    
    # 4- Thread-safe updating of the regions tables
    if len(regions_loci) > 2:
        with lock:
            regions_loci_table[chrom]        = regions_loci
            regions_array_table[chrom]       = regions_array
            regions_relations_table[chrom]   = regions_relations
            regions_transcripts_table[chrom] = regions_transcripts
            jumps_transcripts_table[chrom]   = jumps_transcripts
        
    return chrom, len(regions_loci)

def generate_regions(genome_sequences_table:dict, transcripts_table:dict, exons_table:dict, num_threads:int):
    """Multi-threaded function to generate regions informations 
    Args:
        genome_sequences_table (dict): Dictionary contains the genome seuence for each chromosome
        transcripts_table (dict): Dictionary contains transcripts info of type Transcript
        exons_table (dict): Dictionary contains exons info of type Exon
    Returns:
        tuple : four dictionaries of regions_loci_table , regions_array_table, regions_relations_table, regions_transcripts_table
    """
    regions_loci_table = dict()            # For each chromosome, has an array of the regions boundires 
    regions_array_table = dict()            # For each chromosome, has an array of the size of the chromosome storing the region for that base, act like fast access datastructure to get the region of any locus in the a chromosome
    regions_relations_table   = dict()      # For each chromosome, has a dictionary for each region to store the other regions that has splice-jucntion reglationship with this region
    regions_transcripts_table = dict()      # For each chromosome, has a dictionary for each region to store ALL the transcripts that has this region
    jumps_transcripts_table   = dict()
    
    regions_table_lock = threading.Lock()

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit tasks to the executor and gather the results
        tasks = [executor.submit(get_regions,
                                 chrom, 
                                 len(genome_sequences_table[chrom]),
                                 exons_table,
                                 transcripts_table,
                                 regions_loci_table, 
                                 regions_array_table,
                                 regions_relations_table,  
                                 regions_transcripts_table,
                                 jumps_transcripts_table,
                                 regions_table_lock) for chrom in genome_sequences_table.keys()]
        
        # Retrieve results from the completed tasks
        for future in concurrent.futures.as_completed(tasks):
            result = future.result()
            # print(f"\tChromosome, n_regions: {result}")
    
    return regions_loci_table , regions_array_table, regions_relations_table, regions_transcripts_table, jumps_transcripts_table

def get_region_by_locus(chrom, locus, regions_array_table):
    
    region = None 

    if chrom in regions_array_table.keys():
        try:
            region = regions_array_table[chrom][locus]
        except:
            print(f"locus: {locus}, new locus {locus - 1}")
            region = regions_array_table[chrom][-1]
    
    return region

def get_transcripts_by_region(chrom, region, regions_transcripts_table):
    if chrom in regions_transcripts_table.keys():
        if region in regions_transcripts_table[chrom].keys():
            return regions_transcripts_table[chrom][region]
        
    return set()

def is_spliced_relationship(chrom, donor_region, acceptor_region, regions_relations_table):
    is_spliced = False
    
    if chrom in regions_relations_table.keys():
        if donor_region in regions_relations_table[chrom].keys():
            if acceptor_region in regions_relations_table[chrom][donor_region]:
                is_spliced = True

    return is_spliced

def are_mates_from_same_transcript(chrom, mate_1_is_novel, mate_2_is_novel, mate_1_regions, mate_2_regions, 
                                   regions_relations_table, regions_transcripts_table, jumps_transcripts_table):
    
    # print("mate_1_regions: ", mate_1_regions)
    # print("mate_2_regions: ", mate_2_regions)
    
    if mate_1_is_novel or mate_2_is_novel:
        # print("\nNovel Junctions")
        return True                         # Already penalized
    
    if len(mate_1_regions) == 0 or len(mate_2_regions) == 0: # Chromosome with no regions, no transcripts reported in this chromosome
        # print("\nRegions 0")
        return True 

    if -1 in mate_1_regions or -1 in mate_2_regions:
        # print("\nIntronic")
        return False 
    
    mate_1_regions_lst = sorted(list(mate_1_regions))
    mate_1_regions_trascripts = set()
    
    # print("\nmate_1_regions_lst: ", mate_1_regions_lst)
    
    for i in range(len(mate_1_regions_lst)):
        next_transcripts = get_transcripts_by_region(chrom, mate_1_regions_lst[i], regions_transcripts_table)
        # print("next-transcripts:\t", next_transcripts)
        
        if i == 0:
            mate_1_regions_trascripts = next_transcripts
            # print("start-intersct:\t", mate_1_regions_trascripts)
        else:
            mate_1_regions_trascripts = mate_1_regions_trascripts.intersection(next_transcripts)
            # print("After-intersct:\t", mate_1_regions_trascripts)
        
        if i + 1 < len(mate_1_regions_lst):  # make sure that the next region is not adjacent region 
            donor_region     = mate_1_regions_lst[i] 
            acceptor_region  = mate_1_regions_lst[i + 1]
            jump_transcripts = set()
            is_spliced       = False
            
            if is_spliced_relationship(chrom, donor_region, acceptor_region, regions_relations_table):
                is_spliced = True
                jump_transcripts = jump_transcripts.union(jumps_transcripts_table[chrom][str(donor_region) + "__" + str(acceptor_region)])
                # print("1-jump-transcripts:\t", jump_transcripts)
            
            if is_spliced_relationship(chrom, acceptor_region, donor_region, regions_relations_table):
                is_spliced = True
                jump_transcripts = jump_transcripts.union(jumps_transcripts_table[chrom][str(acceptor_region) + "__" + str(donor_region)])
                # print("2-jump-transcripts:\t", jump_transcripts)
        
            if is_spliced :
                mate_1_regions_trascripts = mate_1_regions_trascripts.intersection(jump_transcripts)
                # print("After-jump-intersct:\t", mate_1_regions_trascripts)
            
    # print("mate_1_regions_trascripts: ", mate_1_regions_trascripts)
         
    
    mate_2_regions_lst = sorted(list(mate_2_regions))
    mate_2_regions_trascripts = set()
    
    # print("\nmate_2_regions_lst: ", mate_2_regions_lst)
    
    for i in range(len(mate_2_regions_lst)):
        next_transcripts = get_transcripts_by_region(chrom, mate_2_regions_lst[i], regions_transcripts_table)
        # print("next-transcripts:\t", next_transcripts)
        
        if i == 0:
            mate_2_regions_trascripts = next_transcripts
            # print("start-intersct:\t", mate_2_regions_trascripts)
        else:
            mate_2_regions_trascripts = mate_2_regions_trascripts.intersection(next_transcripts)
            # print("After-intersct:\t", mate_2_regions_trascripts)
        
        if i + 1 < len(mate_2_regions_lst):  # make sure that the next region is not adjacent region 
            donor_region     = mate_2_regions_lst[i] 
            acceptor_region  = mate_2_regions_lst[i + 1]
            jump_transcripts = set()
            is_spliced       = False
            
            if is_spliced_relationship(chrom, donor_region, acceptor_region, regions_relations_table):
                is_spliced = True
                jump_transcripts = jump_transcripts.union(jumps_transcripts_table[chrom][str(donor_region) + "__" + str(acceptor_region)])
                # print("1-jump-transcripts:\t", jump_transcripts)
            
            if is_spliced_relationship(chrom, acceptor_region, donor_region, regions_relations_table):
                is_spliced = True
                jump_transcripts = jump_transcripts.union(jumps_transcripts_table[chrom][str(acceptor_region) + "__" + str(donor_region)])
                # print("2-jump-transcripts:\t", jump_transcripts)
        
            if is_spliced :
                mate_2_regions_trascripts = mate_2_regions_trascripts.intersection(jump_transcripts)
                # print("After-jump-intersct:\t", mate_2_regions_trascripts)
            
    # print("mate_2_regions_trascripts: ", mate_2_regions_trascripts)
    
        
    if len(mate_1_regions_trascripts.intersection(mate_2_regions_trascripts)) > 0:
        return True
    else:
        return False

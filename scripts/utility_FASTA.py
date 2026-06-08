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


import numpy as np
import subprocess
import gzip
from datetime import datetime

def parse_fasta(path:str):
    """Parsing FASTA file

    Args:
        path (string): Path to FASTA file

    Returns:
        sequences (dictionary): dictionary of seauence name to string of bases
    """
    
    print("[{}]\t[LOG]\tParsing FASTA file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), path))
    
    sequences = dict()
    buffer = ""
    header = ""

    file   = open(path, "r")
    zipped = False
    
    if path[-3:] == ".gz":
        file = gzip.open(path, 'rb')
        zipped = True
        
    while True:
        line = file.readline()
        
        if len(line) == 0:
            break
        
        if zipped:
            line = line.decode()
        
        if line[0] == '>':
            if len(buffer) != 0:
                sequences[header] = buffer
            
            buffer = ""
            if ' ' in line:
                header = line.split(' ')[0][1:]
            else:
                header = line[1:-1]
        else:
            buffer += line[:-1]
    
    if len(buffer) != 0:
        sequences[header] = buffer
                
    file.close()
    return sequences 

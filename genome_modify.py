###
# Copyright 2016 Miler T. Lee, University of Pittburgh
# This file is part of the RESA Suite
#
# RESA Suite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# RESA Suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RESA Suite.  If not, see <http://www.gnu.org/licenses/>.
#
#
# genome_modify.py:
# Inserts leading restriction site sequences corresponding to
# UTR primers that have overhang. Facilitates alignment
# of reads from the edges
#
##

"""
"""

import subprocess
import sys
import re
import resa_util as ru

sides = ['L', 'R']


def case_string_split(input_str):
    """
    splits the string along the first case switch
    and returns a tuple
    aaaaBBBBBxxx  -->  (aaaa, BBBBBxxx)
    """

    index = len(input_str)
    for nt in ['A', 'C', 'G', 'T', 'N']:
        pos = input_str.find(nt)
        if pos > -1:
            index = min(index, pos)

    return input_str[:index], input_str[index:]


def genome_modify(primer_file, utr_bedfile, genome_fasta, fasta_outfile):
    chrs = dict(ru.read_fasta(genome_fasta))
    modified_chrs = {}

    print('ID side strand restriction_site genomic_site primer genomic_primer new_sequence_window')

    #Store the utr positions
    utrs = {}
    f = open(utr_bedfile)
    for line in f:
        fields = line.strip().split()
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        utr_id = fields[3]
        strand = fields[5]
        utrs[utr_id] = (chrom, start, end, strand)
    f.close()

    #Iterate through the primers
    f = open(primer_file)
    header = f.readline()
    for line in f:
        fields = line.strip().split()
        utr_id = fields[0] + '_' + fields[1]
        
        try:
            chrom, start, end, strand = utrs[utr_id]
        except KeyError:
            continue
        
        if chrom not in modified_chrs:
            chr_seq = chrs[chrom].upper()
            modified_chrs[chrom] = 1
        else:
            chr_seq = chrs[chrom]
        
        primers = fields[3], fields[4]
        for i, primer in enumerate(primers):
            restr, primer_seq = case_string_split(primer)
            if not restr:
                continue
            
            if (i == 0 and strand == '+') or (i == 1 and strand == '-'):
                primer_coord = start
                restr_coord = start - len(restr)
                seq_to_replace = chr_seq[restr_coord:start]
                genomic_primer = chr_seq[primer_coord:primer_coord+len(primer_seq)]

                #Make the modification
                chr_seq = chr_seq[:restr_coord] + restr + chr_seq[start:]
                new_seq_window = chr_seq[(restr_coord - 20): (primer_coord+25)]
            elif (i == 0 and strand == '-') or (i == 1 and strand == '+'):
                primer_seq = ru.rc(primer_seq)
                restr = ru.rc(restr)

                primer_coord = end - len(primer_seq)
                restr_coord = end
                seq_to_replace = chr_seq[restr_coord:restr_coord+len(restr)]
                genomic_primer = chr_seq[primer_coord:end]

                #Make the modification
                chr_seq = chr_seq[:restr_coord] + restr + chr_seq[restr_coord+len(restr):]
                new_seq_window = chr_seq[primer_coord:(restr_coord+25)]

            chrs[chrom] = chr_seq
            print(utr_id, sides[i], strand, restr, seq_to_replace, primer_seq, genomic_primer, new_seq_window)

        
    #Output the chromosomes back to a fasta file
    fw = ru.fasta_writer(fasta_outfile)
    for chrom in sorted(modified_chrs.keys()):
        fw.write(chrom, chrs[chrom].upper())
    fw.close()


    
if __name__ == "__main__":
    if sys.argv[1] == '-h': #test code
        print('Usage: python genome_modify.py <primer_file> <utr_bed> <genome_fasta> <out_fasta>')
        sys.exit()
        
    genome_modify(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

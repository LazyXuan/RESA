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
# resa_util.py: Utilities for processing RESA data
##

import bz2
import gzip
import json
import re
import sys
import subprocess


def initialize_loci(utr_bed12_file, utr_fasta_file, test = False):
    """
    Given the UTRs listed in the bed12 file and
    their corresponding sequences in the fasta_file,
    creates a dict of loci[key] = (chr, strand, exon_list, seq)
    """
    seqs = dict(read_fasta(utr_fasta_file))

    loci = {}
    
    f = open(utr_bed12_file)
    for line in f:
        fields = line.strip().split()
        chrom = fields[0]
        start = int(fields[1])
        strand = fields[5]
        feat_id = fields[3]
        
        block_sizes = fields[10].strip(',').split(',')
        block_starts = fields[11].strip(',').split(',')

        exons = []
        for i, (bsize, bstart) in enumerate(zip(block_sizes, block_starts)):
            gstart = start + int(bstart)
            gend = gstart + int(bsize)
            exons.append((gstart, gend))
        loci[feat_id] = (chrom, strand, tuple(exons), seqs[fields[3]].upper())
        if test:
            break
    f.close()
    return loci



###
# UTILITIES
###

nt_mutations = {'C': 'T', 'G': 'A'}
anti_strand_str = {'-': '+', '+': '-'}


###string.maketrans('acgturyACGTURY', 'tgcaayrTGCAAYR')
DNA_TRANS = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQYSAAVWXRZ[\\]^_`tbgdefchijklmnopqysaavwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'


def rc(sequence, reverse = True):
    """
    Reverse complement a DNA sequence, preserving case
    """

    result = sequence.translate(DNA_TRANS)
    if reverse:
        return result[::-1]
    else:
        return result


def bam_entry_is_reverse(samflag):
    """
    Flag is passed in as an integer. Determines
    whether the 0x10 bit is set (16 in base 10),
    which indicates reverse complemented sequence.
    This is done using a binary operator &
    """
    
    return samflag & 16 == 16


def seq_mask(seq, chars = ['A', 'G']):
    """
    Replaces specified characters with N
    """
    for char in chars:
        seq = seq.replace(char, 'N')
    return seq


def load_chr_seqs(genome_fa):
    """
    Loads all chromosome sequences into a dict
    """
    chr_dict = dict(read_fasta(genome_fa))
    return chr_dict


def load_chr_seq(chr_id, chr_dict, genome_fa):
    """
    Loads the chromosome sequence into memory if it's not
    already there
    """
    if chr_id not in chr_dict:
        fasta_file = genome_fa % chr_id
        chr_dict[chr_id] = read_fasta(fasta_file)[0][1]
    return chr_dict[chr_id]


def decode_cigar(cigar):
    """
    Parses the cigar string into integers and letters
    """
    return re.findall('(\d+)([MNDISHPX=])', cigar)


def cigar_span_(cigar):
    """
    Interprets the cigar string as the number of genomic
    positions consumed
    """
    span = 0
    cigar_ops = decode_cigar(cigar)
    for nts, op in cigar_ops:
        nts = int(nts)

        if op != 'I':
            span += nts
    return span


def cigar_span(cigar):
    return sum(int(x) for x in re.findall('(\d+)[MNDSHPX=]', cigar)) #no I


def tx_indexing(exons, minus = False, inverse = False):
    """
    Returns a dict of genomic coordinates -> tx coordinates
    (or the inverse if inverse = True)
    Exons are zero indexed.
    """
    
    positions = []
    for s, e in exons:
        positions += [i for i in range(s, e)]
        
    if minus:
        positions.reverse()
            
    if inverse:
        return {i:x for i, x in enumerate(positions)}
    else:
        return {x:i for i, x in enumerate(positions)}

    
def pretty_str(x, fields = False):
    """
    Handles tuples or lists
    """
    def joined_string(x, sep=','):
        return sep.join(list([str(y) for y in x]))

    if isinstance(x, str):
        return x
    elif isinstance(x, float):
        if abs(x) < 0.001:
            return '%.1E' % x
        else:
            return '%.3f' % x
    elif isinstance(x, tuple) or isinstance(x, list):
        if fields:
            return joined_string(x, '\t')
        elif not x:
            return '.'
        elif isinstance(x[0], tuple) or isinstance(x[0], list):
            return ';'.join([joined_string(y) for y in x])
        else:
            return joined_string(x)
    else:
        return str(x)


#######################
# FASTA file processing
#######################

def read_fasta(filename):
    """
    Returns the contents of a fasta file in a list of (id, sequence)
    tuples. Empty list returned if there are no fasta sequences in the file
    """

    a = fasta_reader(filename)
    seqs = []
    while a.has_next():
        seqs.append(next(a))
    return seqs

                
class fasta_reader:
    """
    Lightweight class for incrementally reading fasta files.

    Supports reading directly from properly named
    gzipped (.gz or .z) or bzip2ed (.bz2) files.
    """

    file = None
    nextheader=''
    
    def __init__(self, filename):
        try:
            if filename.endswith('.gz') or filename.endswith('.z'):
                self.file = gzip.open(filename, 'rb')
            elif filename.endswith('.bz2'):
                self.file = bz2.BZ2File(filename, 'rb')
            else:
                self.file = open(filename, 'r')

            # fast forward to the first entry
            while 1:
                line = self.file.readline()
                if line == '':
                    self.close()
                    return
                elif line[0] == '>':
                    self.nextheader = line[1:].rstrip()
                    return
        except IOError:
            #print('No such file', filename)
            raise
                    
    def has_next(self):
        """
        Returns true if there are still fasta entries
        """
        
        return len(self.nextheader) > 0

    def __next__(self):
        """
        Returns an (id, sequence) tuple, or () if file is finished
        """

        #if global nextheader is empty, return empty
        #otherwise, the header is the nextheader
        try:
            identifier = self.nextheader
            total = []
            while 1:
                line = self.file.readline()
                if line == '' or line[0] == '>':  #EOF, end of entry
                    break
                total.append(line.rstrip())

            sequence = ''.join(total)

            if len(line) > 0:
                self.nextheader = line[1:].rstrip()
            else:
                self.nextheader = ''
                self.close()

            return (identifier, sequence)

        except:
            self.nextheader=''
            self.close()
            return ()

    def close(self):
        """
        Close the fasta file
        """
        
        self.file.close()


def write_fasta(filename, id_or_list, seq='', width=60, gzip_compress = False):
    """
    Writes a fasta file with the sequence(s)
    version 1: write_fasta(myfilename, 'seq1_id', 'AAAAA')
    version 2: write_fasta(myfilename, [('seq1_id', 'AAAAA'),
    ('seq2_id', BBBBB)])
    """

    a = fasta_writer(filename, width=width, gzip_compress = gzip_compress)
    a.write(id_or_list, seq)
    a.close()
    

class fasta_writer:
    """
    Rudimentary fasta file writer

    Supports writing out to a gzipped file. If the passed in filename
    does not end with .gz or .z, .gz is appended.
    """

    file = None
    width = 0
    
    def __init__(self, filename, width=60, gzip_compress = False):
        self.width = width
        try:
            if gzip_compress:
                if not filename.endswith('.gz') and not filename.endswith('.z'):
                    filename += '.gz'
                self.file = gzip.open(filename, 'wb')
            else:            
                self.file = open(filename, 'w')
        except IOError:
            print('Can\'t open file.')

    def write(self, id, seq=''):
        """
        Supports an id and a sequence, an (id, seq) tuple, or
        a list of sequence tuples
        """

        if type(id) == type([]):
            list(map(self.writeone, id))
        else:
            self.writeone(id, seq)

    def writeone(self, id, seq=''):
        """
        Internal method.
        """

        if type(id) == type((0,0)):
            seq = id[1]
            id = id[0]

        line_width = self.width
        if self.width == 0:
            line_width = len(seq)
        self.file.write(">" + id + "\n")
        i = 0

        while i < len(seq):
            self.file.write(seq[i:i+line_width] + "\n")
            i+=line_width
            
    def close(self):
        """
        Closes the fasta file.
        """

        self.file.close()


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
# resa_coverage.py: Data preparation for RESA analysis
##


import resa_util as ru
import pickle
import re
import subprocess
import sys
import numpy as np
import random
import argparse
import os
from datetime import datetime

CHR_SEQS = {} #For lazy loading

SAMFLAGS_PAIRED = {'+': ('99', '147'),
                   '-': ('163', '83')
                  }

SAMFLAGS_SINGLE = {'+': ('0',),
                   '-': ('16',)
                  }

FSAMFLAGS_STRICT_PAIRED = {'+': '3840',
                           '-': '3840'}

FSAMFLAGS_STRICT_SINGLE = {'+': '3856',
                           '-': '3840'}

FSAMFLAGS_LOOSE_PAIRED = {'+': '0',
                           '-': '0'}

FSAMFLAGS_LOOSE_SINGLE = {'+': '16',
                           '-': '0'}

SAMFLAGS = SAMFLAGS_PAIRED
FSAMFLAGS = FSAMFLAGS_STRICT_PAIRED

###############################
# Construct names of bamfiles produced through resa_map
### 

TOPHAT_DIR = '' #'tophat-utr'
WT_BAM = 'wt/accepted_hits.bam'
#C2T,G2A
MUT_BAMS = {'+': ['c2t_reads/c2t_tophat/accepted_hits_converted.bam',
                  'g2a_reads/g2a_tophat/accepted_hits_converted.bam'],
            '-': ['c2t_reads/g2a_tophat/accepted_hits_converted.bam',
                  'g2a_reads/c2t_tophat/accepted_hits_converted.bam']
            }

def process_sample_lists(sample_lists, labels, mut = False, samples_as_bamfiles = False):
    """
    samples[sample_label][strand] = [bam1, bam2, ...]
    """
    
    labels = labels.split(',')
    if len(labels) != len(sample_lists):
        sys.exit('**ERROR: Number of sample lists does not match number of labels, exiting.')

    samples = {}
    for label, sample_list in zip(labels, sample_lists):
        samples[label] = {'+': [], '-': []}
        
        paths = sample_list.split(',')
        for path in paths:
            if mut:
                for strand in '+', '-':
#                    samples[label][strand] += [os.path.join(path, TOPHAT_DIR, x) for x in MUT_BAMS[strand]]
                    samples[label][strand].append(tuple([os.path.join(path, TOPHAT_DIR, x) for x in MUT_BAMS[strand]]))
            else:
                for strand in '+', '-': #Same bamfile for both + and -
                    if samples_as_bamfiles:
                        full_path = path
                    else:
                        full_path = os.path.join(path, TOPHAT_DIR, WT_BAM)
                    samples[label][strand].append((full_path,))
                    
    #Check that all the bamfiles exist
    for bamdict in samples.values():
        for bamlist in bamdict.values():
            for bamtuple in bamlist:
                for bam in bamtuple:
                    if not os.path.exists(bam):
                        sys.exit('**ERROR: Bamfile %s does not exist, exiting.' % bam)
    
    return samples


#################################
# COVERAGE
### 

def bam_read_pairs(bamfile, coord, samflags, Fsamflag = 0):
    """
    Consecutive shell calls for each samflag.
    """
    read_hash = {}
    for samflag in samflags:
        cmd = ['samtools', 'view', '-f', samflag, '-F', Fsamflag, bamfile, coord]
        p = subprocess.Popen(cmd, bufsize = 1, stdout = subprocess.PIPE, universal_newlines = True)
        for line in p.stdout:
            fields = line.split('\t')
            read_id = fields[0]
            start = int(fields[3]) - 1 #Convert to 0-index
            pair_start = int(fields[7])-1
            cigar = fields[5]
            read_seq = fields[9]

            if start < pair_start:
                key = (read_id, start, pair_start)
            else:
                key = (read_id, pair_start, start)
                
            if key not in read_hash:
                read_hash[key] = [(start, cigar, read_seq)]
            else:
                read_hash[key].append((start, cigar, read_seq))
                if start < pair_start:
                    sys.exit('***Error: read pairs not in sorted order, exiting.')
                
    return read_hash


def tx_pileup(bamfile_list, coord, tx_coords, samflags = None, Fsamflag = 0, dup_limit = 0, strict_dup_filtering = False, raw_frags = False, frag_sizes = (0,1000), sampling_rate = 1):
    """
    Positional coverage, including the insert between read pairs
    Performed using feature-relative coordinates rather than genomic.
    As such, a conversion table must be provided tx_coords[genomic_index] = tx_index.

    Read pairs with coordinates beyond the extent of the feature are truncated.
    Coordinates inside the span of the feature but with no equivalency belong
    to a different splice isoform and are not counted.

    If dup_limit > 0, only dup_limit reads with start or end are allowed to be counted

    Return is a dict: coverage[relative_position] = read_count
    """

    if sampling_rate < 1:
        raw_frags = True
    
    min_g_index = min(tx_coords.keys())
    max_g_index = max(tx_coords.keys())
    max_tx_index = max(tx_coords.values())

    coverage = [0] * (max_tx_index+1)
    read_starts = [0] * (max_tx_index+1)
    read_ends = [0] * (max_tx_index+1)
    pcr_dups = {}    
    frags = []
    
    for bamfile in bamfile_list:
        read_hash = bam_read_pairs(bamfile, coord, samflags, Fsamflag)

        #Order matters when strict filtering
        if strict_dup_filtering:
            read_list = sorted(read_hash.values())
        else:
            read_list = read_hash.values()
        for reads in read_list:
            #Sorting is probably unnecessary
            reads.sort()
            
            gstart = reads[0][0]
            gend = reads[-1][0] + ru.cigar_span(reads[-1][1])
            
            if gstart < min_g_index:
                gstart = min_g_index
            if gend > max_g_index:
                gend = max_g_index

            try:
                start = tx_coords[gstart]
                end = tx_coords[gend]
            except KeyError:
                continue

            if end < start:
                start, end = end, start
                
            #Sanity filter: if the inferred span is smaller than the size of a read,
            #then it belongs to a different splice form and should not be counted
            #A buffer is provided to account for indel variants
            tx_span = end - start
            if tx_span < len(reads[0][2]) - 5:
                continue

            #Filter on the specified frag_sizes if provided
            if tx_span < frag_sizes[0] or tx_span > frag_sizes[1]:
                continue

            if raw_frags:
                frags.append((start, end))
            else:
                if not dup_limit or pcr_dups.get((start, end), 0) < dup_limit:
                    if strict_dup_filtering and (read_starts[start] >= dup_limit or read_ends[end] >= dup_limit):
                        continue #ignore the fragment

                    for i in range(start, end+1):
                        coverage[i] += 1
                        
                    #Strict dup filtering depends on the dict order
                    #so results may not match exactly if that changes
                    if dup_limit:
                        pcr_dups[(start, end)] = pcr_dups.get((start, end), 0) + 1
                        if strict_dup_filtering:
                            if start > 0:
                                read_starts[start] += 1
                            if end < max_tx_index:
                                read_ends[end] += 1
                        
    if sampling_rate < 1:
        subset = random.sample(frags, round(len(frags)*sampling_rate))
        for start, end in subset:
            for i in range(start, end+1):
                coverage[i] += 1
        return coverage
    elif raw_frags:
        return frags
    else:
        return coverage


def calculate_coverage(loci, samples, antisense = False, frag_sizes = None, sampling_rate = 1, paired = True, rmdup_strictest = False, dup_limit = 0, strict_dup_filtering = False):
    """
    For each sample, calculates coverage spanning each read pair
    (including the insert). Insert is assumed to be the region
    from the end of read 1 to start of read 2, ignoring any
    potential intron.

    For bisulfite, coverage is calculated separately for each conversion event.

    Coverage is stored in the profiles dict. profiles[sampleID][featID][position] = count_list
    Totals are totals[sampleID] = [+ counts, - counts]
    
    """
    profiles = {}
    total_counts = {}

    for k, (feat_id, feat_data) in enumerate(loci.items()):
        chr, strand, exons = feat_data[:3]
        start = exons[0][0]
        end = exons[-1][1]
        coord = '%s:%d-%d' % (chr, start, end)
        tx_coords = ru.tx_indexing(exons, strand == '-')

        if antisense:
            strand = ru.anti_strand_str[strand]

        for sample_label, bamfile_dict in sorted(samples.items()):
            if sample_label not in profiles:
                profiles[sample_label] = {}
                total_counts[sample_label] = [[], []] #strand separated +, -
            
            if feat_id not in profiles[sample_label]:
                profiles[sample_label][feat_id] = [[] for i in range(1+max(tx_coords.values()))]

            bamfiles = [ bfile for btuple in bamfile_dict[strand] for bfile in btuple ]
            for bamfile in bamfiles:
                sample_cov = tx_pileup([bamfile], coord, tx_coords, SAMFLAGS[strand], frag_sizes = frag_sizes, dup_limit = dup_limit, strict_dup_filtering = strict_dup_filtering, sampling_rate = sampling_rate, Fsamflag = FSAMFLAGS[strand])
                for i, x in enumerate(profiles[sample_label][feat_id]):
                    profiles[sample_label][feat_id][i].append(sample_cov[i])

        #Update total counts
        if strand == '+':
            j = 0
        else:
            j = 1
            
        for sample_label, feat_dict in profiles.items():
            for counts in feat_dict[feat_id]:
                if not counts:
                    continue
                if len(total_counts[sample_label][j]) == 0:
                    total_counts[sample_label][j] = np.array(counts)
                else:
                    total_counts[sample_label][j] = np.add(total_counts[sample_label][j], counts)

        #Status counter
        if k % 50 == 0:
            sys.stderr.write('%d out of %d regions completed\n' % (k, len(loci)))

    return profiles, total_counts



###
# FRAGMENT STATISTICS
### 

def calculate_frag_stats(loci, samples, buffer = 75, antisense = False):
    """
    Calculates the count, mean, median, sd of apparent fragment length
    per locus

    Return is a tuple of dict: frag_stats[sampleID][featID] = (N, mean, median, sd)

    and combined_stats[sampleID] = (N, mean, median, sd)
    """

    frag_lengths = {} #per feature per sample
    combined_lengths = {} #pooled over all features per sample

    frag_stats = {}
    combined_stats = {}
    combined_hist = {}

    frag_counts = {} #per feature per sample
    
    for k, (feat_id, feat_data) in enumerate(loci.items()):
        chr, strand, exons = feat_data[:3]
        start = exons[0][0]
        end = exons[-1][1]
        coord = '%s:%d-%d' % (chr, start, end)
        tx_coords = ru.tx_indexing(exons, strand == '-')
        tx_coords_inv = ru.tx_indexing(exons, strand == '-', inverse=True)    

        feat_len = max(tx_coords_inv.keys())
        if antisense:
            strand = ru.anti_strand_str[strand]
            
        for sample_label, bamfile_dict in sorted(samples.items()):
            if sample_label not in frag_lengths:
                frag_lengths[sample_label] = {}
                combined_lengths[sample_label] = []
                frag_stats[sample_label] = {}
                combined_stats[sample_label] = []
                combined_hist[sample_label] = {}
                frag_counts[sample_label] = {}
                
            if feat_id not in frag_lengths[sample_label]:
                frag_lengths[sample_label][feat_id] = []
                frag_counts[sample_label][feat_id] = {}
                
            bamfiles = [ bfile for btuple in bamfile_dict[strand] for bfile in btuple ]
            length_list = []
            samflags = SAMFLAGS[strand]
            Fsamflag = FSAMFLAGS[strand]
            
            gc = tx_pileup(bamfiles, coord, tx_coords, samflags, raw_frags = True, Fsamflag = Fsamflag)
            for s, e in gc:
                if s > buffer and e < feat_len - buffer:
                    length_list.append(e-s)
                frag_counts[sample_label][feat_id][(s,e)] = frag_counts[sample_label][feat_id].get((s,e), 0) + 1
                    
            #update
            frag_lengths[sample_label][feat_id] += length_list
            combined_lengths[sample_label] += length_list
            
        #Statistics per feature
        for sample_label in frag_lengths.keys():
            a = frag_lengths[sample_label][feat_id]
            if len(a) > 0:
                frag_stats[sample_label][feat_id] = (len(a), np.mean(a), np.median(a), np.std(a))
            else:
                frag_stats[sample_label][feat_id] = (0, float('nan'), float('nan'), float('nan'))

            #Number of times each fragment appears
            a = list(frag_counts[sample_label][feat_id].values())
            frag_counts[sample_label][feat_id] = np.bincount(a)

                
    #Statistics per sample
    for sample_label, numbers in combined_lengths.items():
        if len(numbers) > 0:
            combined_stats[sample_label] = (len(numbers), np.mean(numbers), np.median(numbers), np.std(numbers))
        else:
            combined_stats[sample_label] = (0, float('nan'), float('nan'), float('nan'))
        combined_hist[sample_label] = np.bincount(numbers)
        
    return frag_stats, combined_stats, combined_hist, frag_counts


###
# Bisulfite conversion counting
###

def index_bases(reads, ref_chr_seq = None, keep_insertions = False):
    """
    For a single read pair,
    
    Return: {genomic_coord: (read_base, ref_base)}
    Return: [(read_seq, genomic_seq, genomic_positions), ...]
    
    """
    
    mapped_reads = []
    for start, cigar, read_seq in reads:
        read_segments = []
        genomic_segments = []
        genomic_positions = []
        read_pointer = 0
        genomic_pointer = start
        
        cigar_ops = ru.decode_cigar(cigar)
        
        for nts, op in cigar_ops:
            nts = int(nts)

            if op == 'M' or op == 'X' or op == '=':
                if ref_chr_seq:
                    genomic_segments.append(ref_chr_seq[genomic_pointer:(genomic_pointer + nts)])
                genomic_positions += list(range(genomic_pointer, genomic_pointer + nts))
                genomic_pointer += nts
                read_segments.append(read_seq[read_pointer:(read_pointer + nts)])
                read_pointer += nts

            elif op == 'I':
                if keep_insertions:
                    read_segments.append(read_seq[read_pointer:(read_pointer + nts)])
                    genomic_segments.append(''.join(['-'] * nts))
                    genomic_positions += [genomic_pointer - 0.5] * nts
                    
                read_pointer += nts

            elif op == 'D' or op == 'N' or op == 'S' or op == 'P':
                genomic_pointer += nts
        
        read_seq = ''.join(read_segments)
        genomic_seq = ''.join(genomic_segments)

        mapped_reads.append((read_seq, genomic_seq, genomic_positions))

    return mapped_reads


def conversion_blocks(wt_base, bamfile_list, coord, strand, ref_chr_seq, samflags = None, read_select_fn = None, k=1, Fsamflag = 0):
    """
    base_counts[(genomic_coord)] = (wt_count, mut_count, other_count)
    """

    base_counts = {}

    if strand == '-':
        wt_base = ru.rc(wt_base)
    mut_base = ru.nt_mutations[wt_base]
    
    for bamfile in bamfile_list:
        read_hash = bam_read_pairs(bamfile, coord, samflags, Fsamflag = Fsamflag)

        for key, reads in read_hash.items():
            if len(reads) <= 2:
                mapped_reads = index_bases(reads, ref_chr_seq)
                if read_select_fn:
                    filtered_reads = []                    
                    for read in mapped_reads:
                        if read_select_fn(read, strand):
                            filtered_reads.append(read)
                    mapped_reads = filtered_reads

                bases = {}
                for read_seq, genomic_seq, genomic_positions in mapped_reads:
                    bases.update({x[0]: (x[1], x[2]) for x in zip(genomic_positions, read_seq, genomic_seq)})
                
                filtered_bases = [(pos, base) for pos, (base, ref_base) in bases.items() if ref_base == wt_base]
                if not filtered_bases:
                    continue
                
                filtered_bases.sort()
                if strand == '-':
                    filtered_bases.reverse()
                poss, ref_bases = zip(*filtered_bases)
                
                for i in range(len(poss) - k):
                    if poss[i] not in base_counts:
                        base_counts[poss[i]] = [0,0,0]
                        
                    b = ref_bases[i:i+k]
                    if b.count(wt_base) == k:
                        j = 0
                    elif b.count(mut_base) > 0:
                        j = 1
                    else:
                        j = 2
                    base_counts[poss[i]][j] += 1
                    
    return base_counts
    

def calculate_conversions(loci, samples, k = 1, antisense = False, genome_fa = None, read_select_fn = None):
    """
    conv[sampleID][featID][base][position] = [(wt, mut, other), ...]
    """

    try:
        chr_seqs = ru.load_chr_seqs(genome_fa)
    except IOError:
        sys.exit('***ERROR: Genome fasta not found: %s, exiting' % genome_fa)    
        
    profiles = {}
    bases = ['C', 'G']
    
    for counter, (feat_id, feat_data) in enumerate(loci.items()):
        chr, strand, exons = feat_data[:3]
        start = exons[0][0]
        end = exons[-1][1]
        coord = '%s:%d-%d' % (chr, start, end)
        tx_coords = ru.tx_indexing(exons, strand == '-')        
        tx_coords_inv = ru.tx_indexing(exons, strand == '-', inverse=True)    

        ref_chr_seq = chr_seqs[chr]
        
        if antisense:
            strand = ru.anti_strand_str[strand]

        for sample_label, bamfile_dict in sorted(samples.items()):
            if sample_label not in profiles:
                profiles[sample_label] = {}
                
            if feat_id not in profiles[sample_label]:
                profiles[sample_label][feat_id] = {}
                profiles[sample_label][feat_id]['C'] = [[] for i in range(1+max(tx_coords.values()))]
                profiles[sample_label][feat_id]['G'] = [[] for i in range(1+max(tx_coords.values()))]
                                
            bamfiles = bamfile_dict[strand]
            for btuple in bamfiles:            
                for bamfile, base in zip(btuple, bases):
                    conv = conversion_blocks(base, [bamfile], coord, strand, ref_chr_seq, SAMFLAGS[strand], read_select_fn, k = k, Fsamflag = FSAMFLAGS[strand])

                    for i, x in enumerate(profiles[sample_label][feat_id][base]):
                        genomic_coord = tx_coords_inv.get(i, -1)
                        if genomic_coord > -1:
                            profiles[sample_label][feat_id][base][i].append(conv.get(genomic_coord, [0,0,0]))

        #Status counter
        if counter % 50 == 0:
            sys.stderr.write('%d out of %d regions completed\n' % (counter, len(loci)))
    return profiles



if __name__ == "__main__":
    """
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('task', type=str, help='task to run', choices = ['coverage', 'mut_coverage', 'conversions', 'stats'])
    parser.add_argument('-t', '--test', action='store_true', help='Do a test run with only one UTR (default=False)')    
    parser.add_argument('-o', '--outprefix', type=str, help='Output prefix for pickle file (default=./resa)', default='./resa')
    parser.add_argument('--singleend', action='store_true', help='Sequencing was performed single-end (treating paired-end data as single end is not recommended) (default=False)')    
    parser.add_argument('--read1strand', type=str, default='+', choices = ['+', '-'], help='Strand of read 1 relative to the transcript (default=+)')    
    parser.add_argument('--primaryonly', type=int, help='Use only primary alignments (SAM flags of 0 and 16 for single-end reads, 99/147 and 83/163 for paired-end reads) (default=True)', default=1)
    parser.add_argument('--antisense', action='store_true', help='Operate on antisense strand relative to gene (default=False)')    
    parser.add_argument('--minfraglen', type=int, help='Minimum end-to-end fragment length of read pairs to include (only compatible with coverage task) (default=0)', default=0)
    parser.add_argument('--maxfraglen', type=int, help='Maximum end-to-end fragment length of read pairs to include (only compatible with coverage task) (default=500)', default=500)
    parser.add_argument('--duplimit', type=int, help='Maximum number of times a fragment with identical start and end coordinate can be used. Only compatible with coverage task (default=no limit)', default=0)
    parser.add_argument('--strictdupfilter', action='store_true', help='Impose a stricter policy for read filtering by (non-deterministically) selecting at most <duplimit> fragments that start or end at each coordinate. I.e., a read pair is ignored if one of the read starts has already been seen <duplimit> number of times. However, reads aligning to the ends of UTR amplicons are exempted. Only compatible with coverage task, and --duplimit must also be specified (default=False)')    
    parser.add_argument('--linkedbases', type=int, help='For conversion task, number of consecutive bases to group together (default=1)', default=1)
    
    parser.add_argument('--samples_as_bam', action="store_true", help='Treat sample_list as bamfiles rather than directory names (only valid with coverage and stats tasks) (default=False)')
    req = parser.add_argument_group('required named arguments')
    req.add_argument('-u', '--utrbed', type=str, help='Utr BED12 file', required=True)
    req.add_argument('-f', '--fasta', type=str, help='Utr FASTA file', required = True)    
    req.add_argument('-l', '--group_labels', type=str, help='Group labels (comma separated, e.g. Early,Late,Treated)', required=True)
    req.add_argument('-g', '--genome_fasta', type=str, help='Genome fasta file, only required for conversions task')    
    parser.add_argument('sample_list', type=str, nargs='+', help='Comma-delimited list of condition paths (i.e., replicates) as created by resa_map, one list for each of the group labels')
    
    args = parser.parse_args()
    sys.stderr.write('\nResa_coverage beginning %s\n' % str(datetime.now()))

    if args.task in ['mut_coverage', 'conversions']:
        mut = True
    else:
        mut = False
    
    #Process the UTR bed and fasta files
    loci = ru.initialize_loci(args.utrbed, args.fasta, test=args.test)

    #Organize the samples
    if args.samples_as_bam and args.task not in ['coverage', 'stats']:
        sys.exit('***ERROR: --samples_as_bam incompatible with task %s, exiting.' % args.task)

    samples = process_sample_lists(args.sample_list, args.group_labels, mut = mut, samples_as_bamfiles = args.samples_as_bam)

    #Set SAM flags
    if args.singleend:
        SAMFLAGS = SAMFLAGS_SINGLE
        if args.primaryonly:
            FSAMFLAGS = FSAMFLAGS_STRICT_SINGLE
        else:
            FSAMFLAGS = FSAMFLAGS_LOOSE_SINGLE
    else:
        if not args.primaryonly:
            FSAMFLAGS = FSAMFLAGS_LOOSE_PAIRED
            
    #Interpret strand
    if (args.read1strand == '-' and not args.antisense) or (args.read1strand == '+' and args.antisense):
        flip_strand = True
    else:
        flip_strand = False
        
    #Do tasks
    sys.stderr.write('Performing %s on %d regions with %d samples\n' % (args.task, len(loci), len(samples)))
    if args.task == 'coverage' or args.task == 'mut_coverage':
        result = calculate_coverage(loci, samples, frag_sizes = (args.minfraglen, args.maxfraglen), antisense = flip_strand, dup_limit = args.duplimit, strict_dup_filtering = args.strictdupfilter)
    elif args.task == 'stats':
        result = calculate_frag_stats(loci, samples, antisense = flip_strand)
    elif args.task == 'conversions':
        if not args.genome_fasta:
            sys.exit('***ERROR: --genome_fasta required for conversions task, exiting.')
        result = calculate_conversions(loci, samples, antisense = flip_strand, k = args.linkedbases, genome_fa = args.genome_fasta)
        args.task = 'conversions_%d' % args.linkedbases

        
    #Pickle the output
    file_prefix = '%s.%s' % (args.outprefix, args.task)
    if args.antisense:
        file_prefix += '.antisense'
    outfile = file_prefix + '.pkl'
    sys.stderr.write('Writing results to %s\n' % outfile)
    output = open(outfile, 'wb')
    pickle.dump(result, output, pickle.HIGHEST_PROTOCOL)
    output.close()
    sys.stderr.write('Done. %s\n\n' % str(datetime.now()))

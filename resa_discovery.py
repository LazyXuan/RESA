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
# resa_discovery.py: RESA feature identification
##

import resa_util as ru

import argparse
import math
import operator
import pickle
import re
import sys

import numpy as np
import scipy.stats as sps
import statsmodels.stats.multitest as smm


###
# Positional counting
###

def raw_profile(feat_id, profiles, cond, ngroups = 1, skip = 0):
    """
    Returns the sum per position of raw read counts
    of the pooled samples for the condition

    If ngroups = 0, then the return is trivially a subset of
    the raw_profile

    If ngroups is -1, each pair of raw counts is summed and returned
    as a tuple
    
    If ngroups is -2, every other raw count in returned
    as a tuple. skip = 0 means evens, skip = 1 means odds
    
    If ngroups = 2 and skip = 0, sum is over even samples
    only. If ngroups = 2 and skip = 1, odd samples.
    """
    
    counts, totals = profiles

    if ngroups == 0:
        return [tuple(x) for x in counts[cond][feat_id]]
    elif ngroups == -1:
        return [tuple(x[skip::2]) for x in counts[cond][feat_id]]
    elif ngroups == -1:
        return np.array([tuple([y+z for y,z in zip(x[0::2], x[1::2])]) for x in counts[cond][feat_id]])
    else:
        return np.array([sum(x[skip::ngroups]) for x in counts[cond][feat_id]])

def raw_profile_replicates(feat_id, profiles, cond):
    return raw_profile(feat_id, profiles, cond, ngroups = 0)


def profile_ratio(loci, profiles, cond1, cond2, smooth = 0.5, groups = 1, positive_only = False, negative_only = False):
    """
    Log2 ratio of cond1 counts over cond2 counts
    Profiles is a tuple of (counts, totals)
    
    Default formulation normalizes by total counts in the experiment.
    Alternative formulation normalizes counts per feature individually
    by the total number of counts in that feature

    Return is ratios[feat_id] = value
    
    If groups is 2, then two different sets of ratios are
    returned, using only the odd or even samples for each
    condition (e.g., for C vs G conversion)

    Return is (ratios1, ratios2)
    """
    
    counts, totals = profiles
    
    a = sum(totals[cond1][0]) + sum(totals[cond1][1])
    b = sum(totals[cond2][0]) + sum(totals[cond2][1])
    norm_factor = np.log2(a) - np.log2(b)

    smooth1 = smooth
    smooth2 = smooth * (b/a)
        
    ratios = [{} for i in range(groups)]
    
    for feat_id, feat_data in loci.items():
        for i in range(groups):
            sum_counts1 = np.array([sum(x[i::groups])+smooth1 for x in counts[cond1][feat_id]])
            sum_counts2 = np.array([sum(x[i::groups])+smooth2 for x in counts[cond2][feat_id]])

            r = np.log2(sum_counts1) - np.log2(sum_counts2) - norm_factor

            if positive_only:
                r = np.where(r > 0, r, 0)
            elif negative_only:
                r = np.where(r < 0, r, 0)
                
            ratios[i][feat_id] = r

    if len(ratios) == 1:
        return ratios[0]
    else:
        return ratios

def clip_positive_profile_ratio(loci, profiles, cond1, cond2, smooth = 0.5):
    return profile_ratio(loci, profiles, cond1, cond2, smooth = smooth, positive_only = True)

###
#
###

def bedgraph_entries(chr, numbers, tx_coords):
    """
    Given an array of numbers, convert to bedgraph entries
    indexed by genomic coordinates. Return is a list of tuples
    (chr, start, end, number)
    """

    return [(chr, tx_coords[i], tx_coords[i]+1, x) for i, x in enumerate(numbers)]


###
# Motif location identification
###

def motif_locs(seq, motifs, window = 1):
    """
    Given a sequence, returns the locations of
    motifs that are separated by window number of nts
    start-to-start
    """
    positions = []

    #identify motif positions
    for label, motif in motifs:
        for m in re.finditer(motif, seq):
            positions.append((m.start(), m.end(), motif, label))

    #Filter the ranges to remove any overlapping ones
    positions_filtered = []
    for s, e, motif, label in positions:
        counter = 0
        for p2 in positions:
            if abs(s-p2[0]) < window:
                counter += 1
        if counter == 1 and (label,motif) in motifs:
            positions_filtered.append((s, e, label, motif))
            
    return positions_filtered


###
# Responsive region identification
###

def responsive_regions(numbers, threshold = math.log(1.25,2), search_window = 100, boundary_percentile = 0.25, summit_percentile = 0.9, enrich_percentile = 0.5):
    """
    Identify regions that have threshold fold or greater deviation from 0
    (peaks, valleys)

    Strategy: identify a local min or max in a window. Extend the window
    outward.
    """
    
    #Identify critical points
    critical_pts = []    
    for fn_extrema, fn_op, thresh in (np.argmin,operator.lt,-1*threshold), (np.argmax, operator.gt, threshold):
        i = 0        
        while True:
            right_boundary = min(len(numbers), i+search_window)
            if right_boundary - i < 5:
                break
            critical_pt = fn_extrema(numbers[i:right_boundary])
            if 0 <= critical_pt < (right_boundary - i)-1 or critical_pt + i == len(numbers)-1:
                value = numbers[i+critical_pt]
                if fn_op(value, thresh):
                    critical_pts.append((i+critical_pt, value))
                i = right_boundary
            else:
                i = i+search_window//2

    #For each critical point, extend left and right
    regions = []
    for cp, value in critical_pts:
        r_data = {'value': value, 'crit_pt': cp}
        #establish the sign
        if value > 0:
            cmp_fn = operator.lt
            r_data['label'] = 'P%d' % cp
            r_data['type'] = 'peak'
        else:
            cmp_fn = operator.gt
            r_data['label'] = 'V%d' % cp
            r_data['type'] = 'valley'
            
        #Calculate region boundaries
        for i in range(cp, -1, -1):
            if cmp_fn(numbers[i], boundary_percentile * value):
                break
        left_boundary = i + 1
        
        for i in range(cp, len(numbers)):
            if cmp_fn(numbers[i], boundary_percentile * value):
                break
        right_boundary = i# - 1

        r_data['bound'] = (left_boundary, right_boundary)

        if right_boundary - left_boundary > 1:
            regions.append((left_boundary, right_boundary, r_data))
        
    #Merge regions if necessary, and re-calculate boundaries
    regions.sort(key=lambda x: x[0:2])
    kept = []
    for left, right, region in regions:
        if not kept or left >= kept[-1]['bound'][1]:
            kept.append(region)
        else:
            #overlapping regions are merged
            far_left = min(left, kept[-1]['bound'][0])
            far_right = max(right, kept[-1]['bound'][1])
            
            if kept[-1]['type'] == region['type']:
                if abs(region['value']) > abs(kept[-1]['value']):
                    kept[-1] = region
                kept[-1]['bound'] = (far_left, far_right)

    #Calculate summit boundaries and enriched boundaries, refine region boundary
    for r_data in kept:
        for percentile, zone_label in zip((summit_percentile, enrich_percentile, boundary_percentile), ('summit', 'enriched', 'boundary')):
            summits = []
            for i in range(r_data['bound'][0], r_data['bound'][1]):
                if abs(numbers[i]) > percentile * abs(r_data['value']):
                    if not summits or i - summits[-1][-1] > 3:
                        summits.append([i])
                    else:
                        summits[-1].append(i)
            summits = [ (x[0], x[-1]+1) for x in summits]
            r_data[zone_label] = summits

        r_data['bound'] = (r_data['boundary'][0][0], r_data['boundary'][-1][1])
        r_data['width'] = r_data['bound'][1] - r_data['bound'][0]

    return kept


def output_region_data(regions, outfile):
    fw = open(outfile, 'w')

    header = ['feat_id', 'label', 'coord', 'expt', 'ctrl', 'type', 'value', 'width', 'crit_pt', 'motifs', 'sum_expt_counts', 'sum_ctrl_counts']
        
    fw.write('\t'.join(header) + '\n')
    out_lines = []
    for region_list in regions.values():
        for r in region_list:
            out_line = '\t'.join([ru.pretty_str(r.get(x, '.')) for x in header])
            out_lines.append((r['value'], out_line))
    out_lines.sort()
    for value, out_line in out_lines:
        fw.write(out_line + '\n')
    fw.close()

    
def output_region_bed(regions, outfile):
    fw = open(outfile, 'w')
    fields = ['chr', 'genomic_start', 'genomic_end', 'long_label', 'value', 'strand']
    
    for region_list in regions.values():
        for r in region_list:
            out_line = '\t'.join([ru.pretty_str(r.get(x, '.')) for x in fields])
            fw.write(out_line + '\n')
    fw.close()

    
def output_region_seq(regions, outfile_root, do_valley = False):
    if do_valley:
        fv = open(outfile_root+'.valley.fa', 'w')
    fp = open(outfile_root+'.peak.fa', 'w')
    for region_list in regions.values():
        for r in region_list:
            if r['type'] == 'valley':
                if do_valley:
                    fw = fv
                else: continue
            else:
                fw = fp

            fw.write('>' + r['long_label'] + '\n')                
            fw.write(r['seq'] + '\n')
    if do_valley:
        fv.close()
    fp.close()

    
def output_region_bedgraph(ratio_list, outfile, header_name = ''):
    fw = open(outfile, 'w')
    if header_name:
        fw.write('track type=bedGraph name=%s\n' % header_name)
    
    for chr, start, end, score in sorted(ratio_list):
        fw.write('\t'.join([chr, str(start), str(end), '%.3f' % score]) + '\n')
    fw.close()
    
        
def region_finder(loci, profiles, expt, ctrl, ratio_fn = profile_ratio, count_fn = raw_profile_replicates, motifs = None, lfc_thresh = math.log(1.25, 2), low_conf_lfc_thresh = math.log(1.25, 2), min_peak_width = 5, min_ctrl_coverage = 0, min_expt_coverage = 0, outfile_root = '', peaks_only = False, do_bedgraph = True, write_bedgraph_header = True):
    """
    Iterate through ratios to identify responsive regions
    (peaks/valleys)

    Minimal return: regions[feat_id] = [(start, end, label), ...]
    """
    regions = {}
    ratio_cache = []

    try:
        pr = ratio_fn(loci, profiles, expt, ctrl)
    except:
        sys.exit('***ERROR: Improperly formatted data file? (should be a _coverage_ pickle)')
        
    for feat_id, feat_data in loci.items():
        chr, strand, exons, seq = feat_data
        tx_coords = ru.tx_indexing(exons, strand == '-', inverse = True)
        
        #Find responsive regions
        rr = responsive_regions(pr[feat_id], threshold = low_conf_lfc_thresh)

        #Update the ratio list
        if do_bedgraph:
            ratio_cache += bedgraph_entries(chr, pr[feat_id], tx_coords)
        
        #Identify motif locations
        feat_motifs = motif_locs(seq, motifs, window = 1)

        #Get raw counts
        all_expt_counts = count_fn(feat_id, profiles, expt)
        all_ctrl_counts = count_fn(feat_id, profiles, ctrl)
            
        #For each responsive region, calculate raw reads
        for r_data in rr:
            #Get raw counts at the critical point
            r_data['expt_counts'] = all_expt_counts[r_data['crit_pt']]
            r_data['ctrl_counts'] = all_ctrl_counts[r_data['crit_pt']]
                
            r_data['sum_expt_counts'] = sum(r_data['expt_counts'])
            r_data['sum_ctrl_counts'] = sum(r_data['ctrl_counts'])

            if r_data['sum_ctrl_counts'] < min_ctrl_coverage or \
               r_data['sum_expt_counts'] < min_expt_coverage or \
               r_data['width'] < min_peak_width or \
               abs(r_data['value']) < lfc_thresh or \
               (peaks_only and r_data['type'] == 'valley'):
                continue

            r_data['expt'] = expt
            r_data['ctrl'] = ctrl
            r_data['feat_id'] = feat_id
            r_data['long_label'] = r_data['feat_id'] + '_' + r_data['label'] 
            r_data['dist_from_end'] = min(r_data['bound'][0]-1, len(seq)-r_data['bound'][1]-1)
            r_data['crit_pt_dist_from_end'] = min(r_data['crit_pt'], len(seq) - r_data['crit_pt'])

            g_start, g_end = tx_coords[r_data['bound'][0]], tx_coords[r_data['bound'][1]]
            if strand == '-':
                g_start, g_end = g_end, g_start

            r_data['chr'] = chr
            r_data['genomic_start'] = g_start
            r_data['genomic_end'] = g_end + 1
            r_data['strand'] = strand
            r_data['coord'] = '%s:%d-%d' % (chr, g_start+1, g_end+1) 

            r_data['seq'] = seq[r_data['bound'][0]:(1+r_data['bound'][1])]

            #Get the sequence of the summit
            flank_amt = 10
            summit_seqs = []
            min_s = len(seq)
            max_e = 0
            for s, e in r_data['summit']:
                sseq = seq[s:e]
                l_bound = max(0,s-flank_amt)
                flank_l = seq[l_bound:s].lower()
                r_bound = min(e+flank_amt,len(seq))
                flank_r = seq[e:r_bound].lower()
                summit_seqs.append(flank_l + sseq + flank_r)
                if l_bound < min_s:
                    min_s = l_bound
                if r_bound > max_e:
                    max_e = r_bound

            r_data['summit_seq'] = summit_seqs

            #15nt +/-
            flank_addition = 15
            r_data['spanning_summit_seq'] = seq[max(0, min_s-flank_addition):min(max_e + flank_addition, len(seq))]

            #Annotate motifs
            overlapping_motifs = []                    
            for s, e, short_name, m in feat_motifs:
                if r_data['bound'][0] < s and e < r_data['bound'][1]:
                    overlapping_motifs.append((s, short_name))
            r_data['motifs'] = overlapping_motifs

            if feat_id not in regions:
                regions[feat_id] = []
            regions[feat_id].append(r_data)

    if outfile_root:
        output_region_data(regions, outfile_root + '.txt')
        output_region_bed(regions, outfile_root + '.bed')        
        output_region_seq(regions, outfile_root, do_valley = not peaks_only)

        if do_bedgraph:
            if write_bedgraph_header:
                header_name = header_name = '%s/%s' % (expt, ctrl)
            else:
                header_name = ''
            output_region_bedgraph(ratio_cache, outfile_root + '.bedgraph', header_name = header_name)
    return regions
    

##
# Bisulfite Conversion task
##

def fdr_correct(p_values, alpha=1E-4):
    """
    Returns an FDR correction for the input p-value list, as a dict
    """
    p = [ x for x in p_values if not math.isnan(x) ]
    
    sig_calls, adj_p, alpha_sidak, alpha_bonfer = smm.multipletests(p, alpha = alpha, method='fdr_i')
    return dict(zip(p, adj_p))


def g_test(obs_list, calculate_p = False):
    """
    G-test of independence. Obs_list encodes an rxc
    set of observations as an r-length list of c-tuples.

    Slow method.
    
    Returns the g statistic, or a tuple of (g, pval) if
    calculate_p is True
    """
    
    col_sums = [sum(x) for x in zip(*obs_list)]
    total_sum = sum(col_sums)

    g_stat = 0
    for row in obs_list:
        row_sum = sum(row)
        expected = [row_sum * x/total_sum for x in col_sums]

        for o, e in zip(row, expected):
            try:
                g_stat += o*math.log(o/e)
            except:
                pass
    g_stat *= 2
    
    if calculate_p:
        p_val = sps.chi2.sf(g_stat, (len(obs_list) - 1) * (len(obs_list[0]) - 1))
        return g_stat, p_val
    else:
        return g_stat

    
def repeated_g_test_indep(observation_list, bkgrd_observation_list, calculate_p = False):
    """
    Pooled G test is performed. Then a G test for heterogeneity
    is performed for each of the two observation lists.

    If calculate_p is True, return is a 3-tuple of (g_stat, p_val) tuples
    """
    o1 = sum([ x[0] for x in observation_list ])
    o2 = sum([ x[1] for x in observation_list ])
    
    bo1 = sum([ x[0] for x in bkgrd_observation_list ])
    bo2 = sum([ x[1] for x in bkgrd_observation_list ])

    g_pool = g_test([(o1,o2), (bo1,bo2)], calculate_p = calculate_p)
    g_het = g_test(observation_list, calculate_p = calculate_p)
    g_het_bkgrd = g_test(bkgrd_observation_list, calculate_p = calculate_p)

    return g_pool, g_het, g_het_bkgrd


def conversion_ratio(loci, bisulf_profiles, cond1, cond2 = None, smooth = 0.5, znorm = False):
    """
    Given a bisulfite object, calculates the ratio of
    unconverted nts / unconverted + converted nts per position for cond1
    If cond2 is specified, the same is calculated for cond2,
    then the ratio of cond1/cond2 is calculated
    
    Return is a tuple (C_ratios, G_ratios). C_ratios[feat_id][pos] = ratio

    WARNING: z normalization relies on having a full distribution
    to calculate mean/std.
    """
    bp = bisulf_profiles
    ratios = [{}, {}]
    population = [np.array([]), np.array([])]
    
    for feat_id, feat_data in loci.items():
        seq = np.array(list(feat_data[3]))
        for i, base in enumerate(['C', 'G']):
            rtemp = []            
            for cond in cond1, cond2:
                if cond:
                    sum_unconv = np.array([sum(list(zip(*x))[0])+smooth for x in bp[cond][feat_id][base]])
                    sum_conv = np.array([sum(list(zip(*x))[1]) for x in bp[cond][feat_id][base]])
                    sum_total = sum_unconv + sum_conv
            
                    rtemp.append(np.log2(sum_unconv) - np.log2(sum_total))
            if len(rtemp) == 2:
                r = rtemp[0] - rtemp[1] #ratio of ratios
            else:
                r = rtemp[0]
                    
            r = np.where(seq == base, r, float('nan')) #only keep informative positions
            
            if znorm:
                population[i] = np.append(population[i], r)
                #m = np.nanmean(r)
                #s = np.nanstd(r)
                #r = (r - m) / s
            
            ratios[i][feat_id] = r
    if znorm:
        for ratio_list, popul in zip(ratios, population):
            m = np.nanmean(popul)
            s = np.nanstd(popul)
            
            for feat_id, r in ratio_list.items():
                ratio_list[feat_id] = (r - m) / s
            
    return ratios


def conversion_gtest(loci, bisulf_profiles, cond1, cond2, smooth = 0.5, padj = True):
    """
    Return is a tuple (C_gstats, G_gstats). C_gstats[feat_id][pos] = ((Gpool, pval, padj), (Ghet, pval, padj),
    (Ghet, pval, padj))

    P adjustment is performed using FDR for each test across all positions (3 total)
    """

    def g_fn(do, c1_counts, c2_counts):
        if not do:
            return (())
        c1 = [(x[0]+smooth, x[1]+smooth) for x in c1_counts]
        c2 = [(x[0]+smooth, x[1]+smooth) for x in c2_counts]
        return repeated_g_test_indep(c1, c2, True)

    bp = bisulf_profiles
    gstats = [{}, {}]
    
    for feat_id, feat_data in loci.items():
        seq = np.array(list(feat_data[3]))
        for i, base in enumerate(['C', 'G']):
            pos = np.where(seq == base, True, False)
            g = [g_fn(*x) for x in zip(pos, bp[cond1][feat_id][base], bp[cond2][feat_id][base])]
            gstats[i][feat_id] = g

    if padj:        
        padj_dicts = []
        for i in range(3):
            pvals = [x[i][1] for z in gstats for y in z.values() for x in y if x]
            padj_dicts.append(fdr_correct(pvals, alpha=1E-4))
            
        def g_padj(g_entry):
            if not g_entry:
                return g_entry
            return tuple([(g[0], g[1], padj_dicts[i][g[1]]) for i, g in enumerate(g_entry)])
            
        for g in gstats:
            for feat_id, values in g.items():
                g[feat_id] = [g_padj(x) for x in values]

    return gstats


def conversion_biases(g_stats, z_scores, ratios, g_fdr = 0.05, z_thresh = 1, find_depletions = True):
    """
    Finds positions whose wt/total base ratio is less than (default) or greater than
    (if depleted = False) z_thresh. Then further filters based on the result of the
    G statistic
    """
    if find_depletions:
#        op = operator.lt
        op = operator.le
        z_thresh = -1 * abs(z_thresh)
    else:
#        op = operator.gt
        op = operator.ge        
        z_thresh = abs(z_thresh)
    
    p_data_list = []
    for j, (g, z, r) in enumerate(zip(g_stats, z_scores, ratios)):
        if g and g[0][2] <= g_fdr and op(z, z_thresh):
            p_data = {}
            p_data['pos'] = j
            p_data['g_pool'], p_data['g_pool_p'], p_data['g_pool_padj'] = g[0]
            p_data['g_het_expt'], p_data['g_het_expt_p'], p_data['g_het_expt_padj'] = g[1]
            p_data['g_het_ctrl'], p_data['g_het_ctrl_p'], p_data['g_het_ctrl_padj'] = g[2]
            p_data['z_expt'] = z
            p_data['log_diff'] = r
            p_data_list.append(p_data)
    return p_data_list


def base_lookahead(seq, base, pos, mult = 1):
    poss = [pos]
    for i in range(mult-1):
        poss.append(seq.find(base, poss[-1]+1))
    return tuple(poss)
    
            
def conversions(loci, profiles, expt, ctrl, mult = 1, z_thresh = 1, g_fdr = 0.05, find_depletions = True):
    """
    Iterate through positions, find positions with
    significant asymmetry in expt vs ctrl.
    """

    conv = conversion_gtest(loci, profiles, expt, ctrl, padj = True)
    convz = conversion_ratio(loci, profiles, expt, znorm = True)
    convr = conversion_ratio(loci, profiles, expt, ctrl)
    
    positions = {}
    for feat_id, feat_data in loci.items():
        seq = feat_data[3]
        lc_seq = list(seq.lower())
        positions[feat_id] = []
        for i, base in enumerate(['C', 'G']):
            p_data_list = conversion_biases(conv[i][feat_id], convz[i][feat_id], convr[i][feat_id], g_fdr, z_thresh, find_depletions = find_depletions)
            for p_data in p_data_list:
                p_data['feat_id'] = feat_id
                p_data['expt'] = expt
                p_data['ctrl'] = ctrl
                p_data['base'] = base
                p_data['linked_bases'] = mult
                p_data['pos'] = base_lookahead(seq, base, p_data['pos'], mult)
                p_data['label'] = '%s%s' % (base, ru.pretty_str(p_data['pos']))
                
                #Annotated sequence
                anno_seq = lc_seq[:]
                for index in p_data['pos']:
                    anno_seq[index] = anno_seq[index].upper()
                left = max(0, p_data['pos'][0] - 5)
                right = min(len(seq), p_data['pos'][-1]+6)
                p_data['seq'] = ''.join(anno_seq[left:right])
                
            positions[feat_id] += p_data_list

    return positions


def joint_conversions(loci, profiles, expt, ctrl, outfile = None, z_thresh = 1, g_fdr = 0.05, find_depletions = True):
    """
    E.g., mono, di, tri base alleles simultaneously and non-redundantly.
    """
    if outfile:
        fw = open(outfile, 'w')
    
    convs = []
    try:
        for i, profile in enumerate(profiles):
            convs.append(conversions(loci, profile, expt, ctrl, z_thresh = z_thresh, g_fdr = g_fdr, mult = i+1, find_depletions = find_depletions))
    except:
        sys.exit('***ERROR: Improperly formatted data file? (should be a _conversions_ pickle)')

            
#    fields = ['feat_id', 'label', 'coord', 'strand', 'base', 'mult', 'pos', 'z_expt', 'log_diff', 'g_pool', 'g_pool_p', 'g_pool_padj', 'g_het_expt', 'g_het_expt_p', 'g_het_expt_padj', 'g_het_ctrl', 'g_het_ctrl_p', 'g_het_ctrl_padj', 'seq']
    fields = ['feat_id', 'label', 'coord', 'strand', 'base', 'linked_bases', 'pos', 'z_expt', 'log_diff', 'g_pool', 'g_pool_p', 'g_pool_padj', 'seq']
    if outfile:
        fw.write('\t'.join(fields) + '\n')

    joint_convs = {}
    for feat_id, feat_data in loci.items():
        chr, strand, exons, seq = feat_data        
        tx_coords = ru.tx_indexing(exons, strand == '-', inverse = True)
        
        positions = {}
        for conv_dict in convs:
            for conv in conv_dict[feat_id]:
                for index in conv['pos']:
                    if index not in positions:
                        positions[index] = []
                    positions[index].append((conv['g_pool_padj'], conv))
        filtered_positions = {}
        for index, conv_list in positions.items():
            #Keep the conversion with lowest p value
            conv_list.sort(key = lambda x: x[0])
            best = conv_list[0][1]
            
            if best['label'] not in filtered_positions:
                filtered_positions[best['label']] = (best['pos'], best)
            for pval, conv in conv_list[1:]:
                filtered_positions[conv['label']] = ((-1,), None)

        joint_convs[feat_id] = []
        for index, conv in sorted(filtered_positions.values()):
            if conv:
                #Genomic coordinate
                g_start = tx_coords[conv['pos'][0]]
                g_end = tx_coords[conv['pos'][-1]]
                if g_start > g_end:
                    g_start, g_end = g_end, g_start
                conv['coord'] = '%s:%d-%d' % (chr, g_start+1, g_end+1)
                conv['strand'] = strand

                joint_convs[feat_id].append(conv)
                
                if outfile:
                    out_line = [ru.pretty_str(conv[x]) for x in fields]
                    fw.write('\t'.join(out_line) + '\n')

    if outfile:
        fw.close()

    return joint_convs


def conversion_finder(loci, profiles, expt, ctrl, outfile_root = '', z_thresh = 1, g_fdr = 0.05):
    
    convs = joint_conversions(loci, profiles, expt, ctrl, outfile = outfile_root + '.depletion_biases.txt', z_thresh = z_thresh, g_fdr = g_fdr)
    convs = joint_conversions(loci, profiles, expt, ctrl, outfile = outfile_root + '.enrichment_biases.txt', find_depletions = False, z_thresh = z_thresh, g_fdr = g_fdr)


if __name__ == "__main__":
    """
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('task', type=str, help='task to run', choices = ['regions', 'bisulfite', 'clip'])
    parser.add_argument('-o', '--outprefix', type=str, help='Output prefix for outfiles (default=./resa)', default='./resa')
    
    reg = parser.add_argument_group('options for Regions or Clip task')
    reg.add_argument('--min_region_width', type=int, help='Minimum width of a responsive region to return (default=65)', default=65)
    reg.add_argument('--min_ctrl_coverage', type=int, help='Minimum X coverage of a responsive region to return, as measured by raw read counts in the ctrl condition overlapping the critical point (default=0)', default=0)
    reg.add_argument('--min_expt_coverage', type=int, help='Minimum X coverage of a responsive region to return, as measured by raw read counts in the expt condition overlapping the critical point (default=0)', default=0)
    reg.add_argument('--fold_thresh', type=float, help='Fold difference threshold, expt vs ctrl, must be > 1 (default=1.25)', default=1.25)
    reg.add_argument('--motifs', type=str, help='Quoted comma-delimited list of labels:motifs (e.g., "mir430:GCACTT,are:ATTTTA") to identify in responsive regions. Python-style regular expressions are supported.', default='')
    reg.add_argument('--no_bedgraph', action='store_true', help='Do not write the bedgraph of the expt/ctrl ratio (default=False)')
    reg.add_argument('--no_bedgraph_header', action='store_true', help='Do not include a track header line on the bedgraph (default=False)')

    reg = parser.add_argument_group('options for Bisulfite task')
    reg.add_argument('--z_thresh', type=float, help='Z score threshold for base asymmetry in the expt condition, must be >= 0 (default=1)', default=1)
    reg.add_argument('--p_thresh', type=float, help='FDR adjusted p-value threshold for biased base conversion, expt vs ctrl, must be (0,1] (default=0.05)', default=0.05)
    
    req = parser.add_argument_group('required named arguments')
    req.add_argument('-u', '--utrbed', type=str, help='Utr BED12 file', required = True)
    req.add_argument('-f', '--fasta', type=str, help='Utr FASTA file', required = True)    
    req.add_argument('-d', '--data', type=str, help='Pickle file produced by resa_coverage.py appropriate for the task', required = True)
    req.add_argument('-e', '--expt', type=str, help='Name of experimental group (consistent with name used for resa_coverage.py)', required = True)
    req.add_argument('-c', '--ctrl', type=str, help='Name of control group (consistent with name used for resa_coverage.py)', required = True)
    
    args = parser.parse_args()

    #Check parameters
    if args.fold_thresh <= 1:
        sys.exit('***ERROR: --fold_thresh must be > 1, exiting.')
    if args.z_thresh < 0:
        sys.exit('***ERROR: --z_thresh must be >= 0, exiting.')
    if not (0 < args.p_thresh <= 1):
        sys.exit('***ERROR: --p_thresh must > 0 and <= 1, exiting.')

    #Process the UTR bed and fasta files
    loci = ru.initialize_loci(args.utrbed, args.fasta)

    #Open the pickle file
    profiles = pickle.load(open(args.data, 'rb'))

    #Process the motif list
    if args.motifs:
        motifs = [ tuple(m.split(':')) for m in args.motifs.split(',') ]
    else:
        motifs = []

    #Do tasks
    low_conf_thresh = min(1.25, args.fold_thresh)    
    if args.task == 'regions':
        result = region_finder(loci, profiles, args.expt, args.ctrl, motifs = motifs, lfc_thresh = math.log(args.fold_thresh, 2), low_conf_lfc_thresh = math.log(low_conf_thresh, 2), min_peak_width = args.min_region_width, min_ctrl_coverage = args.min_ctrl_coverage, outfile_root = args.outprefix + '.regions', do_bedgraph = not args.no_bedgraph, write_bedgraph_header = not args.no_bedgraph_header)
    elif args.task == 'clip':
        result = region_finder(loci, profiles, args.expt, args.ctrl, ratio_fn = clip_positive_profile_ratio, motifs = motifs, lfc_thresh = math.log(args.fold_thresh, 2), low_conf_lfc_thresh = math.log(low_conf_thresh, 2), min_peak_width = args.min_region_width, min_ctrl_coverage = args.min_ctrl_coverage, min_expt_coverage = args.min_expt_coverage, outfile_root = args.outprefix + '.clip', peaks_only = True, do_bedgraph = not args.no_bedgraph, write_bedgraph_header = not args.no_bedgraph_header)
    elif args.task == 'bisulfite':
        result = conversion_finder(loci, profiles, args.expt, args.ctrl, outfile_root = args.outprefix + '.conversions', z_thresh = args.z_thresh, g_fdr = args.p_thresh)

        

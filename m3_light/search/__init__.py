"""
Search functions for motifs
"""

import os
import m3_light
from os.path import join as pjoin

def expand_motif(motif):
    """
    Expands motif with beginning and ending nucleotide R,Y,S,W to a list of 4 motifs with alphabet A,T,C,G

    Example: RAAY -> ['AAAC', 'AAAT', 'GAAC', 'GAAT']
    """
    start_n = motif[0]
    end_n = motif[-1]
    motifs = []
    for x in m3_light.genomes.code[start_n]:
        for y in m3_light.genomes.code[end_n]:
            motifs.append(x+motif[1:-1]+y)
    return motifs
       
def all_indices(sequence, sub_sequence, offset=0):
    """
    Returns list of all indices at which sub_sequence is found in sequence
    """
    indices = []
    i = sequence.find(sub_sequence, offset)
    while i >= 0:
        indices.append(i)
        i = sequence.find(sub_sequence, i + 1)
    return indices

def overlap((start_1, stop_1), (start_2, stop_2)):
    """
    If the two given intervals overlap, return overlapping region else return first interval
    """
    if start_2>=start_1+1 and start_2<=stop_1+1:
        return ((start_1, stop_2), True)
    else:
        return ((start_1, stop_1), False)

def find_motif(sequence, motif):
    """
    :param sequence: DNA/RNA sequence
    :param motif: motif to search for with alphabet ATCGRYSW
    
    Returns all intervals at which **motif** is found in **sequence**
    """
    pool = []
    
    if sum([motif.find("R"), motif.find("Y"), motif.find("S"), motif.find("W")]) == -4:
        pool = []
        indices = all_indices(sequence, motif)
        for i in indices:
            pool.append((i, i+len(motif)-1))
    else:
        search_list = expand_motif(motif)
        for motif_i in search_list:
            indices = all_indices(sequence, motif_i)
            for ind in indices:
                pool.append((ind, ind+len(motif_i)-1))

    merges = []
    pool.sort()
    if len(pool)==0:
        return merges
    if len(pool)==1:
        return [pool[0]]
    item_1 = pool[0]
    del pool[0]
    while len(pool)>0:
        item_new = pool[0]
        del pool[0]
        item_1, merged = overlap(item_1, item_new)
        if not merged:
            merges.append(item_1)
            item_1 = item_new
        if not merged and len(pool)==0:
            merges.append(item_new)
        if merged and len(pool)==0:
            merges.append(item_1)
    return merges

def report_chr(data, chr, fw, strand=""):
    """
    Prints out positions from **data** to file **fw**
    """
    positions = data.get(chr, {}).keys()
    if positions==[]:
        return
    positions.sort()
    start = positions[0]
    stop = positions[0]
    for pos_1, pos_2 in zip(positions, positions[1:]):
        if pos_2==pos_1+1:
            stop = pos_2
        else:
            fw.write("%s\t%s\t%s\t%s%s\n" % (chr, start, stop+1, strand, stop-start+1))
            start = pos_2
            stop = pos_2
    
def cluster_threshold(motif):
    """
    Clusters and thresholds **data** from **motif**. The data is read-in from the .bed file
    """
    rlen = 0
    for (chrom, strand), L in m3_light.db.regions_chrom.items():
        for (region_id, start, stop, region_class) in L:
            rlen += stop-start+1
    
    bed_filename = m3_light.config.filename_bed(motif)
    bg = m3_light.utils.Bedgraph(bed_filename)
    bg.cluster(m3_light.config.cluster_hw) # cluster positions (by strand) with half window cluster_w
    
    hc = {}
    # histogram of values by region
    for (chrom, strand), L in m3_light.db.regions_chrom.items():
        for (region_id, start, stop, region_class) in L:
            for i in range (start, stop):
                v = bg.get_value(chrom, strand, i)
                if v>0:
                    hc[v] = hc.setdefault(v, 0) + 1
    
    f_stat = open(m3_light.config.filename_stats(), "at")
    
    h_choosen = {}
    for pth in m3_light.config.pth:
        h_choosen[pth] = {}
    
    f_stat.write("regions_length=%s\n" % rlen)
    f_stat.write("motif=%s\n" % motif)
    
    for key in hc.keys():
        greater_equal = sum([hc[x] for x in hc.keys() if x>=key])
        greater_equal_p = float(greater_equal)/rlen*100
        f_stat.write("|h>=%s|=%s (%.3f %%)\n" % (key, greater_equal, greater_equal_p))
    f_stat.write("\n")

    for pth in m3_light.config.pth:
        distances = []
        for key in hc.keys():
            greater_equal = sum([hc[x] for x in hc.keys() if x>=key])
            greater_equal_p = float(greater_equal)/rlen*100
            distances.append((abs(pth-greater_equal_p), key))
        distances.sort()
        h_choosen[pth] = max(m3_light.config.h_min, distances[0][1])

    for pth in m3_light.config.pth:
        data_extended_plus_filtered = {}
        data_extended_minus_filtered = {}
        for (chrom, strand), L in m3_light.db.regions_chrom.items():
            for (region_id, start, stop, region_class) in L:
                for i in xrange(start, stop):
                    v = bg.get_value(chrom, strand, i)
                    if v>=h_choosen[pth]:
                        if strand=="+":
                            data_extended_plus_filtered.setdefault(chrom, {}).setdefault(i, 1)
                        else:
                            data_extended_minus_filtered.setdefault(chrom, {}).setdefault(i, 1)

        h_str = ["p%s_h%s" % (key, val) for key, val in h_choosen.items()]
        h_str = "_".join(h_str)
    
        fw = open(m3_light.config.filename_pth(pth, "%s" % motif), "wt")
        # write choosen h_min to .bed file
        # fw.write("#%s\n" % h_str)
        chrs = set(data_extended_plus_filtered.keys()).union(data_extended_minus_filtered.keys())
        chrs = list(chrs)
        chrs.sort()
        for chr in chrs:
            report_chr(data_extended_plus_filtered, chr, fw, strand="+")
            report_chr(data_extended_minus_filtered, chr, fw, strand="-")
        fw.close()
        
def make_bed(motif):
    """
    Creates .bed file for motif
    """
    f = open(m3_light.config.filename_bed(motif), "wt")
    chrs = m3_light.db.regions_chrom.keys()
    chrs.sort()
    for (chrom, strand) in chrs:
        results = []
        for (region_id, start, stop, region_class) in m3_light.db.regions_chrom[(chrom, strand)]:
            start_ext = start-m3_light.config.cluster_hw*2
            stop_ext = stop+m3_light.config.cluster_hw*2
            gs = m3_light.genomes.get_chromosome_subseq(m3_light.config.genome, chrom, start_ext, stop_ext)
            if gs==None:
                continue
            if strand=="+":
                motifs = m3_light.search.find_motif(gs, motif)
            if strand=="-":
                motifs = m3_light.search.find_motif(gs, m3_light.genomes.reverse_complement(motif))
            motifs.sort()
            for (mstart, mstop) in motifs:
                results.append((start_ext+mstart, start_ext+mstop))
        # merge results (in case search regions overlapped)
        for start, stop in m3_light.utils.merge_intervals(results):
            f.write("%s\t%s\t%s\t%s1\n" % (chrom, start, stop+1, strand))
    f.close()        

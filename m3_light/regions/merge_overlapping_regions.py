"""
Merges overlapping regions on same strand and chromosome.

Example input:

    chr1	1000	2000	+
    chr1	1500	2500	+
    chr1	5000	6000	+
    
Example output:

    chr1	1000	2500	+
    chr1	5000	6000	+

Output is sorted.

"""

import sys
import os
from operator import itemgetter, attrgetter

input_file = sys.argv[1]
output_file = sys.argv[2]

def sort_file(file_in):
    f = open(file_in, "rt")
    header = f.readline()
    r = f.readline()
    recs = []
    while r:
        r = r.rstrip("\r").rstrip("\n").split("\t")
        r[3] = int(r[3])
        recs.append(r)
        r = f.readline()
    f.close()
    recs = sorted(recs, key=itemgetter(1,3))
    f_temp = open(file_in+"_temp", "wt")
    f_temp.write(header)
    for rec in recs:
        f_temp.write("\t".join(str(el) for el in rec)+"\n")
    f_temp.close()
    os.remove(file_in)
    os.rename(file_in+"_temp", file_in)

def overlap_yesno(item1, item2, type="all"):
    start1 = int(item1[0])
    start2 = int(item2[0])
    stop1 = int(item1[1])
    stop2 = int(item2[1])
    if stop2>=start1 and stop2<=stop1:
        return True
    if start2>=start1 and start2<=stop1:
        return True
    return False
    
def join(item1, item2):
    start1 = int(item1[0])
    start2 = int(item2[0])
    stop1 = int(item1[1])
    stop2 = int(item2[1])
    start_new = min(start1, start2, stop1, stop2)
    stop_new = max(start1, start2, stop1, stop2)
    return (start_new, stop_new)
    
def read_regions(filename):
    region_id = 0
    regions = {}
    f = open(filename, "rt")
    r = f.readline()
    while r:
        r = r.rstrip("\r").rstrip("\n").split("\t")
        chr = r[0]
        list_start = int(r[1])
        list_stop = int(r[2])
        r_start = min(list_start, list_stop)
        r_stop = max(list_start, list_stop)
        strand = r[3]
        regions[region_id] = (chr, strand, r_start, r_stop)
        region_id += 1
        r = f.readline()
    f.close()
    return regions

regions = read_regions(input_file)
region_pool = {}
for rid, (chr, strand, list_start, list_stop) in regions.iteritems():
    reglist = region_pool.get(chr+"_"+strand, [])
    reglist.append((list_start, list_stop))
    region_pool[chr+"_"+strand] = reglist

non_overlapping = {}
for chr_strand, reglist in region_pool.iteritems():
    reglist.sort()
    regions = []
    while len(reglist)>0:
        item1 = reglist[0]
        del reglist[0]
        overlaps = 1
        while overlaps==1:
            overlaps = -1
            for c, item2 in enumerate(reglist):
                if overlap_yesno(item1, item2):
                    overlaps_index = c
                    overlaps = 1
                    break
            if overlaps==1:
                item1 = join(item1, reglist[overlaps_index])
                del reglist[overlaps_index]
        regions = non_overlapping.get(chr_strand, [])
        regions.append(item1)
        non_overlapping[chr_strand] = regions

# write output
f = open(output_file, "wt")
f.write("\t".join(["id", "chrom", "strand", "start", "stop", "class"]) + "\n")
id = 1
for chr_strand, reglist in non_overlapping.iteritems():
    reglist.sort()
    chr = chr_strand.split("_")[0]
    strand = chr_strand.split("_")[1]
    for (start, stop) in reglist:
        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (id, chr, strand, start, stop, "case"))
        id += 1
f.close()

sort_file(sys.argv[2])
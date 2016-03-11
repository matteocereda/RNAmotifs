"""
...exon|a-----b|exon|c-----d|exon...

a;b;c;d

(b-c) = IN position
(a-d) = SKIP position
"""
    
import os
import sys

f = open(sys.argv[1], "rt")
fout = open(sys.argv[2], "wt")
r = f.readline()
r = f.readline()
allowed_pos = {}
while r:
    r = r.rstrip("\r").rstrip("\n").split(";")
    chr = r[2]
    if chr.find("random")!=-1:
        r = f.readline()
        continue
    strand = r[3]
    rskip_from = int(r[4])
    rin_from = int(r[5])
    rin_to = int(r[6])
    rskip_to = int(r[7])
    pos_from = rskip_from-100
    pos_to = rskip_from+500
    fout.write("%s\t%s\t%s\t%s\n" % (chr, pos_from, pos_to, strand))
    pos_from = rskip_to-500
    pos_to = rskip_to+100
    fout.write("%s\t%s\t%s\t%s\n" % (chr, pos_from, pos_to, strand))
    pos_from = rin_from-500
    pos_to = rin_from+100
    fout.write("%s\t%s\t%s\t%s\n" % (chr, pos_from, pos_to, strand))
    pos_from = rin_to-100
    pos_to = rin_to+500
    fout.write("%s\t%s\t%s\t%s\n" % (chr, pos_from, pos_to, strand))
    r = f.readline()
fout.close()

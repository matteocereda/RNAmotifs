"""
Handles genomes information
"""

import m3_light
import os

code = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"]}
revCode = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
revCodeRYSW = {'R' : 'Y', 'Y' : 'R', 'S' : 'S', 'W' : 'W'}
motifs = []
motifs_redundant = []

def get_chromosome_subseq(genome, chromosome, start_pos, end_pos):
    """
    :param genome: mm9, hg19, ...
    :param chromosome: chr1, chrX, ...
    :param start_pos: 0-based
    :param end_pos: 0-based, right closed

    Returns chromosome sequence from start_pos to end_pos. If strand is negative, return reverse complement
    """
    fname = "%s/genomes/%s/%s.string" % (m3_light.config.folder_root, m3_light.config.genome, chromosome)
    if not os.path.exists(fname):
        return None
    f = open(fname, "rt")
    f.seek(start_pos, 0)
    seq = f.read( end_pos - start_pos + 1)
    return seq
 
def reverse_complement(s):
    """
    :param s: sequence to reverse complement
    
    Reverse complement sequence
    """
    reverse_str = lambda s: ''.join([s[i] for i in xrange(len(s)-1, -1, -1)])
    rs = reverse_str(s)
    new_sequence = ""
    for c in rs:
        rc = revCode.get(c.upper(), c.upper())
        if c.islower():
            rc = rc.lower()
        new_sequence += rc
    return new_sequence
    
for x1 in ["A", "C", "T", "G"]:
    for x2 in ["A", "C", "T", "G"]:
        for x3 in ["A", "C", "T", "G"]:
            for x4 in ["A", "C", "T", "G"]:
                motifs.append(x1+x2+x3+x4)

for x1 in ["A", "C", "T", "G"]:
    for x2 in ["A", "C", "T", "G"]:
        for B in ["R", "Y", "S", "W"]:
            for E in ["R", "Y", "S", "W"]:
                motif = B+x1+x2+E
                motifs_redundant.append(motif)

motifs = motifs + motifs_redundant

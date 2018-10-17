"""
Main m3 interface
"""

import config
import regions
import search
import genomes
import utils
import db
import m3_light
import os
import data

def find_tetramers(regions_name, genome, cluster_hw, hmin, perc_th, m=None):
    """
    :param regions_name: name of the regions to process
    :param genome: genome of given regions (e.g. mm9)
    :param m: core to use

    Reads regions list from .tab file in **regions/<regions_name>/<regions_name>.tab** and generate .bed and other files
    """
    m3_light.config.regions_name = regions_name
    m3_light.config.genome = genome
    m3_light.config.cluster_hw = int(cluster_hw)
    m3_light.config.h_min = int(hmin)
    # m3_light.config.pth[0] = float(perc_th)
    # print(m3_light.config.pth)

    m3_light.regions.read()

    # remove stats file
    if os.path.exists(m3_light.config.filename_stats()):
        os.remove(m3_light.config.filename_stats())

    for index, motif in enumerate(genomes.motifs):
        print index+1, motif
        print "\tmotifs.raw"
        m3_light.search.make_bed(motif)
        print "\tmotifs.final"
        search.cluster_threshold(motif)

def fastasplit(filename):
    f = m3_light.data.Fasta(filename)
    while f.read():
        print f.id
        fout = open("%s.string" % f.id, "wt")
        fout.write(f.sequence)
        fout.close()

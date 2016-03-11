"""
Configuration options

================ ============== ======
variable         type           value
================ ============== ======
regions_name     str            name of the regions folder and file (regions/<regions_name>/<regions_name>.tab) from which regions will be read
genome           str            which genome to use (e.g. mm9, hg19,...)
cluster_hw       int            half window size used for clustering, default 15
folder_root      str            root folder of m3, e.g. /home/user/m3
pth              []             list of pth values to compute the thresholding for, default [0.1, 0.25, 0.5, 0.75, 1]
bootstrap        int            number of bootstrapping steps, default 1000
h_min            int            minimal h at the thresholding step, default 4
================ ============== ======
"""

import m3_light
import os
from os.path import join as pjoin

regions_name = ""
genome = ""
cluster_hw = 15
h_min = 4
folder_root = os.getcwd()+"/m3_light"
pth = [0.5]
bootstrap = 1000

def filename_regions():
    return "%s/regions/%s/%s.tab" % (m3_light.config.folder_root, m3_light.config.regions_name, m3_light.config.regions_name)

def folder_bed():
    """
    Returns folder in which to create .bed files
    """
    folder = pjoin(folder_root, "results", regions_name, "motifs.raw")
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder
    
def filename_bed(base):
    """
    Returns .bed filename
    """
    return pjoin(folder_bed(), "%s.bed" % base)
    
def filename_stats():
    folder = pjoin(folder_root, "results", regions_name, "motifs.final")
    if not os.path.exists(folder):
        os.makedirs(folder)
    return pjoin(folder, "stats.txt")
    
def folder_pth(pth):
    folder = pjoin(folder_root, "results", regions_name, "motifs.final", "pth_%s" % pth)
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder
    
def filename_pth(pth, tetra):
    return pjoin(folder_pth(pth), "%s.bed" % tetra)
    
def folder_hits(pth):
    folder = pjoin(folder_root, "results", regions_name, "motifs.final", "pth_%s" % pth)
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def filename_hits(pth, tetra):
    return pjoin(folder_hits(pth), "%s.bed" % tetra)

def filename_fisher(pth):
    return pjoin(folder_root, "results", regions_name, "fisher_pth_%s.tab" % pth)

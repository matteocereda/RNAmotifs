"""
Stores regions .tab files and reads the regions into .db module.
"""

import m3_light

def read():
    """
    Read regions from .tab file, stored in **regions/<m3.config.regions_name>/<m3.config.regions_name>.tab**
    """
    m3_light.db.reset()
    filename = m3_light.config.filename_regions()
    print "reading regions from:", filename
    f = open(filename, "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        region_id = int(r[header.index("id")])
        chrom = r[header.index("chrom")]
        strand = r[header.index("strand")]
        start = int(r[header.index("start")])
        stop = int(r[header.index("stop")])
        region_class = r[header.index("class")]
        m3_light.db.regions.append((region_id, region_class))
        m3_light.db.regions_chrom.setdefault((chrom, strand), []).append((region_id, start, stop, region_class))
        r = f.readline()
    f.close()
    m3_light.db.regions.sort()


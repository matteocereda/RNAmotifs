"""
Database for regions that are read and parsed from a tab delimited input file:

================ ============== ======
variable         type           value
================ ============== ======
regions          []             [(region_id, region_class), ...]
regions_chrom    {}             {(chrom, strand) : [()]}
================ ============== ======

"""

import m3_light

def reset():
    """
    Clears all data
    """
    
    m3_light.db.regions = []
    m3_light.db.regions_chrom = {}
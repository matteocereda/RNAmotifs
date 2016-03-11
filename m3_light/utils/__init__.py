def merge_intervals(data):
    """
    data = [(10,20), (15,30), (100, 200)]
    out = [(10,30), (100,200)]
    """
    if len(data)==0:
        return data
    result = []
    saved = list(data[0])
    for st, en in sorted([sorted(t) for t in data]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            result.append((saved[0], saved[1]))
            saved[0] = st
            saved[1] = en
    result.append(saved)
    return result

class Bedgraph():
    def __init__(self, filename):
        self.data = {}
        self.sites = 0
        self.strand_file = {"+":"", "-": "-"}
        self.load(filename)

    def load(self, filename):
        if filename.endswith(".gz"):
            f = gzip.open(filename, "rb")
        else:
            f = open(filename, "rb")
        r = f.readline()
        while r:
            if r.startswith("#"):
                r = f.readline()
                continue
            r = r.rstrip("\r").rstrip("\n").split("\t")
            chr = r[0]
            start = int(r[1])
            stop = int(r[2])
            cDNA = float(r[3])
            strand = "+" if cDNA>=0 else "-"
            self.data.setdefault(chr, {}).setdefault(strand, {})
            for p in xrange(start, stop):
                self.data[chr][strand][p] = abs(cDNA)
            r = f.readline()
        f.close()

    def chromosomes(self):
        return self.data.keys()

    def get_value(self, chr, strand, pos):
        return self.data.get(chr, {}).get(strand, {}).get(pos, 0)

    def region(self, chr, strand, pos_from, pos_to):
        return sum([self.get_value(chr, strand, i) for i in xrange(pos_from, pos_to+1)])

    def cluster(self, hws):
        data_window = {}
        self.data_sum = 0
        for chr, chr_data in self.data.iteritems():
            data_window[chr] = {}
            for strand, strand_data in chr_data.iteritems():
                data_window[chr][strand] = {}
                for pos, value in strand_data.iteritems():
                    pos_from = max(0, pos-hws)
                    pos_to = max(0, pos+hws)
                    data_window[chr][strand][pos] = self.region(chr, strand, pos_from, pos_to)
        self.data = data_window
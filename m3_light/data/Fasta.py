import m3_light
import gzip
import os

class Fasta:

    """
    This class handles FASTA files.
    """

    def __init__(self, filename):
        self.filename = filename
        self.open_file()
                
    def open_file(self):
        self.id = None
        self.next_id = None
        if self.filename.endswith(".gz"):
            self.f = gzip.open(self.filename, "rt")
        else:
            self.f = open(self.filename, "rt")
        
    def read(self):
        self.sequence = []
        if self.next_id!=None:
            self.id = self.next_id
            self.next_id = None
        r = self.f.readline()
        if not r:
            return False
        while r:
            r = r.replace("\r", "").replace("\n", "")
            if r=="":
                r = self.f.readline()
                continue
            if r[0]==">" and self.id==None:
                self.id = r[1:]
                r = self.f.readline()
                continue
            elif r[0]==">":
                self.next_id = r[1:]
                self.sequence = "".join(self.sequence)
                return True
            self.sequence.append(r)
            r = self.f.readline()
        self.sequence = "".join(self.sequence)
        return True

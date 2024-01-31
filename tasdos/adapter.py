import os
import sys
from tasdos.utils import time_stamp

class Adapter(object):

    def __init__(self, outdir, adapter, adapterfile):
        self.dir = outdir
        self.ad = adapter
        self.adfile = adapterfile
        
    def nextera(self):
        a = ('>adapter/1\n'
             'CTGTCTCTTATACACATCT\n'
             '>adapter/2\n'
             'CTGTCTCTTATACACATCT\n')
        return a
    
    def truseq(self):
        a = ('>adapter/1\n'
             'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n'
             '>adapter/2\n'
             'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n')
        return a
    
    def custom(self):
        if not os.path.isfile(self.adfile):
            print(time_stamp(), 
                  '!!ERROR!! Custom adapter file does not exist.\n', flush=True)
            sys.exit(1)
        with open(self.adfile, 'r') as f:
            a = f.read()
        return a
    
    def run(self):
        filename = '{}/adapter/adapter.fa'.format(self.dir)
        with open(filename, 'w') as f:
            if self.ad == 'NEXTERA':
                f.write(self.nextera())
            elif self.ad == 'TRUSEQ':
                f.write(self.truseq())
            elif self.ad == 'CUSTOM':
                f.write(self.custom())


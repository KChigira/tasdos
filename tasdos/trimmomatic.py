import glob
import sys
import subprocess as sbp
from tasdos.utils import prepare_cmd, call_log, time_stamp


class Trimmomatic(object):

    def __init__(self, fqdir, samplename, outdir, crop, minlen, adapter):
        self.fqdir = fqdir
        self.sn = samplename
        self.dir = outdir
        self.crop = crop
        self.minlen = minlen
        self.ad = adapter

    def run(self):
        #Check the input files are pair end
        tmp = '{}/{}_*.fastq.gz'.format(self.fqdir, self.sn)
        filelist = glob.glob(tmp)
        if len(filelist) != 2:
            print(time_stamp(), 
                  '!!ERROR!! Input fastq are not pair end or more than 2 files exist for one sample name.\n', flush=True)
            sys.exit(1)
            

        if self.ad == 'NONE':
            cmd01 = 'cp {0}/{1}*_R1_*.fastq.gz \
                     {2}/fastq/{1}_R1_trimmed.fastq.gz \
                     2>&1'.format(self.fqdir, self.sn, self.dir)
            cmd02 = 'cp {0}/{1}*_R2_*.fastq.gz \
                     {2}/fastq/{1}_R2_trimmed.fastq.gz \
                     2>&1'.format(self.fqdir, self.sn, self.dir)
            try:
                sbp.run(cmd01,
                        stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                        shell=True, check=True)
                sbp.run(cmd02,
                        stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                        shell=True, check=True)
            except sbp.CalledProcessError:
                print(time_stamp(), 
                      '!!ERROR!! Moving fastq was failed.\n', flush=True)
                sys.exit(1)
        else:
            cmd1 = 'trimmomatic PE -threads 1 -phred33 \
                    {0}/{1}*_R1_*.fastq.gz \
                    {0}/{1}*_R2_*.fastq.gz \
                    {2}/fastq/{1}_R1_trimmed.fastq.gz \
                    {2}/fastq/{1}_R1_unpaired.fastq.gz \
                    {2}/fastq/{1}_R2_trimmed.fastq.gz \
                    {2}/fastq/{1}_R2_unpaired.fastq.gz \
                    ILLUMINACLIP:{2}/adapter/adapter.fa:2:30:10 \
                    LEADING:20 \
                    TRAILING:20 \
                    SLIDINGWINDOW:30:20 \
                    CROP:{3} \
                    MINLEN:{4} \
                    >> {2}/log/trimmmomatic.log \
                    2>&1'.format(self.fqdir, self.sn, self.dir, 
                                self.crop, self.minlen)
            cmd1 = prepare_cmd(cmd1)

            try:
                sbp.run(cmd1,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.dir, 'trimmmomatic', cmd1)
                sys.exit(1)


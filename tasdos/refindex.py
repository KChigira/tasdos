import os
import sys
import subprocess as sbp
from tasdos.utils import time_stamp, prepare_cmd, call_log


class RefIndex(object):

    def __init__(self, ref, dir):
        self.ref = ref
        self.dir = dir

    def run(self):
        print(time_stamp(),
              'Indexing reference fasta.',
              flush=True)

        cmd1 = 'samtools faidx {} \
                >> {}/log/samtools.log \
                2>&1'.format(self.ref, self.dir)
        cmd1 = prepare_cmd(cmd1)
        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.dir, 'samtools', cmd1)
            sys.exit(1)


        if os.path.isfile('{}.dict'.format(os.path.splitext(self.ref)[0])):
            pass
        else:
            cmd2 = 'picard CreateSequenceDictionary R={} \
                    >> {}/log/picard.log \
                    2>&1'.format(self.ref, self.dir)
            cmd2 = prepare_cmd(cmd2)
            try:
                sbp.run(cmd2,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.dir, 'picard', cmd2)
                sys.exit(1)

        if os.path.isfile('{}.bwt'.format(self.ref)):
            pass
        else:
            cmd3 = 'bwa index {} \
                    >> {}/log/bwa.log \
                    2>&1'.format(self.ref, self.dir)
            cmd3 = prepare_cmd(cmd3)
            try:
                sbp.run(cmd3,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.dir, 'bwa', cmd3)
                sys.exit(1)

        print(time_stamp(),
              'Done.',
              flush=True)
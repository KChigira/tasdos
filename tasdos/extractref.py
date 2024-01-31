import sys
import subprocess as sbp
import pandas as pd
from tasdos.refindex import RefIndex
from tasdos.utils import time_stamp, prepare_cmd, call_log


class ExtractRef(object):

    def __init__(self, ref, vcf, dir):
        self.ref = ref
        self.vcf = pd.DataFrame(vcf)
        self.dir = dir

    def run(self):
        print(time_stamp(),
              'Extracting reference sequence data.',
              flush=True)

        #VCF for only target variant.
        #Modify column names of samples
        new_vcf = self.vcf
        for i in range(9, len(new_vcf.columns)):
            cn = new_vcf.columns[i]
            newcn = 'template_{}'.format(cn)
            new_vcf = new_vcf.rename(columns={cn : newcn})

        #Making target refarence
        fn_newref = '{}/ref/extracted_ref.fasta'.format(self.dir)
        for i in range(len(self.vcf)):
            target_chr = self.vcf.at[self.vcf.index[i], '#CHROM']
            target_pos = int(self.vcf.at[self.vcf.index[i], 'POS'])
            #Make strings such as "chr01:100000-100400"
            query = '{}:{}-{}'.format(target_chr,
                                      target_pos - 200,
                                      target_pos + 200)
            cmd1 = 'samtools faidx {} {} \
                    >> {}'.format(self.ref, query, fn_newref)
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

            #Modify the #CHROM and POS of target VCF
            new_vcf.at[new_vcf.index[i], '#CHROM'] = query
            new_vcf.at[new_vcf.index[i], 'POS'] = 201

        #Add index to new reference
        ri = RefIndex(fn_newref, self.dir)
        ri.run()

        #Output new VCF
        fn = '{}/ref/target.vcf'.format(self.dir)
        with open(fn, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
        new_vcf.to_csv(fn, sep='\t', mode='a', index=False)
        #add Index to vcf
        cmd2 = 'gatk IndexFeatureFile \
                -F {} \
                >> {}/log/gatk.log 2>&1'.format(fn, self.dir)
        cmd2 = prepare_cmd(cmd2)
        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.dir, 'gatk', cmd2)
            sys.exit(1)

        print(time_stamp(),
              'Done.',
              flush=True)
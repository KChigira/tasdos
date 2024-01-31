import sys
import pandas as pd
from tasdos.utils import time_stamp, read_vcf


class MergeVcf(object):

    def __init__(self, dir, samplelist):
        self.dir = dir
        self.sl = samplelist

    def run(self):
        print(time_stamp(), 'Merging vcf files.', flush=True)

        #Prepare frame of VCF from template
        tmp = read_vcf('{}/ref/target.vcf'.format(self.dir))
        mergedvcf = pd.DataFrame(tmp[2], columns=tmp[1])

        #Add genotypes for each sample
        for i in range(len(self.sl)):
            fn = '{}/vcf/{}_selected_variants.vcf'.format(self.dir, self.sl[i])
            tmp = read_vcf(fn)
            svcf = pd.DataFrame(tmp[2], columns=tmp[1])
            mergedvcf = pd.concat([mergedvcf, svcf[svcf.columns[9]]], axis=1)
            #row 9 is genotypes

        #Add command info
        command = ' '.join(sys.argv)

        #Output new VCF
        fn = '{}/Result_TAS_A.vcf'.format(self.dir)
        with open(fn, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##command={}\n'.format(command))
        mergedvcf.to_csv(fn, sep='\t', mode='a', index=False)

        print(time_stamp(),
              'Done.',
              flush=True)


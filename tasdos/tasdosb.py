#!/usr/bin/env python3
import os
import sys
import pandas as pd
from tasdos.params import Params
from tasdos.utils import read_vcf, time_for_filename, time_stamp

pm = Params('tasdosb')
args = pm.set_options()

class TasdosB(object):
    def __init__(self, args):
        self.args = args
        self.p1 = args.parent1
        self.p2 = args.parent2
        self.mindep = args.mindep
        self.hchi = args.hetero_chi
        self.noise = args.noise_level
        #This is the name of output directory
        self.dir = 'tasdosB_{}'.format(time_for_filename())
        os.mkdir(self.dir)

        #make input VCF to dataframe
        tmp = read_vcf(args.input)
        self.header = tmp[0] #List
        self.col = tmp[1] #List
        self.vcf = pd.DataFrame(tmp[2], columns=self.col) #dataframe

        self.s_cols = [] #Columun number of samples
        self.p1_col = -1 #Columun number of templates
        self.p2_col = -1

        self.out = [] #output frame
        self.out_dep = [] #output frame

    def command(self):
        #Output command info
        command = ' '.join(sys.argv)
        fn = '{}/command.txt'.format(self.dir)
        with open(fn, 'w') as f:
            f.write('{}\n'.format(command))

    def checkargs(self):
        flag1 = False #Check parents name
        flag2 = False #Check parents name
        for i in range(9, len(self.col)):
            n = self.col[i]
            if n[0:9] == 'template_':
                if n[9:] == self.p1:
                    flag1 = True
                    self.p1_col = i
                if n[9:] == self.p2:
                    flag2 = True
                    self.p2_col = i
            else:
                self.s_cols.append(i)
        if flag1 and flag2:
            pass
        else:
            print(time_stamp(), 
                  '!!ERROR!! Names of parents are wrong.\n', flush=True)
            sys.exit(1)

    def makeframe(self):
        all_list = list(self.vcf['#CHROM'])
        chr_list = []
        pos_list = []
        for l in all_list:
            tmp1 = l.split(':')
            chr = tmp1[0]
            tmp2 = tmp1[1].split('-')
            pos = int((int(tmp2[0]) + int(tmp2[1])) / 2)
            chr_list.append(chr)
            pos_list.append(pos)
        chr_list = pd.Series(data = chr_list, name = 'CHR')
        pos_list = pd.Series(data = pos_list, name = 'POS')
        self.out = pd.concat([chr_list, pos_list], axis=1)
        self.out_dep = self.out.copy()

        #Make the names of markers
        names = []
        for i in range(len(self.out)):
            chr = self.out.iloc[i, 0]
            pos = int(self.out.iloc[i, 1])
            pos_n = f"{round(pos / 10000):04}" #Fill by 0 to be 4 digits
            name = '{}_{}'.format(chr, pos_n)
            names.append(name)
        names = pd.Series(data = names, name = 'NAME')
        self.out = pd.concat([self.out, names], axis=1)

    def convert(self):
        #get Genotypes for each sample
        for i in self.s_cols:
            sn = self.vcf.columns[i]
            res = []
            dep = []
            for j in range(len(self.vcf)):
                #Check haplotype calling format.
                format = self.vcf.iloc[j, 8].split(':')
                gt_index = format.index('GT')
                ad_index = format.index('AD')

                #which genotype A or B
                gt_p1 = self.vcf.iloc[j, self.p1_col].split(':')[gt_index]
                gt_p2 = self.vcf.iloc[j, self.p2_col].split(':')[gt_index]
                if gt_p1 == '0/0' and gt_p2 == '1/1':
                    key = 1 #'0/0' = 'A', '1/1' = 'B'
                elif gt_p1 == '1/1' and gt_p2 == '0/0':
                    key = 2 #'0/0' = 'B', '1/1' = 'A'
                else:
                    res.append('-')
                    dep.append(0)
                    continue

                #analyze for sample genotyping
                data_s = self.vcf.iloc[j, i].split(':') #ex. ['0/0', '12,0', ....]
                if len(data_s) != len(format):
                    res.append('-')
                    dep.append(0)
                    continue
                ad_s = data_s[ad_index].split(',') #ex. [12, 0]
                ref = int(ad_s[0])
                alt = int(ad_s[1])
                dep.append(ref + alt)

                if ref + alt < self.mindep:
                    res.append('-') #Thin depth marker is missing
                    continue

                #Is this marker HOMO?
                minor_freq = min([ref, alt]) / (ref + alt) #ex. [2, 18] --> 0.1
                if minor_freq <= self.noise: #default: noise=0.1
                    #HOMO
                    if key == 1:#'0/0' = 'A', '1/1' = 'B'
                        if ref > alt: #'0/0' = 'A'
                            res.append('A')
                        else: #'1/1' = 'B'
                            res.append('B')
                    if key == 2:#'0/0' = 'B', '1/1' = 'A'
                        if ref > alt: #'0/0' = 'B'
                            res.append('B')
                        else: #'1/1' = 'A'
                            res.append('A')
                else:
                    #Hetero or Missing
                    #chi-squared test
                    tv = (ref + alt) / 2; #theorical varue when segregeted 1:1
                    chi2 = ((ref - tv)**2 + (alt - tv)**2) / tv #chi-squared

                    if ref < 5 or alt < 5:
                        #chi2 test is not suitable for either value is under 5
                        res.append('-')
                    elif chi2 < self.hchi:
                        #Hetero
                        res.append('H')
                    else:
                        res.append('-')
            
            res = pd.Series(data = res, name = sn)
            dep = pd.Series(data = dep, name = sn)
            self.out = pd.concat([self.out, res], axis=1)
            self.out_dep = pd.concat([self.out_dep, dep], axis=1)
        
        #write tsv
        fn = '{}/Result_TAS_B.tsv'.format(self.dir)
        self.out.to_csv(fn, sep='\t', mode='w', index=False)
        fn = '{}/Result_TAS_B_depth.tsv'.format(self.dir)
        self.out_dep.to_csv(fn, sep='\t', mode='w', index=False)

def main():
    print(time_stamp(), 'tasdos B started.', flush=True)

    prog = TasdosB(args)
    prog.command()
    prog.checkargs()
    prog.makeframe()
    prog.convert()

    print(time_stamp(), 'tasdos B successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()

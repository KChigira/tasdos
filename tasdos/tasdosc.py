#!/usr/bin/env python3
import os
import sys
import pandas as pd
from tasdos.params import Params
from tasdos.utils import time_for_filename, time_stamp

pm = Params('tasdosc')
args = pm.set_options()

class TasdosC(object):
    def __init__(self, args):
        self.args = args
        self.p1 = args.parent_sample1
        self.p2 = args.parent_sample2
        self.missrate = args.missing_rate
        self.mfreq = args.minor_freq
       
        #This is the name of output directory
        self.dir = 'tasdosC_{}'.format(time_for_filename())
        os.mkdir(self.dir)

        #make input TSV to dataframe
        self.data = pd.read_csv(args.input, sep='\t') #dataframe

        self.s_cols = [] #Columun number of samples
        self.p1_col = -1 #Columun number of parent sample1
        self.p2_col = -1 #Columun number of parent sample2

        self.out = [] #output frame

        #Below are all for judge remove or remain
        self.A_list = []
        self.B_list = []
        self.H_list = []
        self.M_list = []
        self.A_ratio = []
        self.B_ratio = []
        self.H_ratio = []
        self.M_ratio = []
        self.judge_missrate = [] #boolean
        self.judge_mfreq = [] #boolean
        self.judge_parent = []
        self.judge_all = []

    def command(self):
        #Output command info
        command = ' '.join(sys.argv)
        fn = '{}/command.txt'.format(self.dir)
        with open(fn, 'w') as f:
            f.write('{}\n'.format(command))

    def checkargs(self):
        if self.p1 is None or self.p2 is None:
            self.s_cols = range(3, len(self.data.columns))
            return
        
        flag1 = False #Check parents name
        flag2 = False #Check parents name
        for i in range(3, len(self.data.columns)):
            n = self.data.columns[i]
            if n == self.p1:
                flag1 = True
                self.p1_col = i
            elif n == self.p2:
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

    def countalleles(self):
        for i in range(len(self.data)):
            A_cnt = 0
            B_cnt = 0
            H_cnt = 0
            M_cnt = 0
            for j in self.s_cols:
                allele = self.data.iloc[i, j]
                if allele == 'A':
                    A_cnt = A_cnt + 1
                elif allele == 'B':
                    B_cnt = B_cnt + 1
                elif allele == 'H':
                    H_cnt = H_cnt + 1
                else:
                    M_cnt = M_cnt + 1
            self.A_list.append(A_cnt)
            self.B_list.append(B_cnt)
            self.H_list.append(H_cnt)
            self.M_list.append(M_cnt)
            valid = (A_cnt + B_cnt + H_cnt)
            if valid == 0:
                A_freq = 0
                B_freq = 0
                H_freq = 0
            else:
                A_freq = A_cnt / (A_cnt + B_cnt + H_cnt)
                B_freq = B_cnt / (A_cnt + B_cnt + H_cnt)
                H_freq = H_cnt / (A_cnt + B_cnt + H_cnt)
            Missrate = M_cnt / len(self.s_cols)
            self.A_ratio.append(A_freq)
            self.B_ratio.append(B_freq)
            self.H_ratio.append(H_freq)
            self.M_ratio.append(Missrate)
            if A_freq < self.mfreq or B_freq < self.mfreq:
                self.judge_mfreq.append(False)
            else:
                self.judge_mfreq.append(True)
            if Missrate > self.missrate:
                self.judge_missrate.append(False)
            else:
                self.judge_missrate.append(True)

    def checkparents(self):
        for i in range(len(self.data)):
            if self.p1 is None or self.p2 is None or not(self.args.check_parents):
                self.judge_parent.append(True)
                continue
            
            p1_allele = self.data.iloc[i, self.p1_col]
            p2_allele = self.data.iloc[i, self.p2_col]
            if p1_allele == 'A' and p2_allele == 'B':
                self.judge_parent.append(True)
            else:
                self.judge_parent.append(False)

    def judge(self):
        cnt_true = 0
        for i in range(len(self.data)):
            if self.judge_mfreq[i] and self.judge_missrate[i] and self.judge_parent:
                self.judge_all.append(True)
                cnt_true = cnt_true + 1
            else:
                self.judge_all.append(False)
        print(time_stamp(), 
              '{} / {} markers are passed.\n'.format(cnt_true, len(self.data)), flush=True)

    def makeoutput(self):
        add = [self.A_list, self.B_list, self.H_list, self.M_list,
               self.A_ratio, self.B_ratio, self.H_ratio, self.M_ratio,
               self.judge_missrate, self.judge_mfreq,
               self.judge_parent, self.judge_all]
        addcol = ['Count_A', 'Count_B', 'Count_H', 'Count_Missing',
                  'Ratio_A', 'Ratio_B', 'Ratio_H', 'Missing_Rate',
                  'JUDGE_missrate', 'JUDGE_minor_allele',
                  'JUDGE_parents', 'JUDGE']
        addframe = pd.DataFrame(add)
        addframe = addframe.T #inverse
        addframe.columns = addcol
        output1 = pd.concat([self.data, addframe], axis=1)
        #write tsv
        fn = '{}/Result_TAS_C.tsv'.format(self.dir)
        output1.to_csv(fn, sep='\t', mode='w', index=False)

        #format for R/qtl
        #Make dataset of number of chrom.
        chrnum =  []
        chr = 0
        current = ''
        for l in list(output1['CHR']):
            if l != current:
                chr = + chr + 1
                current = l
            chrnum.append(chr)
        chrnum = pd.Series(data = chrnum, name = 'chr')

        #Col[0]= Marker name, Col[1]=number of chrom., Col[2]~=Only sample genotype
        output2 = pd.concat([output1['NAME'], chrnum, output1[output1.columns[self.s_cols]]], axis=1)
        #delete removed markers
        output2 = output2[self.judge_all]

        #inverse
        output2 = output2.T.copy()
        output2 = output2.reset_index() #new column 'index' generated
        output2.iloc[0, 0] = 'id'
        output2.iloc[1, 0] = ''
        fn = '{}/Result_TAS_C_formated_for_Rqtl.csv'.format(self.dir)
        output2.to_csv(fn, sep=',', mode='w', index=False, header=False)

def main():
    print(time_stamp(), 'tasdos C started.', flush=True)

    prog = TasdosC(args)
    prog.command()
    prog.checkargs()
    prog.countalleles()
    prog.checkparents()
    prog.judge()
    prog.makeoutput()

    print(time_stamp(), 'tasdos C successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()

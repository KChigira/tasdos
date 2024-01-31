import sys
import subprocess as sbp

import pandas as pd
from tasdos.utils import read_vcf, time_stamp, prepare_cmd, call_log

class Mapping(object):

    def __init__(self, dir, samplename):
        self.dir = dir
        self.sn = samplename

    def run(self):
        #Mapping reads using bwa mem
        cmd1 = 'bwa mem -t 1 \
                -R "@RG\\tID:{0}\\tLB:{0}\\tPL:ILLUMINA\\tSM:{0}" \
                {1}/ref/extracted_ref.fasta \
                {1}/fastq/{0}_R1_trimmed.fastq.gz \
                {1}/fastq/{0}_R2_trimmed.fastq.gz \
                > {1}/bam/{0}_aligned_reads.sam'.format(self.sn, self.dir)
        cmd1 = prepare_cmd(cmd1)
        try:
            sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 1 was failed.\n', flush=True)
            sys.exit(1)

        #Sort sam and convert to bam
        cmd2 = 'samtools sort -@ 1 -O bam \
                -o {1}/bam/{0}_aligned_reads.bam \
                {1}/bam/{0}_aligned_reads.sam'.format(self.sn, self.dir)
        cmd2 = prepare_cmd(cmd2)
        try:
            sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 2 was failed.\n', flush=True)
            sys.exit(1)

        #Remove SAM
        cmd3 = 'rm {1}/bam/{0}_aligned_reads.sam'.format(self.sn, self.dir)
        cmd3 = prepare_cmd(cmd3)
        try:
            sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 3 was failed.\n', flush=True)
            sys.exit(1)

        #Add Index to BAM
        cmd4 = 'samtools index {1}/bam/{0}_aligned_reads.bam'.format(self.sn, self.dir)
        cmd4 = prepare_cmd(cmd4)
        try:
            sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 4 was failed.\n', flush=True)
            sys.exit(1)

        #Extract only the primary mapped reads
        cmd5 = 'samtools view -@ 1 -bh -F 256 {1}/bam/{0}_aligned_reads.bam \
                > {1}/bam/{0}_aligned_reads_primary.bam'.format(self.sn, self.dir)
        cmd5 = prepare_cmd(cmd5)
        try:
            sbp.run(cmd5, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 5 was failed.\n', flush=True)
            sys.exit(1)

        #Add Index to BAM
        cmd6 = 'samtools index {1}/bam/{0}_aligned_reads_primary.bam'.format(self.sn, self.dir)
        cmd6 = prepare_cmd(cmd6)
        try:
            sbp.run(cmd6, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 6 was failed.\n', flush=True)
            sys.exit(1)

        #Make VCF
        cmd7 = 'gatk HaplotypeCaller \
                -R {1}/ref/extracted_ref.fasta \
                -I {1}/bam/{0}_aligned_reads_primary.bam \
                -O {1}/vcf/{0}_raw_variants.vcf \
                --alleles {1}/ref/target.vcf \
                --genotyping-mode GENOTYPE_GIVEN_ALLELES'.format(self.sn, self.dir)
        cmd7 = prepare_cmd(cmd7)
        try:
            sbp.run(cmd7, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! Mapping command 7 was failed.\n', flush=True)
            sys.exit(1)
        

        #Extract only target variants from VCF
        tmp = read_vcf('{}/ref/target.vcf'.format(self.dir))
        tmpvcf = pd.DataFrame(tmp[2], columns=tmp[1])

        #Raw VCF has many non-target variants
        raw = read_vcf('{1}/vcf/{0}_raw_variants.vcf'.format(self.sn, self.dir))
        rawvcf = pd.DataFrame(raw[2], columns=raw[1])

        #Make Extracted vcf from template VCF
        newvcf = tmpvcf.iloc[:, 0:10].copy()
        newvcf.iloc[:, 5:10] = '.'
        #Modify name of sample column
        newvcf = newvcf.rename(columns={newvcf.columns[9] : self.sn})

        for i in range(len(newvcf)):
            for j in range(len(rawvcf)):
                #Compare values of #CHROM,POS,REF,ALT
                row_new = ','.join(list(newvcf.iloc[i, 0:4]))
                row_raw = ','.join(list(rawvcf.iloc[j, 0:4]))
                if row_new == row_raw:
                    newvcf.iloc[i, :] = rawvcf.iloc[j, :]
                    break

        #Output new VCF
        fn = '{1}/vcf/{0}_selected_variants.vcf'.format(self.sn, self.dir)
        with open(fn, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
        newvcf.to_csv(fn, sep='\t', mode='a', index=False)

        #add Index to vcf
        cmd8 = 'gatk IndexFeatureFile \
                -F {} \
                >> {}/log/gatk.log 2>&1'.format(fn, self.dir)
        cmd8 = prepare_cmd(cmd8)
        try:
            sbp.run(cmd8,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.dir, 'gatk', cmd8)
            sys.exit(1)

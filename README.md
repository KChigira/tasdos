# tasdos
Down stream analysis of genotyping method using Targeted Amplicon Sequencing
#### version 0.0.1 (2023.1.31)

## Outline

'tasdos' is a tool for down stream analysis of genotyping method using Targeted Amplicon Sequencing.
Input files are pair end fastq and VCF containing target SNP which was made by 'mkdesigner'.
We can get genotype file which is suitable to QTL analysis using R/qtl software just in 3 steps.

#### Citation
- Comming soon.

## Install

### Install via PyPI
```
pip install tasdos
```

#### Dependencies
 - python >=3.8,<4.0
 - pandas >=2.0.2,<3.0.0
 - samtools >=1.6,<2.0
 - gatk4 >=4.4.0.0,<5.0.0.0
 - picard >=2.18.29,<3.0.0
 - bwa >=0.7.17

These must be installed manually.

#### Check about packages in dependency
Dependent packages often get errors about shared libraries.    
Please refer to the README of 'mkdesigner' if you got error when you typed like below.
```
samtools --version
```

## Usage
#### Tutorial
(1) Haplotype calling   
```
tasdosA -I directory_of_input_fastq \
        -R reference_genome.fasta \
        -V target_VCF_made_by_mkdesigner.vcf \
        --cpu 6 --adapter NEXTERA \
        --seqlen 150 --minlen 60
```
Result file will be located as   
'./tasdosA_00000000000000/Result_TAS_A.vcf'    

(2) Analize genotype of each samples    
```
tasdosB -I tasdosA_00000000000000/Result_TAS_A.vcf \
        -p1 (Name of parent A in the input VCF of 'tasdosA') \
        -p2 (Name of parent B in the input VCF of 'tasdosA') \
        --mindep 6
```
Result file will be located as    
'./tasdosB_00000000000000/Result_TAS_B.tsv'    

(3) Filter markers by designaetd thresholds    
```
tasdosC -I tasdosB_00000000000000/Result_TAS_B.tsv \
        --parent_sample1 (Name of a sample of parent A, if incruded) \
        --parent_sample2 (Name of a sample of parent B, if incruded) \
        --missing_rate 0.2 \
        --minor_freq 0.1
```
Result file will be located as     
'./tasdosC_00000000000000/Result_TAS_C_formated_for_Rqtl.csv'    
This file can be used as input file of R/qtl.    

#### Commands
tasdosA    
usage: tasdosA -I <Directory containing input FASTQ>   
        -R <File of reference FASTA>    
        -V <File of target VCF>    
        ...    
options:   
  -h, --help      show this help message and exit   
  -I , --input    Directory containing input FASTQ.   
                  This directory must contain only fastq file used in genotyping.   
                  gzip (fastq.gz) also supported.   
  -R , --ref      File of reference FASTA.   
  -V , --vcf      File of target VCF.   
                  (VCF made by mkselect is recommended.)   
  --cpu           Number of CPUs to use.   
  --adapter       Adapter sequences used for trimming fastq.  
                  NONE means the input fastq has already trimmed.   
                  When CUSTOM designated, --adapterfile must be specified.   
  --adapterfile   This is valid when --adapter = CUSTOM.   
  --seqlen        Sequence length of fastq.   
  --minlen        Ignore reads which are shorter than this value after trimming.   
  -v, --version   show program's version number and exit   
  
  
tasdosB   
usage: tasdosB -I <VCF file which is the output of tasdosA>   
        -p1 <Parent name genotyped as A>   
        -p2 <Parent name genotyped as B>   
        ...   
options:   
  -h, --help        show this help message and exit   
  -I , --input      VCF file which is the output of tasdosA.   
  -p1 , --parent1   Parent name genotyped as A.   
                    Use the name of vcf column in the input file of tasdosA.   
  -p2 , --parent2   Parent name genotyped as B.  
                    Use the name of vcf column in the input file of tasdosA.  
  --mindep          Minimum depth to genotype.  
                    Variants with depth lower than this  
                    will be genotyped as missing.   
  --hetero_chi      Threshold value of chi-square when genotyping as hetero.   
                    Default value is the threshold for p=0.05   
  --noise_level     When genotyping as homo, differrent variants below this ratio will be  accepted.  
  -v, --version     show program's version number and exit    


tasdosC    
usage: tasdosC -I <TSV file which is the output of tasdosB>    
        --parent_sample1 <Parent sample expected to be A>    
        --parent_sample2 <Parent sample expected to be B>    
        ...   
options:   
  -h, --help         show this help message and exit   
  -I , --input       TSV file which is the output of tasdosB.   
  --parent_sample1   Parent sample expected to be genotype A.   
                     This must be specified if parental lines are included in your samples.   
  --parent_sample2   Parent sample expected to be genotype B.  
                     This must be specified if parental lines are included in your samples.  
  --missing_rate     Markers with more missing than this  
                     value will be removed  
  --check_parents    Test the genotype of the parent line.   
                     If they are inconsistent with the predicted genotype, the marker will be removed.  
                     This is invalid if -p1 and -p2 are not specified.  
  --minor_freq       Threshold of minor allele frequency (MAF).  
                     Markers whose MAF are lower than this,  
                     they are removed.  
  -v, --version      show program's version number and exit  
    
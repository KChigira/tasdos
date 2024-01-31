import argparse
import sys
from tasdos.__init__ import __version__

class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'tasdosa':
            parser = self.tasdosa_options()
        elif self.program_name == 'tasdosb':
            parser = self.tasdosb_options()
        elif self.program_name == 'tasdosc':
            parser = self.tasdosc_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args
    
    def tasdosa_options(self):
        parser = argparse.ArgumentParser(description='tasdos version {}'.format(__version__),formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('tasdosA -I <Directory containing input FASTQ>\n'
                        '        -R <File of reference FASTA>\n'
                        '        -V <File of target VCF>\n'
                        '        ... \n')

        # set options
        parser.add_argument('-I', '--input',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory containing input FASTQ.\n'
                                  'This directory must contain only fastq file used in genotyping.\n'
                                  'gzip (fastq.gz) also supported.\n'),
                            metavar='')
        
        parser.add_argument('-R', '--ref',
                            action='store',
                            required=True,
                            type=str,
                            help='File of reference FASTA.',
                            metavar='')
        
        parser.add_argument('-V', '--vcf',
                            action='store',
                            required=True,
                            type=str,
                            help=('File of target VCF.\n'
                                  '(VCF made by mkselect is recommended.)\n'),
                            metavar='')
        
        parser.add_argument('--cpu',
                            action='store',
                            default=2,
                            type=int,
                            help=('Number of CPUs to use.\n'),
                            metavar='')
        
        parser.add_argument('--adapter',
                            action='store',
                            default='NONE',
                            choices=['NONE', 'NEXTERA', 'TRUSEQ', 'CUSTOM'],
                            help=('Adapter sequences used for trimming fastq.\n'
                                  'NONE means the input fastq has already trimmed.\n'
                                  'When CUSTOM designated, --adapterfile must be specified.\n'),
                            metavar='')
        
        parser.add_argument('--adapterfile',
                            action='store',
                            type=str,
                            help=('This is valid when --adapter = CUSTOM.\n'),
                            metavar='')
        
        parser.add_argument('--seqlen',
                            action='store',
                            default=150,
                            type=int,
                            help=('Sequence length of fastq.\n'),
                            metavar='')
        
        parser.add_argument('--minlen',
                            action='store',
                            default=60,
                            type=int,
                            help=('Ignore reads which are shorter than this value after trimming.\n'),
                            metavar='')
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def tasdosb_options(self):
        parser = argparse.ArgumentParser(description='tasdos version {}'.format(__version__),formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('tasdosB -I <VCF file which is the output of tasdosA>\n'
                        '        -p1 <Parent name genotyped as A>\n'
                        '        -p2 <Parent name genotyped as B>\n'
                        '        ... \n')

        # set options
        parser.add_argument('-I', '--input',
                            action='store',
                            required=True,
                            type=str,
                            help='VCF file which is the output of tasdosA.',
                            metavar='')
        
        parser.add_argument('-p1', '--parent1',
                            action='store',
                            required=True,
                            type=str,
                            help=('Parent name genotyped as A.\n'
                                  'Use the name of vcf column in the input file of tasdosA.\n'),
                            metavar='')
        
        parser.add_argument('-p2', '--parent2',
                            action='store',
                            required=True,
                            type=str,
                            help=('Parent name genotyped as B.\n'
                                  'Use the name of vcf column in the input file of tasdosA.\n'),
                            metavar='')
        
        parser.add_argument('--mindep',
                            action='store',
                            default=6,
                            type=int,
                            help=('Minimum depth to genotype.\n'
                                  'Variants with depth lower than this\n'
                                  'will be genotyped as missing.\n'),
                            metavar='')

        parser.add_argument('--hetero_chi',
                            action='store',
                            default=3.84,
                            type=float,
                            help=('Threshold value of chi-square when genotyping as hetero.\n'
                                  'Default value is the threshold for p=0.05\n'),
                            metavar='')
        
        parser.add_argument('--noise_level',
                            action='store',
                            default=0.1,
                            type=float,
                            help=('When genotyping as homo, differrent variants below this ratio will be  accepted.\n'),
                            metavar='')

        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def tasdosc_options(self):
        parser = argparse.ArgumentParser(description='tasdos version {}'.format(__version__),formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('tasdosC -I <TSV file which is the output of tasdosB>\n'
                        '        --parent_sample1 <Parent sample expected to be A>\n'
                        '        --parent_sample2 <Parent sample expected to be B>\n'
                        '        ... \n')

        # set options
        parser.add_argument('-I', '--input',
                            action='store',
                            required=True,
                            type=str,
                            help='TSV file which is the output of tasdosB.',
                            metavar='')
        
        parser.add_argument('--parent_sample1',
                            action='store',
                            default=None,
                            type=str,
                            help=('Parent sample expected to be genotype A.\n'
                                  'This must be specified if parental lines are included in your samples.\n'),
                            metavar='')
        
        parser.add_argument('--parent_sample2',
                            action='store',
                            default=None,
                            type=str,
                            help=('Parent sample expected to be genotype B.\n'
                                  'This must be specified if parental lines are included in your samples.\n'),
                            metavar='')
        
        parser.add_argument('--missing_rate',
                            action='store',
                            default=0.2,
                            type=float,
                            help=('Markers with more missing than this\n'
                                  'value will be removed\n'),
                            metavar='')

        parser.add_argument('--check_parents',
                            action='store_true',
                            help=('Test the genotype of the parent line.\n'
                                  'If they are inconsistent with the predicted genotype, the marker will be removed.\n'
                                  'This is invalid if -p1 and -p2 are not specified.\n'))

        parser.add_argument('--minor_freq',
                            action='store',
                            default=0,
                            type=float,
                            help=('Threshold of minor allele frequency (MAF).\n'
                                  'Markers whose MAF are lower than this,\n'
                                  'they are removed.\n'),
                            metavar='')


        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

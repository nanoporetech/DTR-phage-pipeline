import os,sys
from pathlib import Path
from snakemake import shell
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('raw_reads', help='FASTA containing the reads used for polishing', type=str)
    parser.add_argument('draft', help='FASTA containing the single read to serve as the draft genome', type=str)

    # Optional arguments
    parser.add_argument('-t', '--threads', help='Threads to use [4]', type=int, default=4)
    parser.add_argument('-n', '--iterations', help='Number of polishing iterations [3]', type=int, default=3)
    parser.add_argument('-q', '--min_q', help='Minimum Qscore to include in polishing [9]', type=int, default=9)
    parser.add_argument('-o', '--output', help='Output file containing polished genome [output.racon_3x.fasta]', type=str, default='output.racon_3x.fasta')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):

    final_output_ext = args.output.split(".")[-1]

    if os.path.getsize(args.draft)==0:
        # Create empty polished file
        shell("touch {}".format(args.output))
    else:
        
        final_output = args.output.replace('{}x'.format(args.iterations), '{nrepeat}x')

        for repeat in range(1,args.iterations+1):
            polished_contigs        = final_output.format(nrepeat=repeat)
            previous_polished_draft = final_output.format(nrepeat=repeat-1)
            alignment               = polished_contigs.replace(final_output_ext, '.aln.sam')
            if repeat == 1: 
                previous_polished_draft = args.draft

            shell('minimap2 -t {args.threads} -ax map-ont '
                '{previous_polished_draft} {args.raw_reads} > {alignment}; '
                      
                'racon ' 
                '--include-unpolished '
                '--quality-threshold={args.min_q} '
                '-t {args.threads} '
                '{args.raw_reads} '
                '{alignment} '
                '{previous_polished_draft} > {polished_contigs}')

if __name__=='__main__':
    args = parse_args()

    main(args)

import os
import subprocess
import pdb
import sys
from util import parent_dir, run_cmd

def aligner_tophat2(args):

    pdb.set_trace()

    sam_file = args[args.index('-o')+1]
    genome_file = args[args.index('-g')+1]

    genome_index = parent_dir(sam_file)+'/genome/'
    run_cmd('mkdir -p %s'%genome_index)
    genome_index = '%s/genome_index'%(genome_index)

    if '-r1' in args:
        readFiles = []
        readFiles.append(args[args.index('-r1')+1])
        readFiles.append(args[args.index('-r2')+1])
    elif '-r' in args:
        readFiles = []
        readFiles.append(args[args.index('-r')+1])

    if '-N' in args:
        numThread = int(args[args.index('-N')+1])
    else:
        numThread = 1

    cmd = 'bowtie2-build %s %s'%(genome_file, genome_index)

    run_cmd(cmd)

    cmd = 'tophat2 -p %d --no-sort-bam -o %s %s %s'%(numThread, parent_dir(sam_file), genome_index, ' '.join(readFiles))

    run_cmd(cmd)

    return

'''
usage:

-- tophat2

SE:

python aligner.py -a tophat2 -o sam_file -g genome_file -r readFile [-N num_Thread]

PE:

python aligner.py -a tophat2 -o sam_file -g genome_file -r1 readFile1 -r2 readFile2 [-N num_Thread]

-- (other aligners)


'''

if __name__ == '__main__':

    args = sys.argv

    choice = args[args.index('-a')+1]
    if choice == 'tophat2':
        aligner_tophat2(args)
from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import os
import subprocess
import pdb
import sys
from RefShannon.util import parent_dir, run_cmd
from RefShannon.dep_path import tool_paths

ROOT = os.path.dirname(__file__)

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

    cmd = '%s %s %s'%(tool_paths["bowtie2-build"], genome_file, genome_index)

    run_cmd(cmd)

    cmd = '%s -p %d --no-sort-bam -o %s %s %s'%(
        tool_paths["tophat2"], numThread, parent_dir(sam_file), genome_index, ' '.join(readFiles))

    run_cmd(cmd)

    return

def aligner_hisat2(args):

    genome_file = args[args.index('-g')+1]

    sam_file = args[args.index('-o')+1]
    sam_file_dir = parent_dir(sam_file)
    run_cmd('mkdir -p %s'%sam_file_dir)

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

    if '-sorted_bam' in args:
        sorted_bam = True
    else:
        sorted_bam = False

    if '-fasta' in args:
        read_type = '-f'
    else:
        read_type = '-q'
    # index
    index_hisat2_dir = parent_dir(genome_file)+'/index_hisat2/'
    run_cmd('mkdir -p %s'%index_hisat2_dir)
    index_hisat2 = '%s/index'%(index_hisat2_dir) # [genome_file_dir]/index_hisat2/index

    if os.path.exists(index_hisat2 + '.1.ht2') == False:
      cmd = '%s -f -p %d %s %s'%(
        tool_paths["hisat2-build"], numThread, genome_file, index_hisat2)
      run_cmd(cmd)

    # alignment
    if len(readFiles) == 1:
      cmd = '%s -p %d %s -x %s -U %s -S %s'%(
        tool_paths["hisat2"], numThread, read_type, index_hisat2, readFiles[0], sam_file)
    elif len(readFiles) == 2:
      cmd = '%s -p %d %s -x %s -1 %s -2 %s -S %s'%(
        tool_paths["hisat2"], numThread, read_type, index_hisat2, readFiles[0], readFiles[1], sam_file)        
    else:
      print('unexpected number of read files at aligner_hisat2')
      pdb.set_trace()
    run_cmd(cmd)

    # sam to bam
    if sorted_bam is True:
      bam_file = sam_file[:-3] + "bam"
      sorted_bam_file = sam_file[:-3] + "sorted.bam"

      cmd = '%s view -@ %d -bS %s > %s'%(
        tool_paths["samtools"], numThread, sam_file, bam_file)
      run_cmd(cmd)

      cmd = '%s sort -@ %d %s -o %s'%(
        tool_paths["samtools"], numThread, bam_file, sorted_bam_file)
      run_cmd(cmd)

      cmd = 'rm %s'%(sam_file)
      run_cmd(cmd)

    return

'''
usage:

-- tophat2

SE:

python aligner.py -a tophat2 -o sam_file -g genome_file -r readFile [-N num_Thread]

PE:

python aligner.py -a tophat2 -o sam_file -g genome_file -r1 readFile1 -r2 readFile2 [-N num_Thread]

-- hisat2

SE:

python aligner.py -a hisat2 -o sam_file -g genome_file -r readFile [-N num_Thread] [-fasta] [-sorted_bam]

PE:

python aligner.py -a hisat2 -o sam_file -g genome_file -r1 readFile1 -r2 readFile2 [-N num_Thread] [-fasta] [-sorted_bam]

-- (other aligners)


'''

if __name__ == '__main__':

    args = sys.argv

    choice = args[args.index('-a')+1]
    if choice == 'tophat2':
        aligner_tophat2(args)
    elif choice == 'hisat2':
        aligner_hisat2(args)
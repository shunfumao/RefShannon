
import os
import subprocess
import pdb
import sys
from util import parent_dir, run_cmd

'''
usage:

SE: 

python aligner_star.py -o sam_file -g genome_file -r read_file [-N num_Thread]

e.g.
python aligner_star.py -o /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/hits.sam 
                       -g /data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/genome/chr15.fa
                       -r /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reads_1.fa
                       -N 10

PE:

python aligner_star.py -o sam_file -g genome_file -r1 read_file1 -r2 read_file2 [-N num_Thread]

e.g.
python aligner_star.py -o /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/hits.sam 
                       -g /data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/genome/chr15.fa
                       -r1 /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reads_1.fa
                       -r2 /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reads_2.fa
                       -N 10
'''

def parse_args_aligner_star():

    args = sys.argv

    args_error = False 

    if '-o' in args:
        sam_file = args[args.index('-o')+1]
        #data_fld = '/'.join(sam_file.split('/')[:-1])+'/'
        data_fld = parent_dir(sam_file)+'/'
        sam_name = sam_file.split('/')[-1]
    else:
        args_error = True

    if '-g' in args:
        genomeFile = args[args.index('-g')+1]
        genomeDir = parent_dir(genomeFile)+'/'
    else:
        args_error = True

    if '-r' in args:
        readFile = args[args.index('-r')+1]
    elif '-r1' in args:
        readFile = args[args.index('-r1')+1]
        if '-r2' in args:
            readFile += ' ' + args[args.index('-r2')+1]
        else:
            args_error = True
    else:
        args_error = True

    if '-N' in args:
        num_Thread = int(args[args.index('-N')+1])
    else:
        num_Thread = 1

    return [args_error, data_fld, sam_name, genomeDir, genomeFile, readFile, num_Thread]

def aligner_star():

    [args_error, data_fld, sam_name, genomeDir, genomeFile, readFile, num_Thread]=parse_args_aligner_star()
    if args_error==True:
        print('aligner_star arguments error')
        return

    #1. mapping to reference
    cmd =   "STAR --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir+" "+\
            "--genomeFastaFiles "+genomeFile+" "+\
            "--runThreadN %d"%num_Thread
    run_cmd(cmd, shell=True)
    
    #2.
    runDir=data_fld+"/1pass/"
    cmd = "mkdir -p "+runDir
    run_cmd(cmd, shell=True)
    cmd =   "STAR --genomeDir "+genomeDir+" "\
            "--readFilesIn "+readFile+" "\
            "--outFileNamePrefix "+runDir+" "\
            "--runThreadN %d "%num_Thread+\
            "--outSAMstrandField intronMotif"
    run_cmd(cmd, shell=True)
    
    #3. mapping to reference (2pass)
    genomeDir_2pass=data_fld+"/genome_2pass/"
    cmd = "mkdir -p "+genomeDir_2pass
    run_cmd(cmd, shell=True)
    cmd =   "STAR --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir_2pass+" "+\
            "--genomeFastaFiles "+ genomeFile+" "+\
            "--sjdbFileChrStartEnd "+ runDir +"/SJ.out.tab "+\
            "--sjdbOverhang 99 --runThreadN %d"%num_Thread #75
    run_cmd(cmd, shell=True)
    
    #4.
    runDir_2pass=data_fld+"/2pass/"
    cmd = "mkdir -p "+runDir_2pass
    run_cmd(cmd, shell=True)
    cmd =   "STAR --genomeDir "+genomeDir_2pass+" "+\
            "--readFilesIn "+readFile+" "+\
            "--outFileNamePrefix "+runDir_2pass+" "+\
            "--runThreadN %d "%num_Thread+\
            "--outSAMstrandField intronMotif"
    run_cmd(cmd, shell=True)

    cmd = 'mv %s/Aligned.out.sam %s/%s'%(runDir_2pass, data_fld, sam_name)
    run_cmd(cmd, shell=True)

    cmd = 'rm -r %s'%(runDir)
    run_cmd(cmd, shell=True)

    cmd = 'rm -r %s'%(genomeDir_2pass)
    run_cmd(cmd, shell=True)

    cmd = 'rm -r %s'%(runDir_2pass)
    run_cmd(cmd, shell=True)

    '''
    itms = os.listdir(genomeDir)
    for itm in itms:
        if itm==genomeName: continue
        path = genomeDir + '/' + itm
        if os.path.isfile(path): os.remove(path)
        if os.path.isdir(path): shutil.rmtree(path)
    '''

    return

if __name__ == '__main__':
    aligner_star()

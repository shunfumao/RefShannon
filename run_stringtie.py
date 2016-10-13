from util import *
import sys, os, pdb

'''
dependencies:

samtools

usage:

#sam --> gtf and fasta
#sam_file should has (1) header info (2) xs tag to be used by stringtie

python run_stringtie.py -i sam_file -g genomeFile [-O out_dir] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]

#chrs_dir/chr_i/hits.sam --> out_dir/chr_i/algo_output/name_tag.gtf & name_tag.fasta

python run_stringtie.py -I chrs_dir -g multi_genomeFile [-O out_dir] [-chrs chr_a[,chr_b,...]] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]
python run_stringtie.py -I chrs_dir -G genomeDir [-O out_dir] [-chrs chr_a[,chr_b,...]] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]

'''

def do_stringtie_i_g(args):

    sam_file = args[args.index('-i')+1]
    genomeFile = args[args.index('-g')+1]
    files_to_clear = []

    if '-O' in args:
        res_dir = args[args.index('-O')+1]
    else:
        res_dir = parent_dir(sam_file)+'/algo_output/'
    run_cmd('mkdir -p %s'%res_dir)


    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-n' in args:
        name_tag = args[args.index('-n')+1]
    else:
        name_tag = 'stringtie'

    if '-addHead' in args:
        addHead = True
    else:
        addHead = False

    if '-clear' in args:
        clear = True
    else:
        clear = False

    if '-NoSeperatedLines':
        NoSeperatedLines = True #remove \n of sequences for each transcript
    else:
        NoSeperatedLines = False

    if addHead==True:

        bamFile = sam_file[:-4]+'.bam'
        bamFileSorted = sam_file[:-4]+'.sorted.bam'
        genomeFileIndex = genomeFile + '.fai'

        if os.path.exists(genomeFileIndex)==False:
            cmd = 'samtools faidx %s'%genomeFile
            run_cmd(cmd)
            files_to_clear.append(genomeFileIndex)

        if os.path.exists(bamFile)==False:
            cmd = 'samtools view -bt %s %s > %s -@ %d'%(genomeFileIndex, sam_file, bamFile, N_jobs) #header purpose
            run_cmd(cmd)
            files_to_clear.append(bamFile)

        if os.path.exists(bamFileSorted)==False:
            cmd = 'samtools sort -o %s %s -@ %d'%(bamFileSorted, bamFile, N_jobs)
            run_cmd(cmd)
            files_to_clear.append(bamFileSorted)

    gtfFile = res_dir + '/%s.gtf'%name_tag
    fastaFile = res_dir + '/%s.fasta'%name_tag

    if addHead==True:
        FileToUse = bamFileSorted
    else:
        FileToUse = sam_file

    cmd = 'stringtie %s -o %s -p %d -f 0.0 -c 0.001'%(FileToUse, gtfFile, N_jobs)
    run_cmd(cmd)

    cmd = 'gffread -w %s -g %s %s'%(fastaFile, genomeFile, gtfFile)
    run_cmd(cmd)

    if NoSeperatedLines==True:
        fastaFile2 = res_dir + '/%s2.fasta'%name_tag
        combineSeperateLines(fastaFile, fastaFile2)
        cmd = 'mv %s %s'%(fastaFile2, fastaFile)
        run_cmd(cmd)

    #clear files
    if clear==True:
        for f_path in files_to_clear:
            run_cmd('rm %s'%f_path)

    return

def do_stringtie_I_g(args):

    multi_genomeFile = args[args.index('-g')+1]
    
    tmpFolder = parent_dir(multi_genomeFile)+'/tmp_do_stringtie_I_g/'
    run_cmd('mkdir -p %s'%tmpFolder)

    cmd = 'python util.py --splitMultiFasta -i %s -O %s'%(multi_genomeFile, tmpFolder)
    run_cmd(cmd)

    args.append('-G')
    args.append(tmpFolder)
    do_stringtie_I_G(args)

    run_cmd('rm -r %s'%tmpFolder)

    return

def do_stringtie_I_G(args):

    chrs_dir = args[args.index('-I')+1]
    genomeDir = args[args.index('-G')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir)

    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-n' in args:
        name_tag = args[args.index('-n')+1]
    else:
        name_tag = 'stringtie'

    if '-addHead' in args:
        addHead_str = '-addHead'
    else:
        addHead_str = ''

    if '-clear' in args:
        clear_str = '-clear'
    else:
        clear_str = ''

    if '-NoSeperatedLines' in args:
        NoSeperatedLines_str = '-NoSeperatedLines'
    else:
        NoSeperatedLines_str = ''

    for target in target_list:
        sam_file = '%s/%s/hits.sam'%(chrs_dir, target)
        genomeFile = '%s/%s.fa'%(genomeDir, target)
        res_dir = '%s/%s/algo_output/'%(out_dir, target)
        target_args = '-i %s -g %s -O %s -N %d -n %s %s %s %s'% \
                      (sam_file, genomeFile, res_dir, N_jobs, name_tag, addHead_str, clear_str, NoSeperatedLines_str)
        do_stringtie_i_g(target_args.split())
        print('%s processed'%target)

    return

if __name__ == '__main__':

    args = sys.argv

    if '-i' in args and '-g' in args:
        do_stringtie_i_g(args)

    if '-I' in args and '-g' in args:
        do_stringtie_I_g(args)

    if '-I' in args and '-G' in args:
        do_stringtie_I_G(args)
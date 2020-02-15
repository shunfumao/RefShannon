from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

from RefShannon.util import *
import sys, os, pdb

from RefShannon.dep_path import extAssembler_paths, tool_paths

ROOT = os.path.dirname(__file__)

'''
dependencies:

samtools

usage:

#sam --> gtf and fasta
#sam_file should has (1) header info (2) xs tag to be used by ext assemblers (e.g. stringtie, cufflinks, scallop)
#other assemblers:
#TransComb: need tophat2 bam
#CLASS2: need path/to/run_class.pl; BAM format sorted by chromosome and position
#scallop: need sorted bam file; need XS tag or library type.
#strawberry: need sorted bam file
#ryuto: need to move its default gtf to another location (-o does not work)

python run_extAssembler.py [--assembler assemblerName] [--maxSens] (e.g. stringtie/cufflinks/TransComb/CLASS2/scallop. default stringtie) -i sam_file -g genomeFile [-O out_dir] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]

#chrs_dir/chr_i/hits.sam --> out_dir/chr_i/algo_output/name_tag.gtf & name_tag.fasta

python run_extAssembler.py [--assembler assemberName] [--maxSens]  -I chrs_dir -g multi_genomeFile [-O out_dir] [-chrs chr_a[,chr_b,...]] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]
python run_extAssembler.py [--assembler assemberName] [--maxSens]  -I chrs_dir -G genomeDir [-O out_dir] [-chrs chr_a[,chr_b,...]] [-N N_jobs] [-n name_tag] [-addHead] [-clear] [-NoSeperatedLines]

'''

def do_extAssembler_i_g(args):

    sam_file = args[args.index('-i')+1] #could be bam
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

    if '-NoSeperatedLines' in args:
        NoSeperatedLines = True #remove \n of sequences for each transcript
    else:
        NoSeperatedLines = False

    if '--assembler' in args:
        assembler = args[args.index('--assembler')+1]
    else:
        assembler = 'stringtie'

    if '--maxSens' in args:
        maxSens = True
    else:
        maxSens = False

    if addHead==True:

        bamFile = sam_file[:-4]+'.bam'
        bamFileSorted = sam_file[:-4]+'.sorted.bam'
        genomeFileIndex = genomeFile + '.fai'

        if os.path.exists(genomeFileIndex)==False:
            cmd = '%s faidx %s'%(tool_paths["samtools"], genomeFile)
            run_cmd(cmd)
            files_to_clear.append(genomeFileIndex)

        if os.path.exists(bamFile)==False:
            cmd = '%s view -bt %s %s > %s -@ %d'%(
              tool_paths["samtools"], genomeFileIndex, sam_file, bamFile, N_jobs) #header purpose
            run_cmd(cmd)
            files_to_clear.append(bamFile)

        if os.path.exists(bamFileSorted)==False:
            cmd = '%s sort -o %s %s -@ %d'%(
                tool_paths["samtools"], bamFileSorted, bamFile, N_jobs)
            run_cmd(cmd)
            files_to_clear.append(bamFileSorted)

    if assembler=='cufflinks':
        gtfFile = res_dir + '/transcripts.gtf'
        fastaFile = res_dir + '/%s.fasta'%name_tag
    elif assembler=='strawberry':
        gtfFile = res_dir + '/assembled_transcripts.gtf'
        fastaFile = res_dir + '/%s.fasta'%name_tag
    elif assembler=='trinity':
        gtfFile = None
        fastaFile_original = res_dir + '/Trinity-GG.fasta' 
        fastaFile = res_dir + '/%s.fasta'%name_tag
        # pdb.set_trace()
    else:
        gtfFile = res_dir + '/%s.gtf'%name_tag
        fastaFile = res_dir + '/%s.fasta'%name_tag

    if addHead==True:
        FileToUse = bamFileSorted
    else:
        FileToUse = sam_file

    if assembler=='stringtie':
        if maxSens==True:
            cmd = '%s %s -o %s -p %d -f 0.0 -c 0.001'%(
              extAssembler_paths["stringtie"], FileToUse, gtfFile, N_jobs)
        else:
            cmd = '%s %s -o %s -p %d '%(
              extAssembler_paths["stringtie"], FileToUse, gtfFile, N_jobs) # default setting
    elif assembler=='cufflinks':
        if maxSens==True:
            cmd = '%s -o %s -p %d -F 0.001 %s'%(
              extAssembler_paths["cufflinks"], parent_dir(gtfFile), N_jobs, FileToUse)
        else:
            cmd = '%s -o %s -p %d %s'%(
              extAssembler_paths["cufflinks"], parent_dir(gtfFile), N_jobs, FileToUse)
    elif assembler=='TransComb':
        if maxSens==True:
            print('TransComb max sens TBD')  #-f: default 0 -f 1 sens decreases
            pdb.set_trace()
        else:
            cmd = '%s -b %s -s unstranded -o %s -l 200'%(
              extAssembler_paths["TransComb"] ,FileToUse, gtfFile)
    elif assembler=='CLASS2':
        if maxSens==True:
            cmd = 'perl %s -a %s -o %s -p %d -F 0.0'%(
              extAssembler_paths["CLASS2"], FileToUse, gtfFile, N_jobs)        

        else:
            cmd = 'perl %s -a %s -o %s -p %d'%(
              extAssembler_paths["CLASS2"], FileToUse, gtfFile, N_jobs)   
        cmd = cmd + ' --wd %s/class_tmp/ --clean'%(parent_dir(gtfFile))     
    elif assembler=='scallop':
        cmd = '%s -i %s -o %s'%(
          extAssembler_paths["scallop"] ,FileToUse, gtfFile)
    elif assembler == 'strawberry':
        cmd = '%s -o %s --no-quant -p %d %s'%(
            extAssembler_paths["strawberry"], parent_dir(gtfFile), N_jobs, FileToUse)
    elif assembler == 'ryuto':
        ryuto_dir = parent_dir(gtfFile)
        cmd = 'mkdir -p %s'%ryuto_dir
        run_cmd(cmd)

        cmd = '%s index %s'%(tool_paths["samtools"], FileToUse)
        run_cmd(cmd)

        # tune
        cmd = '%s '%extAssembler_paths["ryuto"]
        if '--no-trimming' in args:
          cmd += '--no-trimming '
        if '--mean-filter' in args:
          cmd += '--mean-filter %s '%args[args.index('--mean-filter')+1]
        if '--score-filter' in args:
          cmd += '--score-filter %s '%args[args.index('--score-filter')+1]
        cmd += '%s '%FileToUse
        print(cmd)
        # pdb.set_trace()
    elif assembler == 'trinity':
        trinity_out_dir = parent_dir(fastaFile)
        cmd = 'mkdir -p %s'%trinity_out_dir
        run_cmd(cmd)

        cmd = '%s --genome_guided_bam %s '%(extAssembler_paths["trinity"], FileToUse)+\
              '--genome_guided_max_intron 10000 '+\
              '--max_memory 10G --CPU %d --output %s --full_cleanup'%(N_jobs, trinity_out_dir)

    # pdb.set_trace()
    print(cmd)
    run_cmd(cmd)
    # pdb.set_trace()

    if assembler == 'ryuto':
      cmd = 'mv transcripts.gtf %s'%gtfFile
      # pdb.set_trace()
      run_cmd(cmd)
    elif assembler == 'trinity':
      cmd = 'cp %s %s'%(fastaFile_original, fastaFile)
      run_cmd(cmd)
      # pdb.set_trace()

    if assembler != 'trinity':  # trinity excluded
      cmd = '%s -w %s -g %s %s'%(tool_paths["gffread"], fastaFile, genomeFile, gtfFile)
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

def do_extAssembler_I_g(args):

    multi_genomeFile = args[args.index('-g')+1]
    
    tmpFolder = parent_dir(multi_genomeFile)+'/tmp_do_extAssembler_I_g/'
    run_cmd('mkdir -p %s'%tmpFolder)

    cmd = 'python %s/util.py --splitMultiFasta -i %s -O %s'%(ROOT, multi_genomeFile, tmpFolder)
    run_cmd(cmd)

    args.append('-G')
    args.append(tmpFolder)
    do_extAssembler_I_G(args)

    run_cmd('rm -r %s'%tmpFolder)

    return

def do_extAssembler_I_G(args):

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

    if '--assembler' in args:
        assembler = args[args.index('--assembler')+1]
    else:
        assembler = 'stringtie'

    if '--maxSens' in args:
        maxSensStr = '--maxSens'
    else:
        maxSensStr = ''

    for target in target_list:
        sam_file = '%s/%s/hits.sam'%(chrs_dir, target)
        genomeFile = '%s/%s.fa'%(genomeDir, target)
        res_dir = '%s/%s/algo_output/'%(out_dir, target)
        target_args = '--assembler %s %s -i %s -g %s -O %s -N %d -n %s %s %s %s'% \
                      (assembler, maxSensStr, sam_file, genomeFile, res_dir, N_jobs, name_tag, addHead_str, clear_str, NoSeperatedLines_str)
        do_extAssembler_i_g(target_args.split())
        print('%s processed'%target)

    return

if __name__ == '__main__':

    args = sys.argv

    if '-i' in args and '-g' in args:
        do_extAssembler_i_g(args)

    if '-I' in args and '-g' in args:
        do_extAssembler_I_g(args)

    if '-I' in args and '-G' in args:
        do_extAssembler_I_G(args)

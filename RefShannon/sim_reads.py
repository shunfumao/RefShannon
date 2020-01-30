import sys, os, pdb, re, random
from RefShannon.util import *
import RefShannon.run_parallel_cmds

RNASeqReadSimulatorPath='RNASeqReadSimulator/' #folder path to store external software codes

def process_new_transcript(out_dir, name_prefix, target_files, tr_info):

    if tr_info[0]=='_': return

    chr_id = tr_info[0]
    tr_id = tr_info[1]
    strand = tr_info[2]
    exons = tr_info[3]
    exons.sort(key=lambda x:x[0])

    chrom_stt = exons[0][0]-1 #0-based
    chrom_stp = exons[-1][1] #0-based, exclusive
    blocks_cnt = len(exons)
    blocks_size = [(exons[i][1]-exons[i][0]+1) for i in range(blocks_cnt)]
    blocks_stt = [(exons[i][0]-1-chrom_stt) for i in range(blocks_cnt)]

    #pdb.set_trace()

    bed_line = ''
    bed_line += '%s\t'%chr_id
    bed_line += '%d\t'%chrom_stt
    bed_line += '%d\t'%chrom_stp
    bed_line += '%s\t'%tr_id
    bed_line += '0\t' #score
    bed_line += '%s\t'%strand
    bed_line += '%d\t'%chrom_stt #thick start
    bed_line += '%d\t'%chrom_stt #thick end
    bed_line += '0\t' #itemRGB
    bed_line += '%d\t'%blocks_cnt
    st = ','.join([str(i) for i in blocks_size])+','
    bed_line += '%s\t'%st
    st = ','.join([str(i) for i in blocks_stt])+','
    bed_line += '%s\n'%st

    if chr_id not in target_files:
        out_path = '%s/%s.target.%s.bed'%(out_dir, name_prefix, chr_id) #path/to/tr.target.chr_id.bed
        target_files[chr_id] = [out_path, open(out_path, 'w'), 0, 0]
    
    target_files[chr_id][1].write(bed_line)
    target_files[chr_id][2]+=1 #update trNum
    target_files[chr_id][3]+=sum(blocks_size) #update accTrLen
    #pdb.set_trace()

    return

'''
gtf format: http://uswest.ensembl.org/info/website/upload/gff.html (1-based, inclusive)
bed format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 (0-based, exclusive)
'''
def gtf2bed(args):

    in_gtf_file = args[args.index('-i')+1] #tr.gtf
    name_prefix = in_gtf_file.split('/')[-1][:-4] #tr

    out_dir = args[args.index('-O')+1]
    run_cmd('mkdir -p %s'%out_dir)

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = []

    num_lines = sum([1 for line in open(in_gtf_file)]); i=0; j=0; T=num_lines/100;
    print('%d lines at %s'%(num_lines, in_gtf_file))

    target_files = {} #key - chr_id, val - [file_path, bed file handler, trNum, accTrLen]
    tr_info = ['_', '_', '_', []] #chr, tr, strand, exon blocks (1-based, inclusive)

    with open(in_gtf_file, 'r') as f:

        for line in f:

            i+=1
            if i>T:
                i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of gtf processed'%j); sys.stdout.flush()
           
            if line[0]=='#': continue
 
            tokens = line.split()
            
            if target_list!=[] and tokens[0] not in target_list: continue            
            #try:
            if tokens[2]!='exon': continue
            #except:
            #pdb.set_trace()

            #pdb.set_trace()
            chr_id = tokens[0]
            tr_id = tokens[tokens.index('transcript_id')+1][1:-2]
            strand = tokens[6]
            exon_stt = int(tokens[3])
            exon_stp = int(tokens[4])

            if tr_info[0]==chr_id and tr_info[1]==tr_id and tr_info[2]==strand:
                tr_info[3].append([exon_stt, exon_stp])
                #pdb.set_trace()
            else:
                #pdb.set_trace()
                process_new_transcript(out_dir, name_prefix, target_files, tr_info)
                tr_info = [chr_id, tr_id, strand, []]
                tr_info[3].append([exon_stt, exon_stp])
                #pdb.set_trace()

        process_new_transcript(out_dir, name_prefix, target_files, tr_info)

    for chr_id, file_path_handler in target_files.items():
        file_path_handler[1].close()
        
        f1 = file_path_handler[0] #path/to/tr.target.chr_id.bed
        f2 = file_path_handler[0][:-4]+'.sorted.bed.tmp' #path/to/tr.target.chr_id.sorted.bed
        f3 = file_path_handler[0][:-4]+'.sorted.bed' #path/to/tr.target.chr_id.sorted.bed

        cmd = 'sort -k 1,1 -k 2,2n %s > %s'%(f1,f2) #sort
        run_cmd(cmd)

        comment_str = '# trNum %d avgTrLen %d'%(file_path_handler[2], file_path_handler[3]/file_path_handler[2])
        
        cmd = 'echo \'%s\' > %s'%(comment_str, f3)
        run_cmd(cmd)
        
        cmd = 'cat %s >> %s'%(f2, f3)
        run_cmd(cmd)

        cmd = 'rm %s %s'%(f1,f2) #remove unsorted
        run_cmd(cmd)

    print('')
    return


def sim_reads_g_b(args):

    chr_genome = args[args.index('-g')+1]
    chr_tr_bed = args[args.index('-b')+1]
    out_dir = args[args.index('-O')+1]

    readLength = int(args[args.index('-l')+1])

    if '-n' in args:
        numReads = int(args[args.index('-n')+1])
    elif '-c' in args:
        readCoverage = int(args[args.index('-c')+1])
        comment_str = open(chr_tr_bed, 'r').readline() # "# trNum (trNum) avgTrLen (avgTrLen)"
        comment_tokens = comment_str.split()
        trNum = int(comment_tokens[comment_tokens.index('trNum')+1])
        avgTrLen = int(comment_tokens[comment_tokens.index('avgTrLen')+1])
        numReads = int(float(readCoverage * trNum * avgTrLen) / readLength)
        print('%s (cov=%d) requires numReads=%s'%(chr_tr_bed, readCoverage, numReads))
        #pdb.set_trace()
    else:
        print('error: need -n or -c in sim_reads_g_b')
        return
    
    errRate = float(args[args.index('-r')+1])

    if '-p' in args:
        paired = True
        st = args[args.index('-p')+1]
        st = st.split(',')
        frag_mean = int(st[0])
        frag_sigma = int(st[1])
    else:
        paired = False
        frag_mean = -1
        frag_sigma = -1

    if '--removeIntermediate' in args:
        removeIntermediate = True
    else:
        removeIntermediate = False

    #pdb.set_trace()

    # exp
    out_dir_int = out_dir + '/intermediate/'
    run_cmd('mkdir -p %s'%out_dir_int)

    exp_path = out_dir_int + '/exp.txt'    
    cmd = 'python %s/genexplvprofile.py %s > %s'%(RNASeqReadSimulatorPath, chr_tr_bed, exp_path)
    run_cmd(cmd)

    # read bed
    read_bed_path = out_dir_int + '/reads.bed'
    if paired == True:
        cmd = 'python %s/gensimreads.py -e %s -n %d -l %d -p %d,%d --stranded %s > %s'% \
          (RNASeqReadSimulatorPath, exp_path, numReads, readLength, frag_mean, frag_sigma, chr_tr_bed, read_bed_path)
        run_cmd(cmd)
    else:
        cmd = 'python %s/gensimreads.py -e %s -n %d -l %d --stranded %s > %s'% \
          (RNASeqReadSimulatorPath, exp_path, numReads, readLength, chr_tr_bed, read_bed_path)
        run_cmd(cmd)

    # read fasta
    if paired == True:
        cmd = 'python %s/getseqfrombed.py -r %.2f -l %d %s %s | python %s/splitfasta.py -o %s/reads'% \
              (RNASeqReadSimulatorPath, errRate, readLength, read_bed_path, chr_genome, RNASeqReadSimulatorPath, out_dir)
        run_cmd(cmd)
    else:
        cmd = 'python %s/getseqfrombed.py -r %.2f -l %d %s %s > %s/reads.fa'% \
              (RNASeqReadSimulatorPath, errRate, readLength, read_bed_path, chr_genome, out_dir)
        run_cmd(cmd)

    # remove intermediate exp and read bed files
    if removeIntermediate == True:
        cmd = 'rm -r %s'%out_dir_int
        run_cmd(cmd)

    return

#genome_trBed_dic = {} #key - chr_i, val - [chr_genome_path, chr_tr_sorted_bed_path]
def prepare_genome_trBed_dic(genomesDir, trBedsDir, target_list):

    genome_trBed_dic = {}

    genomeFilesList = os.listdir(genomesDir)
    trBedFilesList = os.listdir(trBedsDir)

    for trBedFile in trBedFilesList:
        if len(trBedFile)>=4 and trBedFile[-4:]!='.bed': continue #non bed file

        itm = re.search(r'\.target\.(.+)\.sorted\.bed',trBedFile)
        if itm==None: continue #no chr info

        chr_i = itm.group(1)
        genome_trBed_dic[chr_i] = ['', trBedsDir+'/'+trBedFile]

    #pdb.set_trace()

    for genomeFile in genomeFilesList:
        if len(genomeFile)>=3 and genomeFile[-3:]!='.fa': continue #non fa file

        itm = re.search(r'(.+)\.fa',genomeFile)
        if itm==None: continue #no chr info

        chr_i = itm.group(1)
        if chr_i in genome_trBed_dic:
            genome_trBed_dic[chr_i][0] = genomesDir + '/' + genomeFile
        else:
            genome_trBed_dic[chr_i] = [genomesDir + '/' + genomeFile, '']

    #pdb.set_trace()

    #keep chr_i that has both genomeFile and trBed
    for key, val in genome_trBed_dic.items():
        if val[0]=='' or val[1]=='':
            del genome_trBed_dic[key]
 
    #pdb.set_trace()

    #keep chr_i only in target_list (if target_list nonempty)
    if target_list != []:
        for target in genome_trBed_dic.keys():
            if target not in target_list:
                del genome_trBed_dic[target]

    #pdb.set_trace()
    return genome_trBed_dic

def sim_reads_G_B(args):

    genomesDir = args[args.index('-G')+1]
    trBedsDir = args[args.index('-B')+1]

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = []

    #genome_trBed_dic = {} #key - chr_i, val - [chr_genome_path, chr_tr_sorted_bed_path]
    genome_trBed_dic = prepare_genome_trBed_dic(genomesDir, trBedsDir, target_list)
    genome_trBed_list = genome_trBed_dic.items()

    if '-S' in args:
        numSamples = int(args[args.index('-S')+1])
    else:
        numSamples = 1

    OutDir = args[args.index('-O')+1]

    if '-n' in args:
        numReads = int(args[args.index('-n')+1])
        readQuantityStr = '-n %d'%numReads
    elif '-c' in args:
        readCoverage = int(args[args.index('-c')+1])
        readQuantityStr = '-c %d'%readCoverage
    else:
        print('error: need -n or -c in sim_reads_G_B')
        return
    #pdb.set_trace()

    readLength = int(args[args.index('-l')+1])
    errRate = float(args[args.index('-r')+1])

    if '-p' in args:
        pairedStr = ' '.join(args[args.index('-p'):args.index('-p')+2])
        #pdb.set_trace()
    else:
        pairedStr = ''

    if '--removeIntermediate' in args:
        removeIntermediateStr = '--removeIntermediate'
    else:
        removeIntermediateStr = ''

    if '-N' in args:
        nJobs = int(args[args.index('-N')+1])
    else:
        nJobs = 1

    if nJobs==1:
        for j in range(numSamples):
            for i in range(len(genome_trBed_list)):
                curr_output_dir = '%s/sample_%d/%s/'% \
                                  (OutDir, j, genome_trBed_list[i][0])
                cmd = 'python sim_reads.py --simReads'
                cmd += ' -g %s -b %s -O %s %s -l %d -r %f %s %s'% \
                       (genome_trBed_list[i][1][0], \
                        genome_trBed_list[i][1][1], \
                        curr_output_dir, \
                        readQuantityStr, \
                        readLength, \
                        errRate, \
                        pairedStr, \
                        removeIntermediateStr)
                run_cmd(cmd)

    else:

        #pdb.set_trace()
        r_str = ''.join(random.choice('0123456789') for _ in range(10))
        jobs_dir = '%s/tmp_%s/'%(OutDir, r_str)
        run_cmd('mkdir -p %s'%jobs_dir)

        #pdb.set_trace()
        cnt = 0
        for j in range(numSamples): #samples
            for i in range(len(genome_trBed_list)): #chrs
                curr_output_dir = '%s/sample_%d/%s/'% \
                                  (OutDir, j, genome_trBed_list[i][0])
                cmd = 'python sim_reads.py --simReads'
                cmd += ' -g %s -b %s -O %s %s -l %d -r %f %s %s'% \
                       (genome_trBed_list[i][1][0], \
                        genome_trBed_list[i][1][1], \
                        curr_output_dir, \
                        readQuantityStr, \
                        readLength, \
                        errRate, \
                        pairedStr, \
                        removeIntermediateStr)

                ith_job_file = jobs_dir + '/%d.txt'%cnt
                with open(ith_job_file, 'w') as ijf:
                    ijf.write(cmd)
                cnt += 1

        #pdb.set_trace()
        #pstring = ' '.join([str(i) for i in range(cnt)])
        #cmd = 'time parallel --jobs %d python sim_reads.py --simReads -job %s/{}.txt ::: %s'% \
        #       (nJobs, jobs_dir, pstring)
        #run_cmd(cmd)

        cmds = []
        for i in range(cnt):
            cmd = 'python sim_reads.py --simReads -job %s/%d.txt'%(jobs_dir, i)
            cmds.append(cmd)
        run_parallel_cmds.run_cmds(cmds, nJobs)

        run_cmd('rm -r %s'%jobs_dir)
        #pdb.set_trace()

    return

# python sim_reads.py --simReads -job path/to/job_file
# load commands from path/to/job_file to run python sim_reads.py --simReads -g -b
def sim_reads_job(args):

    job_file = args[args.index('-job')+1]
    with open(job_file) as f:
        for line in f:
            cmd = line.strip()
            run_cmd(cmd)

    return


def poolChrs(args):

    sampleDir = args[args.index('-I')+1]
    OutDir = args[args.index('-O')+1]

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = []

    if '-paired' in args:
        paired = True
    else:
        paired = False

    if paired==True:
        dstFilePaths = [OutDir+'/reads_1.fa', OutDir+'/reads_2.fa']
        dstFiles = [open(dstFilePaths[0], 'w'), open(dstFilePaths[1], 'w')]
    else:
        dstFilePaths = [OutDir+'/reads.fa']
        dstFiles = [open(dstFilePaths[0], 'w')]

    for chr_i in os.listdir(sampleDir):
        if target_list!=[] and chr_i not in target_list:
            continue

        if paired==True:
            srcFilePaths = [sampleDir+'/'+chr_i+'/reads_1.fa', sampleDir+'/'+chr_i+'/reads_2.fa']
            if os.path.exists(srcFilePaths[0])==False or os.path.exists(srcFilePaths[1])==False:
                continue
        else:
            srcFilePaths = [sampleDir+'/'+chr_i+'/reads.fa']
            if os.path.exists(srcFilePaths[0])==False:
                continue

        nLines = sum([1 for line in open(srcFilePaths[0])]); 

        for iFile in range(len(dstFilePaths)):
            f = open(srcFilePaths[iFile],'r')
            i=0; j=0; T=nLines/100;
            for line in f:
                i+=1
                if i>=T:
                    i=0; j+=1; sys.stdout.write('\r');
                    sys.stdout.write('%s (%d-th reads): %d %% processed'%(chr_i, iFile, j)); sys.stdout.flush()
            
                if line[0]=='>':
                    st = '>'+'target_'+chr_i+'_'+line[1:]
                    dstFiles[iFile].write(st)
                else:
                    dstFiles[iFile].write(line)
            #pdb.set_trace()
            f.close()
            print('')

    for iFile in range(len(dstFiles)):
        dstFiles[iFile].close()

    print('')

    return

def poolSamples(args):

    InputDir = args[args.index('-I')+1]
    OutDir = args[args.index('-O')+1]

    if '-samples' in args:
        sample_str = args[args.index('-samples')+1]
        sample_list = [('sample_'+str(itm)) for itm in sample_str.split(',') if itm != '']
    else:
        sample_list = []

    if '-paired' in args:
        paired = True
    else:
        paired = False

    if paired==True:
        dstFilePaths = [OutDir+'/reads_1.fa', OutDir+'/reads_2.fa']        
    else:
        dstFilePaths = [OutDir+'/reads.fa']
    dstFiles = [open(df,'w') for df in dstFilePaths]

    for sampleName in os.listdir(InputDir):
        #pdb.set_trace()
        if sample_list!=[] and sampleName not in sample_list: continue

        if paired==True:
            srcFilePaths = [InputDir+'/'+sampleName+'/reads_1.fa', InputDir+'/'+sampleName+'/reads_2.fa']
        else:
            srcFilePaths = [InputDir+'/'+sampleName+'/reads.fa']

        toProcess = True
        for src_file in srcFilePaths:
            if os.path.exists(src_file)==False:
                toProcess = False; break;

        if toProcess == False: continue

        nLines = sum([1 for line in open(srcFilePaths[0])])

        for iFile in range(len(srcFilePaths)):
            with open(srcFilePaths[iFile], 'r') as f:
                i=0; j=0; T=nLines/100;
                for line in f:
                    i+=1
                    if i>=T:
                        i=0; j+=1; sys.stdout.write('\r');
                        sys.stdout.write('%s (%d-th reads): %d %% processed'% \
                                        (sampleName, iFile, j))
                        sys.stdout.flush()

                    if line[0]=='>':
                        st = '>'+sampleName+'_'+line[1:]
                        dstFiles[iFile].write(st)
                    else:
                        dstFiles[iFile].write(line)
                print('')

    for iFile in range(len(dstFiles)):
        dstFiles[iFile].close()

    print('')

    return

'''
usage:

# tr.gtf --> out_dir/tr.target.chr_a.sorted.bed, out_dir/tr.target.chr_b.sorted.bed, ...
#
# required naming:
# input: tr.gtf: 'exon' lines contain transcript_id "(tr_id_to_be_extracted)" 
# output: (tr).target.(chr_a).sorted.bed
# 
# for each output (tr).target.(chr_a).sorted.bed, a line of "# trNum (trNum) avgTrLen (avgTrLen)" will be added.
# this will be helpful to determine numReads to be generated if readCoverage is given
#

python sim_reads.py --gtf2bed -i path/to/tr.gtf -O out_dir [-chrs chr_a[,chr_b,...]]

# (chr_i).fa, (tr).target.(chr_i).sorted.bed
#                               --> out_dir/intermediate/exp.txt, reads.bed 
#                               --> out_dir/reads.fa (SE) or reads_1.fa + reads_2.fa (PE)
#                               (e.g. out_dir=path/to/sample_j/chr_i/)
# 
# (tr).target.(chr_i).sorted.bed needs 1st line to be "# trNum (trNum) avgTrLen (avgTrLen)" for '-c readCoverage'
#

python sim_reads.py --simReads -g genomeFile -b trBedFile -O out_dir
                               (-n numReads or -c readCoverage) -l readLength -r errRate
                               [-p m,d] [--removeIntermediate]

# genomesDir/(chr_i).fa,... + trBedsDir/(tr).target.(chr_i).sorted.bed,... restricted by -chrs list of chr_a,...
#                                [if no --removeIntermediate]
#                                --> outDir/sample_j/chr_i/intermediate/exp.txt + reads.bed
#                                --> outDir/sample_j/chr_i/reads.fa (SE) or reads_1.fa & reads_2.fa (PE)
#
# required naming for -G -B mode: (chr_i).fa, (tr).target.(chr_i).sorted.bed
# reads are strand specific (e.g. PE reads, rA represents transcript strand)
#
# ! nJobs > 1 to be checked
# 

python sim_reads.py --simReads -G genomesDir -B trBedsDir [-chrs chr_a[,chr_b,...]]
                               [-S numSamples]
                               -O outDir
                               (-n numReads or -c readCoverage) -l readLength -r errRate
                               [-p m,d] [--removeIntermediate]
                               [-N nJobs]

# load schedule/job_i.txt to run python sim_reads.py --simReads -g -b
#
# only to be used by --simReads -G -B

python sim_reads.py --simReads -job job_file

# pool reads from different chrs of the same sample
#
# sampleDir/chr_i/reads.fa (SE) or reads_1.fa & reads_2.fa (PE) (constrained by -chrs)
#      --> OutDir/reads.fa (SE) or reads_1.fa & reads_2.fa (PE)
#
# old_read_name --> target_(chr_i)_old_read_name
#

python sim_reads.py --poolChrs -I sampleDir -O OutDir [-chrs chr_a[,chr_b,...]] [-paired]

# pool reads (containing different chrs) from different samples
#
# InputDir/sample_j/reads.fa (SE) or reads_1.fa & reads_2.fa (PE) (constrained by -samples)
#      --> OutDir/reads.fa (SE) or reads_1.fa & reads_2.fa (PE)
#
# target_(chr_i)_old_read_name --> sample_(j)_target_(chr_i)_old_read_name 
# s_a, s_b, ... are sample index numbers (e.g. 0,1,2,...)
#

python sim_reads.py --poolSamples -I InputDir -O OutDir [-samples s_a[,s_b,...]] [-paired]

'''
if __name__ == '__main__':

    args = sys.argv

    if '--gtf2bed' in args:
        gtf2bed(args)
    
    elif '--simReads' in args:
        if '-g' in args and '-b' in args:
            sim_reads_g_b(args)
        elif '-G' in args and '-B' in args:
            sim_reads_G_B(args)
        elif '-job' in args:
            sim_reads_job(args)
    
    elif '--poolChrs' in args:
        poolChrs(args)

    elif '--poolSamples' in args:
        poolSamples(args)




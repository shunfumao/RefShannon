import sys, pdb, math
#from filter_FP.filter_FP import *
from util import run_cmd
import run_parallel_cmds

def process_oneThreshold(args):

    rec_fasta = args[args.index('-i')+1]
    fp_rec_file = args[args.index('-o')+1]
    depth_file = args[args.index('-d')+1]
    log_file = args[args.index('-l')+1]
    threshold = float(args[args.index('-t')+1])

    write_filtered_tr(depth_file, rec_fasta, fp_rec_file, log_file, threshold)
    print('%s and %s written (t=%f)'%(fp_rec_file, log_file, threshold))

    return

def process_multiThresholds(args):

    out_dir = args[args.index('-O')+1]
    run_cmd('mkdir -p %s'%out_dir)

    thr_str = args[args.index('-t')+1]
    thr_list = [float(itm) for itm in thr_str.split(',') if itm != '']

    fp_rec_files = []

    #pdb.set_trace()

    for thr in thr_list:

        args2 = []
        
        args2 += args[args.index('-i'):args.index('-i')+2]

        fp_rec_file = '%s/rec_t_%.2f.fasta'%(out_dir, thr)
        fp_rec_files.append((fp_rec_file, 'rec_t_%.2f'%thr))

        log_file = '%s/rec_t_%.2f.log'%(out_dir, thr)

        args2.append('-o'); args2.append(fp_rec_file)

        args2 += args[args.index('-d'):args.index('-d')+2]

        args2.append('-l'); args2.append(log_file)

        args2.append('-t'); args2.append(str(thr))

        process_oneThreshold(args2)

    if '--eval' in args:

        #pdb.set_trace()

        ref_file = args[args.index('-r')+1]

        out_dir2 = out_dir + '/eval/'
        run_cmd('mkdir -p %s'%out_dir2)      

        for fp_rec_file, name_tag in fp_rec_files:

            cmd = 'python refShannon.py --eval -i %s -r %s -O %s -n %s'%(fp_rec_file, ref_file, out_dir2, name_tag)
            run_cmd(cmd)
    return

'''
filter inPer file (blat rec onto ref) to get outPer

restrict ref tr by refScope (e.g. oracle set, abso path & in .fasta format) & rec tr by recScope

no restriction on ref or rec if refScope is none or recScope is none
'''
def filterBlatRes(inPer,outPer,refScope=None,recScope=None):

    refScope = getTrIDs(refScope)
    recScope = getTrIDs(recScope)

    nLines = sum([1 for l in open(inPer)]); T=nLines/100; p=0; q=0;
    #cnt = 0; cnt_max = 100000;

    fi = open(inPer, 'r')
    fo = open(outPer, 'w')

    print('filterBlatRes from "%s":'%(inPer))

    #pdb.set_trace()
    for line in fi:
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed'%(q)); sys.stdout.flush()

        #cnt += 1;
        #if cnt > cnt_max: break;

        tokens = line.split()

        refName = tokens[9]; 
        recName = tokens[13];

        if (refScope != [] and refName not in refScope) or (recScope != [] and recName not in recScope):
            continue

        fo.write(line)

    print('')

    fi.close()
    fo.close()

    return

def getTrIDs(fa_file):

    print('getTrIDs from "%s":'%(fa_file))

    if fa_file is None: return []

    trIDs = []
    
    nLines = sum([1 for l in open(fa_file)]); T=nLines/100; p=0; q=0;
    for line in open(fa_file,'r'):
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed'%(q)); sys.stdout.flush()
        if line[0]=='>':
            trIDs.append(line.split()[0][1:])
    print('')
    return trIDs

#python filter_FP_batch.py --filterBlatRes -i inPer -o outPer [-t refScope -tLab refTag] [-r recScope -rLab recTag]
def process_filterBlatRes(args):

    inPer = args[args.index('-i')+1]
    
    outPer = args[args.index('-o')+1] #assume includes info about refTag and recTag
    
    if '-t' in args:
        refScope = args[args.index('-t')+1]
        refTag = args[args.index('-tLab')+1]
    else:
        refScope = None
        refTag = ''

    if '-r' in args:
        recScope = args[args.index('-r')+1]
        recTag = args[args.index('-rLab')+1]
    else:
        recScope = None
        recTag = ''

    filterBlatRes(inPer,outPer,refScope,recScope)

    return

# inPerIndex.0 ... inPerIndex.(numPerFiles-1) filtered to out_dir/inPerIndex.0_ref(Tag)_rec(Tag)
#python filter_FP_batch.py --filterBlatResBatch -p nJobs -i inPerIndex -n numPerFiles -O out_dir [-t refScope -tLab refTag] [-r recScope -rLab recTag]
def process_filterBlatRes_batch(args):

    #filter even one inPerIndex (sub part of Per file) is slow
    #maybe we need to generate fp file seperatedly

    return

#line_start: 1-based of in_name file
def cut_file(in_name,out_name,line_start,line_end):
    run_cmd('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )
    return

#python filter_FP_batch --cutFile1Job -i in_name -O out_dir -c chunkIndex (0,1,2,...) -s splitSize
def cut_file_1job(args):

    #pdb.set_trace()
    in_name = args[args.index('-i')+1]
    out_dir = args[args.index('-O')+1]
    chunkIndex = int(args[args.index('-c')+1])
    splitSize = int(args[args.index('-s')+1])

    tokens = [i for i in in_name.split('/') if i != '']
    out_name = '%s/%s.%d'%(out_dir, tokens[-1], chunkIndex)
    line_start = chunkIndex*splitSize + 1
    line_end = (chunkIndex+1)*splitSize

    cut_file(in_name, out_name, line_start, line_end)

    print('cut_file_1job %d-th part done'%chunkIndex)

    return

#python filter_FP_batch.py --cutFileNJobs -i in_name -O out_dir -n nJobs -s numSplit
#cut in_name into out_dir/in_name.0 ... in_name.(numSplit-1) using nJobs
def cut_file_parallel(args):

    in_name = args[args.index('-i')+1]
    out_dir = args[args.index('-O')+1]
    nJobs = int(args[args.index('-n')+1])
    numSplit = int(args[args.index('-s')+1])

    #pdb.set_trace()

    nLines = sum([1 for l in open(in_name)])
    splitSize = int(math.ceil(float(nLines)/numSplit))

    #x = [str(i) for i in range(numSplit)]
    #pstring = ' '.join(x)
    #pdb.set_trace()
    #cmd = 'time parallel --jobs %d python filter_FP_batch.py --cutFile1Job -i %s -O %s -c {} -s %d ::: %s'%(nJobs, in_name, out_dir, splitSize, pstring)
    #run_cmd(cmd)

    cmds = []
    for i in range(numSplit):
        cmd = 'python filter_FP_batch.py --cutFile1Job -i %s -O %s -c %d -s %d'%(in_name, out_dir, i, splitSize)
        cmds.append(cmd)
    run_parallel_cmds.run_cmds(cmds, nJobs)


    return

'''
similar as in analyze_simData_ROC.py
'''
MIN_TR_LEN = 200 #trLen below MIN_TR_LEN are not considered for FP/REC judge; typically 200; -1: disable

FP_FRAC = 0.9
def isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2, fpThreshold=FP_FRAC):

    #if float(match1)/recLen < FP_FRAC and float(match2)/refTrLen2 < FP_FRAC: #orig
    if float(match1)/recLen < fpThreshold:
    #if float(match1)/recLen < FP_FRAC or float(match1)/refTrLen1 < FP_FRAC:
        return True 
    else:
        return False

REC_FRAC = 0.9
def isReconstructed(match, recLen, refLen, recThreshold=REC_FRAC):

    if match >= recThreshold * refLen: #orig
    #if match >= REC_FRAC * refLen and match >= REC_FRAC * recLen:
        return True
    else:
        return False

'''
usage:
python filter_FP_batch.py --calcSensWrapper --logFile lf
'''
def calcSensWrapper(args): #(logFile, mtl=MIN_TR_LEN):

    #pdb.set_trace()

    lf = args[args.index('--logFile')+1]

    [num_ref_recovered, sens_ratio]  = calcSens(lf)
    
    print('%d\t%s'%(num_ref_recovered, lf))
    return [num_ref_recovered, sens_ratio] 

'''
input: logFile
output: num_ref_tr_recovered, sens ratio

logFile format: unique T_ref, T_rec, max m, l_ref, l_rec
'''
def calcSens(logFile, mtl=MIN_TR_LEN, recThreshold=REC_FRAC):
    num_ref_recovered = 0
    num_ref = 0

    nLines=sum([1 for l in open(logFile,'r')]); T=nLines/100; p=0; q=0;
    for line in open(logFile,'r'):
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (calcSens)'%q); sys.stdout.flush()

        if line[0]=='#': continue #skip comments

        tokens = line.split()
        match = int(tokens[2])
        
        refName = tokens[0]#unique
        refLen = int(tokens[3])

        recName = tokens[1]
        recLen = int(tokens[4])
        
        if mtl>0 and refLen < mtl:
            continue
        else:
            num_ref +=1
            if isReconstructed(match, recLen, refLen, recThreshold):
                num_ref_recovered += 1
    print('\n{} out of {} reconstructed (MIN_TR_LEN={})'.format(num_ref_recovered, num_ref, mtl))
    return [num_ref_recovered, float(num_ref_recovered)/num_ref]

'''
input: fp_logFile
output: num_rec_tr_fp, fp ratio

if not None, update tr_labels={} key - tr_id, val - -1: skip 0: non-fp 1: fp

fp_logFile format: unique T_rec, max m1, l_ref1, t_ref1, m2 (max m/l_ref), l_ref2, t_ref2
'''
def calcFP(fp_logFile, mtl=MIN_TR_LEN, tr_labels=None, fpThreshold=FP_FRAC):
    num_rec_tr_fp = 0
    num_rec = 0

    nLines=sum([1 for l in open(fp_logFile,'r')]); T=nLines/100; p=0; q=0;
    for line in open(fp_logFile):

        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (calcFP)'%q); sys.stdout.flush()

        if line[0]=='#': continue

        fields = line.split()

        if len(fields)<9:#rec tr not covered by any ref tr
            #pdb.set_trace()
            recTrName = fields[0]#unique
            recLen = int(fields[3])
            if mtl>0 and recLen < mtl:
                if tr_labels is not None: tr_labels[recTrName]=-1
                continue
            else:
                if tr_labels is not None: tr_labels[recTrName]=1
                num_rec += 1
                num_rec_tr_fp += 1            
                continue
        
        recTrName = fields[0]#unique
        recLen = int(fields[3])

        match1 = int(fields[1]) # max_t m_t,r
        refTrLen1 = int(fields[2])
        refTrName1 = fields[4]

        match2 = int(fields[5]) # max_t m_t,r/len(t)
        refTrLen2 = int(fields[6])
        refTrName2 = fields[8]

        if mtl>0 and recLen < mtl:
            if tr_labels is not None: tr_labels[recTrName]=-1
            continue
        else:
            num_rec +=1
            if isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2, fpThreshold):
                if tr_labels is not None: tr_labels[recTrName]=1
                num_rec_tr_fp+= 1
            else:
                if tr_labels is not None: tr_labels[recTrName]=0

    print('\n{} out of {} false positives (MIN_TR_LEN={})'.format(num_rec_tr_fp, num_rec, mtl))

    return [num_rec_tr_fp, float(num_rec_tr_fp)/num_rec]

'''
python filter_FP_batch.py --eval1Job -t Tref -r Trec -O out_dir
'''
def eval_1Job(args): #Tref & Trec --> per, log, fplog --> sens & fp

    #pdb.set_trace()

    Tref = args[args.index('-t')+1]
    Trec = args[args.index('-r')+1]
    out_dir = args[args.index('-O')+1]

    #out_dir1 = out_dir + '/tmp/'; run_cmd('mkdir -p %s'%out_dir1)
    #out_dir2 = out_dir + '/res/'; run_cmd('mkdir -p %s'%out_dir2)

    resNameStem = [i for i in Trec.strip().split('/') if i!='']
    resNameStem = resNameStem[-1][:-6] #skip .fasta

    #pdb.set_trace()

    perFile = '%s/%s_per.txt'%(out_dir, resNameStem)
    logFile = '%s/%s_log.txt'%(out_dir, resNameStem)
    fpLogFile = '%s/%s_fp_log.txt'%(out_dir, resNameStem)
    resFile = '%s/%s_res.txt'%(out_dir, resNameStem) 

    #gen log & per
    cmd = 'python analyze_simData_ROC.py --genPerLog --ref %s '%Tref + \
                                          '--rec %s '%Trec + \
                                          '--per %s '%perFile + \
                                          '--log %s'%logFile
    run_cmd(cmd)


    #gen fpLog
    cmd = 'python analyze_simData_ROC.py --genFpLog --rec %s '%Trec + \
                                         '--per %s '%perFile + \
                                         '-o   %s'%fpLogFile
    run_cmd(cmd)

    #cal sens
    [num_ref_recovered, ref_recovered_ratio] = calcSens(logFile)

    #cal fp
    [num_rec_fp, fp_ratio] = calcFP(fpLogFile)

    with open(resFile, 'w') as f:
        st = ''
        st += '#ref file: %s\n'%Tref
        st += '#rec file: %s\n'%Trec
        st += '#per file: %s\n'%perFile
        st += '#log file: %s\n'%logFile
        st += '#fp log file: %s\n'%fpLogFile
        st += '#num_ref_recovered\tref_recovered_ratio\tnum_rec_fp\tfp_ratio\n'
        st += '%d\t%f\t%d\t%f\n'%(num_ref_recovered, ref_recovered_ratio, num_rec_fp, fp_ratio)
        f.write(st)

    return

'''
analysis ROC can be done by analyze_simData_ROC or filter_FP_batch
- for analyze_simData_ROC, rec transcripts are filtered by abundance (Kallisto)
- for filter_FP_batch, rec transcripts are filtered by read coverage on rec transcripts (hisat2 --> depth)

usage:

python filter_FP_batch.py --oneThreshold -i rec_fasta -o fp_rec_file -d depth_file -l log_fie -t threshold

###if a transcript is re-covered by over 90%, output it

python filter_FP_batch.py --multiThresholds -i rec_fasta -O out_dir -d depth_file -t [threshold1,[threshold2,...]] [--eval -r ref_file]

###filtered fp files are rec_t_thr.fasta (and .log) stored at out_dir,
###eval files rec_t_thr_log.txt (and _per.txt) stored at out_dir/eval/ 

python filter_FP_batch.py --filterBlatRes -i inPer -o outPer [-t refScope -tLab refTag] [-r recScope -rLab recTag]

python filter_FP_batch --cutFile1Job -i in_name -O out_dir -c chunkIndex (0,1,2,...) -s splitSize
- called by --cutFileNJobs

python filter_FP_batch.py --cutFileNJobs -i in_name -O out_dir -n nJobs -s numSplit


'''
if __name__ == '__main__':

    args = sys.argv

    #filter fasta
    if '--oneThreshold' in args:
        process_oneThreshold(args)
    elif '--multiThresholds' in args:
        process_multiThresholds(args)

    #filter per
    elif '--filterBlatRes' in args:
        process_filterBlatRes(args)

    #cut file (e.g. per)
    elif '--cutFile1Job' in args:
        cut_file_1job(args)
    elif '--cutFileNJobs' in args:
        cut_file_parallel(args)

    elif '--eval1Job' in args:
        eval_1Job(args) #Tref & Trec --> per, log, fplog --> sens & fp
    elif '--calcSensWrapper' in args:
        calcSensWrapper(args)
    else:
        print('unexpected mode')
        pdb.set_trace() 


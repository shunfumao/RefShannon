import sys, tester, random
from util import *

def genPerLog(args):

    refFa = args[args.index('--ref')+1] #2
    recFa = args[args.index('--rec')+1] #1
    perFile = args[args.index('--per')+1] #3

    print('genPerLog: to generate perFile')

    run_cmd('python blat.py {} {} {}'.format(recFa, refFa, perFile))

    logFile = args[args.index('--log')+1]

    print('genPerLog: to generate logFile')

    tester.analyzer_blat_noExp(perFile, logFile, '', 0) 

    return

def genFpLog(args):

    recFa = args[args.index('--rec')+1]
    perFile = args[args.index('--per')+1]
    FpLog = args[args.index('-o')+1]

    #pdb.set_trace() # check with sim 0912
    tester.false_positive(recFa,perFile,FpLog)

    return

def quantByKallisto(args):

    out_prefix = args[args.index('-o')+1]
    recFa = args[args.index('--rec')+1]

    if '-r2' in args:
        reads = args[args.index('-r1')+1]
        reads += ' '+args[args.index('-r2')+1]
        SE_prefix = ' '
    else:
        reads = args[args.index('-r1')+1]
        SE_prefix = ' --single '

    if '-t' in args:
        nThreads = int(args[args.index('-t')+1])
    else:
        nThreads = 1

    cmd = 'kallisto index -i %s_kal.index %s --make-unique'%(out_prefix, recFa)
    run_cmd(cmd)

    #pdb.set_trace()

    if SE_prefix != ' ':
        cmd = 'kallisto quant -i %s_kal.index -o %s_kal_out -t %d %s -l 200 -s 20 %s '% \
          (out_prefix, out_prefix, nThreads, SE_prefix, reads)
    else:
        cmd = 'kallisto quant -i %s_kal.index -o %s_kal_out -t %d %s %s'% \
          (out_prefix, out_prefix, nThreads, SE_prefix, reads)

    run_cmd(cmd)        

    return

MIN_TR_LEN = 200 #trLen below MIN_TR_LEN are not considered for FP/REC judge

FP_FRAC = 0.9
def isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2):

    #if float(match1)/recLen < FP_FRAC and float(match2)/refTrLen2 < FP_FRAC: #orig
    if float(match1)/recLen < FP_FRAC:
    #if float(match1)/recLen < FP_FRAC or float(match1)/refTrLen1 < FP_FRAC:
        return True 
    else:
        return False

REC_FRAC = 0.9
def isReconstructed(match, recLen, refLen):

    if match >= REC_FRAC * refLen: #orig
    #if match >= REC_FRAC * refLen and match >= REC_FRAC * recLen:
        return True
    else:
        return False

def calc_fp(fpLogFile):

    false_positives = set()

    total = 0
    with open(fpLogFile) as f:
        #pdb.set_trace()
        for line in f:

            fields = line.split()

            if len(fields)<9:
                if len(fields)==8:
                    pdb.set_trace() #fpLogFile has refTrName1 and refTrName2 missing
                continue
            
            recTrName = fields[0]
            recLen = int(fields[3])

            match1 = int(fields[1]) # max_t m_t,r
            refTrLen1 = int(fields[2])
            refTrName1 = fields[4]

            match2 = int(fields[5]) # max_t m_t,r/len(t)
            refTrLen2 = int(fields[6])
            refTrName2 = fields[8]

            if recLen < MIN_TR_LEN: continue

            total += 1
            #if (float(fields[1])<FP_FRAC*float(fields[3])) and (float(fields[5])<FP_FRAC*float(fields[6])):
            if isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2):
                false_positives.add(recTrName)

    print('{} false positives of {}'.format(len(false_positives), total))
    
    return false_positives

def calc_rec2(rec_log, good_names):
    rec_set = set()
    tot = 0
    for line in open(rec_log, 'r'):
        if line[0]=='#': continue

        tokens = line.split()
        match = int(tokens[2])
        
        refName = tokens[0]
        refLen = int(tokens[3])

        recName = tokens[1]
        recLen = int(tokens[4])
        if recName not in good_names: continue

        if recLen < MIN_TR_LEN: continue
        tot +=1

        if isReconstructed(match, recLen, refLen):
            rec_set.add(refName) # rec_set is a set, we need to add refName (a tr may be covered by more than 1 rec tr)
    print('{} reconstructed'.format(len(rec_set)))
    return len(rec_set)

def calc_rec(rec_per, good_names):
    rec_set = set()
    tot = 0

    nLines=sum([1 for l in open(rec_per,'r')]); T=nLines/100; p=0; q=0;
    for line in open(rec_per,'r'):
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (perfile)'%q); sys.stdout.flush()

        tokens = line.split()
        match = int(tokens[0]);

        refName = tokens[9]; #need to select for the tokens[9][29:]  when having larger prefix
        refLen = int(tokens[10])

        recName = tokens[13];
        recLen = int(tokens[14]);
        if recName not in good_names: continue

        if recLen < MIN_TR_LEN: continue
        tot +=1

        if isReconstructed(match, recLen, refLen):
            rec_set.add(refName) # rec_set is a set, we need to add refName (a tr may be covered by more than 1 rec tr)
    print('{} reconstructed'.format(len(rec_set)))
    return len(rec_set)

def calculate_curve(cutoff, weights, fp_set, refFile, format): #rec_per --> rec_loc
    #pdb.set_trace()
    weights = weights[int(0.1*cutoff*len(weights)):]
    good_names = set([x[0] for x in weights])
    #pdb.set_trace()
    
    no_fp = len(fp_set.intersection(good_names))
    fr_fp = float((no_fp))/len(good_names)

    if format=='-p':
        no_rec = calc_rec(refFile, good_names) #per file; len(rec_set.intersection(good_names))
    elif format=='-l':
        no_rec = calc_rec2(refFile, good_names) #log file

    return [no_rec,fr_fp]

#from code_sim_0912a/fp_analysis/
def roc_curve(ab_file, fp_file, rec_file, format):
    weights = []
    with open(ab_file) as f:
        f.readline()
        for line in f:
            name, l1, _, ec, weight = line.split()
            weights.append((name, float(ec)/float(l1)*100)) #expected count at kallisto abundance output, 100 is scaling factor, does not matter

    weights.sort(key = lambda x: x[1])

    fp_set = calc_fp(fp_file)
    #pdb.set_trace()

    #rec_set = calc_rec(rec_per)
    #pdb.set_trace()

    rec_fp = []
    for cutoff in range(10):
        [rec,fp] = calculate_curve(cutoff, weights, fp_set, rec_file, format)
        rec_fp.append([float(cutoff)/10, rec,fp])
    #print rec_fp
    return rec_fp

def calcROC(args):

    rocFile = args[args.index('-o')+1]
    abFile = args[args.index('-a')+1]
    fpLogFile = args[args.index('-f')+1]
    if '-p' in args:
        recFile = args[args.index('-p')+1]
        rec_fp = roc_curve(abFile, fpLogFile, recFile, '-p') #list of (cutoff ratio, rec num, fp ratio)
    elif '-l' in args:
        recFile = args[args.index('-l')+1]
        rec_fp = roc_curve(abFile, fpLogFile, recFile, '-l') #list of (cutoff ratio, rec num, fp ratio)
    
    with open(rocFile, 'w') as f:
        st = '#%s\n#%s\n#%s\n'%(abFile, recFile, fpLogFile)
        f.write(st)

        st = '#cutoff_ratio\trec_num\tfp_ratio\n'
        f.write(st)

        for cutoff_ratio, rec_num, fp_ratio in rec_fp:
            st = '%.5f\t%d\t%.5f\n'%(cutoff_ratio, rec_num, fp_ratio) 
            f.write(st)

    return

def calcROC_1Point(args):

    cutoff = int(args[args.index('-c')+1])
    tempRocFile = (args[args.index('-o')+1])
    ab_file = (args[args.index('-a')+1])
    if '-p' in args:
        recFile = (args[args.index('-p')+1])
        f_format = '-p'
    elif '-l' in args:
        recFile = (args[args.index('-l')+1])
        f_format = '-l'

    fpLogFile = (args[args.index('-f')+1])

    weights = []
    with open(ab_file) as f:
        f.readline()
        for line in f:
            name, l1, _, ec, weight = line.split()
            weights.append((name, float(ec)/float(l1)*100)) #expected count at kallisto abundance output, 100 is scaling factor, does not matter
    weights.sort(key = lambda x: x[1])

    fp_set = calc_fp(fpLogFile)

    [rec,fp] = calculate_curve(cutoff, weights, fp_set, recFile, f_format)

    with open(tempRocFile, 'w') as f:
        st = '%.5f\t%d\t%.5f\n'%(float(cutoff)/10, rec, fp)
        f.write(st)

    return

def calcROC_Parallel(args):

    rocFile = args[args.index('-o')+1]
    abFile = args[args.index('-a')+1]
    fpLogFile = args[args.index('-f')+1]
    if '-p' in args:
        recFile = args[args.index('-p')+1]
        f_format = '-p'
    elif '-l' in args:
        recFile = args[args.index('-l')+1]
        f_format = '-l'

    nJobs = int(args[args.index('-N')+1])

    r_str = ''.join(random.choice('0123456789') for _ in range(10))
    tmpFld = parent_dir(rocFile) + '/temp_%s/'%r_str
    run_cmd('mkdir -p %s'%tmpFld)

    if 0:
        pdb.set_trace()
        for cutoff in range(10):
            tmpFile = '%s/%d.txt'%(tmpFld, cutoff)
            cmd = 'python analyze_simData_ROC.py --calcROC_1Point '+ \
                  '-c %d '%cutoff+ \
                  '-o %s '%tmpFile+ \
                  '-a %s '%abFile+ \
                  '%s %s '%(f_format, recFile)+ \
                  '-f %s '%fpLogFile
            run_cmd(cmd)
    else:
        pstring = ' '.join([str(i) for i in range(10)])
        cmd = 'time parallel --no-notice --jobs %d '%nJobs + \
              'python analyze_simData_ROC.py --calcROC_1Point '+ \
                  '-c {} '+ \
                  '-o %s/{}.txt '%tmpFld+ \
                  '-a %s '%abFile+ \
                  '%s %s '%(f_format, recFile)+ \
                  '-f %s '%fpLogFile+ \
                  '::: %s'%pstring
        run_cmd(cmd)

    #merge
    rec_fp = []
    tmpRocFileList = os.listdir(tmpFld)
    for rF in tmpRocFileList:
        with open(tmpFld + '/%s'%rF, 'r') as f:
            fields = f.readline().split()
            cutoff_ratio = float(fields[0])
            rec_num = int(fields[1])
            fp_ratio = float(fields[2])
            rec_fp.append([cutoff_ratio, rec_num, fp_ratio])
    rec_fp.sort(key=lambda x:x[0])

    #clear
    run_cmd('rm -r %s'%tmpFld)
    
    with open(rocFile, 'w') as f:
        st = '#%s\n#%s\n#%s\n'%(abFile, recFile, fpLogFile)
        f.write(st)

        st = '#cutoff_ratio\trec_num\tfp_ratio\n'
        f.write(st)

        for cutoff_ratio, rec_num, fp_ratio in rec_fp:
            st = '%.5f\t%d\t%.5f\n'%(cutoff_ratio, rec_num, fp_ratio) 
            f.write(st)

    return

#modified from sim_0912a/fp_analysis plot_roc
def plotROC(data):

    import numpy as np
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    markers = ['o', '>', '<', '+', '*']

    for i in range(len(data)):
        label, roc = data[i]
        marker = markers[i%len(markers)]
        x = [i[1] for i in roc] #rec num
        y = [i[2] for i in roc] #fp
        plt.plot(x, y, marker=marker, label=label)
    
    plt.legend(loc='upper center')
    ax.set_xlabel('number of recovered transcripts')
    ax.set_ylabel('false positive rate')

    plt.show()
    
    return

def getROCContents(ith_file):

    res = []

    with open(ith_file, 'r') as f:
        for line in f:
            if line[0]=='#': continue
            fields = line.split()
            cutoff = float(fields[0])
            rec = int(fields[1])
            fp = float(fields[2])
            res.append([cutoff, rec, fp])

    return res

def drawROC(args):

    data = [] #[label_i, list of [cutoff, rec, fp]]

    cnt = 1
    while '-i%d'%cnt in args:
        ith_file = args[args.index('-i%d'%cnt)+1]
        ith_label = args[args.index('-n%d'%cnt)+1]
        ith_content = getROCContents(ith_file)
        data.append([ith_label, ith_content])
        cnt += 1

    plotROC(data)

    return

'''
Usage:

python analyze_simData_ROC.py --genPerLog --ref path/to/refFa \
                                          --rec path/to/recFa \
                                          --per path/to/perFile \
                                          --log path/to/logFile

python analyze_simData_ROC.py --genFpLog --rec path/to/recFa \
                                         --per path/to/perFile \
                                         -o   path/to/FpLog

python analyze_simData_ROC.py --kallisto -o path/to/prefix \
                                         --rec path/to/recFa \
                                         -r1 path/to/reads_1.fa \
                                         [-r2 path/to/reads_2.fa] \
                                         [-t nThreads] \

python analyze_simData_calcROC.py --calcROC -o path/to/rocFile \
                                            -a path/to/abFile \
                                            [-p path/to/recPerFile or -l path/to/recLogFile] \
                                            -f path/to/fpLogFile 
                                            [-N nJobs]

python analyze_simData_calcROC.py --calcROC_1Point -c cutoff \
                                                   -o path/to/tempRocFile \
                                                   -a path/to/abFile \
                                                   [-p path/to/recPerFile or -l path/to/recLogFile] \
                                                   -f path/to/fpLogFile 

python analyze_simData_calcROC.py --drawROC -i1 rocFile1 -n1 label1 \
                                            [-i2 rocFile2 -n2 label2 [-i3 rocFile3 -n3 label3 ...]] #N/A on server

'''

if __name__ == '__main__':

    args = sys.argv

    if '--genPerLog' in args:
        genPerLog(args)
    elif '--genFpLog' in args:
        genFpLog(args)
    elif '--kallisto' in args:
        quantByKallisto(args)
    elif '--calcROC' in args:
        if '-N' in args:
            calcROC_Parallel(args)
        else:
            calcROC(args)
    elif '--calcROC_1Point' in args:
        calcROC_1Point(args)
    elif '--drawROC' in args:
        drawROC(args)



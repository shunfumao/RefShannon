from util import *
import pdb, subprocess

#tune FLOW_RATIO_THRESHOLD to see sens/fp change of refShannon

#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSimChr15/'
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySimChr15/'

sg_dir = '%s/refShannon/intermediate/'%dataDir
genomeFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/chr15.fa'
Tref = '%s/reference.fasta'%dataDir
rootResDir = '%s/refShannon_test/'%dataDir
rootResFile = '%s/res.txt'%rootResDir

nJobs = 5
#T = ['0.01', '0.02', '0.03', '0.04'] #['0.05','0.15', '0.2', '0.25', '0.3']#['0.1', '0.09', '0.08', '0.07', '0.06', '0.05']
#T = ['0.4', '0.5', '0.6', '0.7', '0.8', '0.9']
#T = ['0.92', '0.94', '0.96', '0.98']
T = ['0.0', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '0.95']
def gen_res():
    for t in T:
        resDir = '%s/%s/'%(rootResDir, t)
        run_cmd('mkdir -p %s'%resDir)

        #pdb.set_trace()

        #sg --> transcripts
        cmd = 'python refShannon.py --sf -I %s '%sg_dir+ \
              '-g %s '%genomeFile+ \
              '-O %s '%resDir+ \
              '-N %d '%nJobs+ \
              '-target chr15'
        run_cmd(cmd)

        #eval
        Trec = '%s/reconstructed.fasta'%(resDir)

        cmd = 'python filter_FP_batch.py --eval1Job ' + \
              '-t %s '%Tref + \
              '-r %s '%Trec + \
              '-O %s'%resDir
        run_cmd(cmd)

        print('%s done'%resDir)
    return

#extract res from gen_res()
def check_res():
    results = {}
    with open(rootResFile, 'w') as f:
        for t in T:
            resFile = '%s/%s/reconstructed_res.txt'%(rootResDir, t)
            try:
                res = subprocess.check_output('tail %s'%resFile, shell=True)
                res = res.strip().split()[-4:]
                results[t] = [int(res[0]), float(res[3])] #sens, fp
                f.write('t\t%s\tsens\t%d\tfp\t%f\n'%(t, int(res[0]), float(res[3])))
            except:
                results[t] = []
                f.write('t\t%s\n'%t)
    print(results)
    return

#to verify if tr_rec<200 (more likely to be non-fp) will affect fp
def check_trec_fp():

    from filter_FP_batch import isFalsePositive

    fp_logFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon_test/0.1/reconstructed_fp_log.txt'

    num_rec_tr_fp = 0
    num_rec = 0

    nLines=sum([1 for l in open(fp_logFile,'r')]); T=nLines/100; p=0; q=0;
    trec_ge200_fp = []; trec_ge200_nonfp = []
    trec_lt200_fp = []; trec_lt200_nonfp = []
    for line in open(fp_logFile):

        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (calcFP)'%q); sys.stdout.flush()

        if line[0]=='#': continue

        fields = line.split()

        if len(fields)<9:#rec tr not covered by any ref tr
            #pdb.set_trace()
            recLen = int(fields[3])
            if recLen < 200:
                trec_lt200_fp.append(fields[0])
                continue
            else:
                trec_ge200_fp.append(fields[0])            
                continue
        
        recTrName = fields[0]#unique
        recLen = int(fields[3])

        match1 = int(fields[1]) # max_t m_t,r
        refTrLen1 = int(fields[2])
        refTrName1 = fields[4]

        match2 = int(fields[5]) # max_t m_t,r/len(t)
        refTrLen2 = int(fields[6])
        refTrName2 = fields[8]

        if recLen < 200:
            if isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2):
                trec_lt200_fp.append(fields[0])
            else:
                trec_lt200_nonfp.append(fields[0])
        else:
            if isFalsePositive(recLen, match1, refTrLen1, match2, refTrLen2):
                trec_ge200_fp.append(fields[0])
            else:
                trec_ge200_nonfp.append(fields[0])
    pdb.set_trace()
    return

if __name__ == '__main__':
    gen_res()
    check_res()
    #check_trec_fp()


from util import *
import pdb, subprocess

#run stringtie or cufflinks (maxSens or Default modes)

#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSimChr15/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySimChr15/'

#large datasets:
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/'
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/'

alginment = '%s/hits.sorted.bam'%dataDir

cases = []
#cases.append(['stringtie', '', 'stringtie_DefaultParam'])
#cases.append(['stringtie', '--maxSens', 'stringtie_f_0_c_0.001'])
cases.append(['cufflinks', '', 'cufflinks_DefaultParam'])
#cases.append(['cufflinks', '--maxSens', 'cufflinks_F_0.001'])
#cases.append(['stringtie', '', 'stringtie_DefaultParam'])

#genomeFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/chr15.fa'
genomeFile = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/genome/hg19.fa'

#Tref = '%s/reference.fasta'%dataDir
Tref = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/reference.fasta'

nJobs = 20

def gen_res():
    for case in cases:
        resDir = '%s/%s/'%(dataDir, case[2])
        run_cmd('mkdir -p %s'%resDir)

        #pdb.set_trace()
        cmd = 'python run_extAssembler.py --assembler %s '%case[0] + \
              '%s '%case[1] + \
              '-i %s '%alginment + \
              '-g %s '%genomeFile + \
              '-O %s '%resDir + \
              '-N %d '%nJobs + \
              '-n %s '%case[0] + \
              '-NoSeperatedLines'
        print(cmd)
        run_cmd(cmd)

        #eval
        Trec = '%s/%s.fasta'%(resDir, case[0])

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
    for case in cases:
        resFile = '%s/%s/%s_res.txt'%(dataDir, case[2], case[0])
        try:
            res = subprocess.check_output('tail %s'%resFile, shell=True)
            res = res.strip().split()[-4:]
            results[(case[0], case[2])] = [int(res[0]), float(res[3])] #sens, fp
        except:
            results[(case[0], case[2])] = []

    for k, v in results.items():
        print('%s %30s \tsens %d fp %f'%(k[0], k[1], v[0], v[1]))
    return

if __name__ == '__main__':
    gen_res()
    #check_res()

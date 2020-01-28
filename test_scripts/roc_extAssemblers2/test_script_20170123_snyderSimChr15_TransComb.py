from util import *
import pdb, subprocess

dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/TransCombTest/SnyderSimChr15/tophat2/'

alginment = '%s/accepted_hits.bam'%dataDir

cases = []
#cases.append(['stringtie', '', 'stringtie_DefaultParam'])
#cases.append(['stringtie', '--maxSens', 'stringtie_f_0_c_0.001'])
#cases.append(['cufflinks', '', 'cufflinks_DefaultParam'])
#cases.append(['cufflinks', '--maxSens', 'cufflinks_F_0.001'])
#cases.append(['TransComb', '', 'TransComb_DefaultParam'])
cases.append(['TransComb', '', 'TransComb_f_1'])

genomeFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/chr15.fa'
#genomeFile = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/genome/hg19.fa'

Tref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/reference.fasta'

nJobs = 20


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

    Tref1 = '%s/reference.fasta'%(resDir)
    cmd = 'cp %s %s'%(Tref, Tref1)
    run_cmd(cmd)

    cmd = 'python filter_FP_batch.py --eval1Job ' + \
          '-t %s '%Tref1 + \
          '-r %s '%Trec + \
          '-O %s'%resDir
    run_cmd(cmd)

    cmd = 'rm %s'%(Tref1)
    run_cmd(cmd)

    print('%s done'%resDir)
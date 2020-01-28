import os

#Output
rsemExpIndex = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/simWW/exp'
#outSam = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/wwSim/simReadsOnRefTr.sam'
outSam = rsemExpIndex + '.transcript.bam'
#outSamSorted = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/wwSim/simReadsOnRefTrSorted.sam'
outSamSorted = rsemExpIndex + '.transcript.sorted.sam'
covFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/wwSim/coverage.txt'
OracleTref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/wwSim/wwSim_oracle_et25_ct100_cut0.fa'
ResDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/sens_fp_oracle/'

#Input
Tref = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/reference.fasta'
rsemTrefIndex = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref'
simReads = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/ww_snyderRef/sim_reads.fa'
TrecFiles = []
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon/reconstructed.fasta')
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/stringtie/stringtie.fasta')
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/cufflinks/cufflinks.fasta')

'''
memory error

cmd = 'python aligner_star.py -o %s '%outSam+\
      '-g %s '%Tref+\
      '-r %s '%simReads+\
      '-N 20'
os.system(cmd)
'''

cmd = 'rsem-calculate-expression --no-qualities --star -p 20 '+\
      '%s '%simReads+\
      '%s '%rsemTrefIndex+\
      '%s '%rsemExpIndex
os.system(cmd)

cmd = 'samtools sort -o %s -@ 20 %s'%(outSamSorted, outSam)
os.system(cmd)

cmd = 'samtools depth %s > %s'%(outSamSorted, covFile)
os.system(cmd)

cmd = 'python analyze_realData_sensitivity.py '+\
      '--findCoveredTr --et 25 --ct 100 --cut 0 '+\
      '--trCov %s '%covFile+\
      '--trFa %s '%Tref+\
      '-o %s'%OracleTref
os.system(cmd)


for Trec in TrecFiles:
    cmd = 'python filter_FP_batch.py --eval1Job '+\
          '-t %s '%OracleTref+\
          '-r %s '%Trec+\
          '-O %s '%ResDir
    os.system(cmd)





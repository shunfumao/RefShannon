import os

#Output
rsemExpIndex = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/simKidney/exp'
#outSam = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/kidneySim/simReadsOnRefTr.sam'
outSam = rsemExpIndex + '.transcript.bam'
#outSamSorted = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/kidneySim/simReadsOnRefTrSorted.sam'
outSamSorted = rsemExpIndex + '.transcript.sorted.sam'
covFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/kidneySim/coverage.txt'
OracleTref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/kidneySim/kidneySim_oracle_et25_ct100_cut0.fa'
ResDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/sens_fp_oracle/'

#Input
Tref = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/reference.fasta'
rsemTrefIndex = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref'
simReads1 = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_1.fq'
simReads2 = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_2.fq'
TrecFiles = []
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon/reconstructed_all.fasta')
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/stringtie/stringtie.fasta')
TrecFiles.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/cufflinks/cufflinks.fasta')

'''
memory error

cmd = 'python aligner_star.py -o %s '%outSam+\
      '-g %s '%Tref+\
      '-r1 %s '%simReads1+\
      '-r2 %s '%simReads2+\
      '-N 20'
os.system(cmd)
'''

cmd = 'rsem-calculate-expression --star --paired-end -p 20 '+\
      '%s '%simReads1+\
      '%s '%simReads2+\
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





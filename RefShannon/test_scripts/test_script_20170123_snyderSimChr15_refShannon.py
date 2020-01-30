from util import *

#Dir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2_unsorted/'
#bam_file = '%s/accepted_hits.bam'%Dir
#sam_file = '%s/accepted_hits.sam'%Dir 

Dir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/star/'
sam_file = '%s/hits.sam'%Dir 

#cmd = 'samtools view -h -o %s -@ 20 %s'%(sam_file, bam_file)
#run_cmd(cmd)

genome_file = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/chr15.fa'
target = 'chr15'

Tref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/reference.fasta'

outDir = '%s/RefShannon_force_SE_ROC/'%Dir
run_cmd('mkdir -p %s'%outDir)

pairedStr = '-paired'

nJobs = 20

#Fvals = ['0', '0.25', '0.4', '0.96']
Fvals = ['0', '0.25', '0.4', '0.96']

'''
cmd = 'python refShannon.py --batch '+\
      '-i %s '%sam_file+\
      '-g %s '%genome_file+\
      '-O %s '%outDir+\
      '%s '%pairedStr+\
      '-target %s '%target+\
      '-N 20'
run_cmd(cmd)
'''

for fval in Fvals:
      resDir = outDir+'/%s/'%fval
      run_cmd('mkdir -p %s'%resDir)

      genomeFile = genome_file

      #sam --> sg
      sg_dir = '%s/intermediate/'%(resDir)
      cmd = 'python refShannon.py --sg '+ \
            '-i %s '%sam_file+ \
            '-g %s '%genomeFile+ \
            '-O %s '%resDir+ \
            '%s -target chr15 '%(pairedStr)+ \
            '-F %s'%fval
      run_cmd(cmd)

      #pdb.set_trace()

      #sg --> transcripts
      cmd = 'python refShannon.py --sf -I %s '%sg_dir+ \
            '-g %s '%genomeFile+ \
            '-O %s '%resDir+ \
            '-N %d '%nJobs+ \
            '-target chr15 '+ \
            '-F %s'%fval
      run_cmd(cmd)

      #eval
      Trec = '%s/reconstructed.fasta'%(resDir)
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

from RefShannon.util import *
import pdb, subprocess, os

PATH = os.path.dirname(__file__)
ROOT = parent_dir(PATH)

def exAssembler_roc(args):
  alignment, genomeFile, cases, resDir, reference, nJobs = args
  
  for case in cases:
    resDir2 = '%s/%s/%s/'%(resDir, case[0], case[2])
    run_cmd('mkdir -p %s'%resDir2)
    cmd = 'python %s/run_extAssembler.py '%ROOT + \
          '--assembler %s '%case[0] + \
          '%s '%case[1] + \
          '-i %s '%alignment + \
          '-g %s '%genomeFile + \
          '-O %s '%resDir2 + \
          '-N %d '%nJobs + \
          '-n %s '%case[0] + \
          '-NoSeperatedLines'
    print(cmd)
    run_cmd(cmd)

    #eval
    Trec = '%s/%s.fasta'%(resDir2, case[0])

    Tref1 = '%s/reference.fasta'%(resDir2)
    cmd = 'cp %s %s'%(reference, Tref1)
    run_cmd(cmd)

    cmd = 'python %s/filter_FP_batch.py --eval1Job '%ROOT + \
          '-t %s '%Tref1 + \
          '-r %s '%Trec + \
          '-O %s'%resDir2
    run_cmd(cmd)

    cmd = 'rm %s'%(Tref1)
    run_cmd(cmd)

    print('%s done'%resDir)
  return
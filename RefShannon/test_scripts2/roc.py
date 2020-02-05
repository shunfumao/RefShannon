from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

from RefShannon.util import *
import pdb, subprocess, os
from RefShannon.run_parallel_cmds import run_cmds

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
    # pdb.set_trace()
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

def gen_res(args):
  rootResDir, T, sam_file, genomeFile, pairedStr, chrom, Tref, nJobs = args

  sg_cmds = []
  sf_cmds = []
  cp_ref_cmds = []
  ev_cmds = []
  rm_ref_cmds = []
  for t in T:
    resDir = '%s/%s/'%(rootResDir, t)
    run_cmd('mkdir -p %s'%resDir)

    #sam --> sg
    sg_dir = '%s/intermediate/'%(resDir)
    cmd = 'python %s/refShannon.py --sg '%ROOT + \
      '-i %s '%sam_file+ \
      '-g %s '%genomeFile+ \
      '-O %s '%resDir+ \
      '%s -target %s'%(pairedStr, chrom)
    sg_cmds.append(cmd)

    #sg --> transcripts
    cmd = 'python %s/refShannon.py --sf -I %s '%(ROOT, sg_dir)+ \
          '-g %s '%genomeFile+ \
          '-O %s '%resDir+ \
          '-N %d '%nJobs+ \
          '-target %s '%chrom+ \
          '--outputFasta '
    sf_cmds.append(cmd)

    #eval
    Trec = '%s/reconstructed.fasta'%(resDir)

    Tref1 = '%s/reference.fasta'%(resDir)
    cmd = 'cp %s %s'%(Tref, Tref1)
    cp_ref_cmds.append(cmd)

    cmd = 'python %s/filter_FP_batch.py --eval1Job '%ROOT + \
          '-t %s '%Tref1 + \
          '-r %s '%Trec + \
          '-O %s'%resDir
    ev_cmds.append(cmd)

    cmd = 'rm %s'%(Tref1)
    rm_ref_cmds.append(cmd)

  pdb.set_trace()
  run_cmds(sg_cmds,noJobs=nJobs)
  # pdb.set_trace()
  run_cmds(sf_cmds,noJobs=nJobs)
  # pdb.set_trace()
  run_cmds(cp_ref_cmds,noJobs=nJobs)
  # pdb.set_trace()
  run_cmds(ev_cmds,noJobs=nJobs)
  # pdb.set_trace()
  run_cmds(rm_ref_cmds,noJobs=nJobs)
  # pdb.set_trace()
  
  return

#extract res from gen_res()
def check_res(args):
  rootResDir, rootResFile, T = args

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
  #print(results)
  results = results.items()
  results.sort(key=lambda x:x[0])
  for t, sens_fp in results:
    print('t=%f\t%d\t%f\t'%(float(t), int(sens_fp[0]), float(sens_fp[1])))
  return

def refShannon_roc(args):
  rootResDir, T, sam_file, genomeFile, pairedStr, chrom, Tref, nJobs = args
  run_cmd('mkdir -p %s'%rootResDir)

  args1 = (rootResDir, T, sam_file, genomeFile, pairedStr, chrom, Tref, nJobs)
  gen_res(args1)

  rootResFile = '%s/res.txt'%rootResDir
  args2 = (rootResDir, rootResFile, T)
  check_res(args2)

  return
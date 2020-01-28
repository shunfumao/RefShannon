from util import *
import pdb, subprocess
from run_parallel_cmds import *

#tune FLOW_RATIO_THRESHOLD to see sens/fp change of refShannon
#onto large dataset

#gen_sg2 merge type II OFF
#gen_sg2 merge type II tune cov approx 4 + SF lessP

#sim
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/'

#real - for sens analysis
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/'

chrs_dir = '%s/chrs_hits/'%dataDir #chrs_dir/:target/hits.sam split from dataDir/hits.sam sim ww
#chrs_dir = '%s/chrs/'%dataDir #chrs_dir/:target/hits.sam split from dataDir/hits.sam

#readl ww
#chrs_dir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_tmp/' #chrs_dir/:target/hits.sam split from dataDir/hits.sam

genomeDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/'

Tref = '%s/reference.fasta'%dataDir

#condition = 'mergeII_approx4_SF_lessP'
#condition = 'OrigMergeII_SF_lessP'
condition = 'mergeII_approx4_SF_lessP_code_merge_test'

rootResDir = '%s/refShannon_test_%s/'%(dataDir, condition)
rootResFile = '%s/res.txt'%rootResDir

pairedStr = '' #'-paired' #'-paired' or ''

test_chrs = '-chrs chr1_gl000191_random,chr1_gl000192_random' #test or ''

nJobs = 20

#T = ['0', '0.1']# -- test
#T = ['0', '0.1', '0.25', '0.4', '0.96']
#T = ['0', '0.25', '0.96']
T = ['0', '0.25', '0.4', '0.96']
#T = ['0']
T = [[condition,t] for t in T]

def gen_res():

    #pdb.set_trace()
    
    '''resDir = '%s/%s/chrs/'%(rootResDir, T[0][0]) #refShannon_test/(cond)/chrs/:target/intermediate/
    run_cmd('mkdir -p %s'%resDir)

    cmd = 'python refShannon.py --sg '+ \
          '-I %s '%chrs_dir+ \
          '%s '%test_chrs+ \
          '-G %s '%genomeDir+ \
          '-O %s '%resDir+ \
          '%s -N %d'%(pairedStr, nJobs)
    run_cmd(cmd)'''

    sg_cmds = []
    sf_cmds = []
    pool_cmds = []
    cp_ref_cmds = []
    ev_cmds = []
    rm_ref_cmds = []
        
    for cond, t in T:

        resDir = '%s/%s/%s/chrs/'%(rootResDir, cond, t) #refShannon_test/(cond)/(t)/chrs/:target/intermediate/
        run_cmd('mkdir -p %s'%resDir)

        cmd = 'python refShannon.py --sg '+ \
        '-I %s '%chrs_dir+ \
        '%s '%test_chrs+ \
        '-G %s '%genomeDir+ \
        '-O %s '%resDir+ \
        '%s -N %d -F %s'%(pairedStr, nJobs, t)
        #run_cmd(cmd)
        sg_cmds.append(cmd)

        algoResDir = '%s/%s/%s/chrs/'%(rootResDir, cond, t) #refShannon_test/(cond)/(t)/chrs/:target/algo_output
        pooledAlgoResDir = '%s/%s/%s/'%(rootResDir, cond, t) #refShannon_test/(cond)/(t)/: reconstructed_all.fasta
        run_cmd('mkdir -p %s'%algoResDir)
        run_cmd('mkdir -p %s'%pooledAlgoResDir)

        # gen Trec for each target
        cmd = 'python refShannon.py --sf '+ \
              '-Is %s '%resDir+ \
              '%s '%test_chrs+ \
              '-O %s '%algoResDir+ \
              '-N %d -F %s'%(nJobs, t)
        #run_cmd(cmd)
        sf_cmds.append(cmd)

        #pool res of fasta
        cmd = 'python refShannon.py --eval '+ \
              '-I %s '%algoResDir+ \
              '%s '%test_chrs+ \
              '-r %s '%Tref+ \
              '-O %s '%pooledAlgoResDir+ \
              '--poolFastaOnly'
        #run_cmd(cmd)
        pool_cmds.append(cmd)

        #eval
        Trec = '%s/reconstructed_all.fasta'%(pooledAlgoResDir)
        Tref1 = '%s/reference.fasta'%(pooledAlgoResDir)
        
        cmd = 'cp %s %s'%(Tref, Tref1)
        #run_cmd(cmd)
        cp_ref_cmds.append(cmd)

        cmd = 'python filter_FP_batch.py --eval1Job ' + \
              '-t %s '%Tref1 + \
              '-r %s '%Trec + \
              '-O %s'%pooledAlgoResDir
        #run_cmd(cmd)
        ev_cmds.append(cmd)

        cmd = 'rm %s*'%(Tref1)
        #run_cmd(cmd)
        rm_ref_cmds.append(cmd)

    run_cmds(sg_cmds,noJobs=nJobs)
    run_cmds(sf_cmds,noJobs=nJobs)
    run_cmds(pool_cmds,noJobs=nJobs)
    run_cmds(cp_ref_cmds,noJobs=nJobs)
    run_cmds(ev_cmds,noJobs=nJobs)
    run_cmds(rm_ref_cmds,noJobs=nJobs)

    return

#extract res from gen_res()
def check_res():
    results = {}
    with open(rootResFile, 'w') as f:
        for cond, t in T:
            resFile = '%s/%s/%s/reconstructed_all_res.txt'%(rootResDir, cond, t)
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

if __name__ == '__main__':
    gen_res()
    check_res()
    #check_trec_fp()


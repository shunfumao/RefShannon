from util import *
import pdb, subprocess
from run_parallel_cmds import *
import os

#tune FLOW_RATIO_THRESHOLD to see sens/fp change of refShannon
#onto large dataset
#gen_sg2 merge type II OFF

#sim
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/'
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/'

#real - for sens analysis
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/'
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/'

#chrs_dir = '%s/chrs_hits/'%dataDir #chrs_dir/:target/hits.sam split from dataDir/hits.sam
#chrs_dir = '%s/chrs/'%dataDir #chrs_dir/:target/hits.sam split from dataDir/hits.sam

#readl ww
#chrs_dir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_tmp/' #chrs_dir/:target/hits.sam split from dataDir/hits.sam
#real snyder
chrs_dir = '/data1/shunfu1/ref_shannon_modi/data//snyder_0729a/chrs_0807b/'

genomeDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/'

Tref = '%s/reference.fasta'%dataDir

condition = 'mergeII_partial_tune_cov_approx4_trial_cov_unlimit_dump'

rootResDir = '%s/refShannon_test_%s/'%(dataDir, condition)
rootResFile = '%s/res.txt'%rootResDir

pairedStr = '-paired' #'-paired' #'-paired' or ''

test_chrs = '' #'-chrs chr1_gl000191_random,chr1_gl000192_random' #test or ''

nJobs = 1

#T = ['0', '0.1']# -- test
#T = ['0', '0.1', '0.25', '0.4', '0.96']
#T = ['0', '0.25', '0.96']
#T = ['0', '0.25', '0.96']
T = ['0']
T = [[condition,t] for t in T]

def merge_cov_dmp():

    for cond, t in T:

        resCovFile = 'dmp/resCovFile.txt'
        run_cmd('touch %s'%resCovFile)

        resDir = '%s/%s/%s/chrs/'%(rootResDir, cond, t) #refShannon_test/nomergetypeII/(t)/chrs/:target/intermediate/

        target_list = os.listdir(chrs_dir)

        for target in target_list:
            print('target: %s'%target)
            cov_file = '%s/%s/cov_dmp/cov_dmp.txt'%(resDir, target)
            run_cmd('cat %s >> %s'%(cov_file, resCovFile))

merge_cov_dmp()
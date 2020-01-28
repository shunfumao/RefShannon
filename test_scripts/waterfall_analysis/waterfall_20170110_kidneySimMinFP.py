from util import run_cmd
import pdb

choice = 1 # 0 - sens 1 - fp

files = []
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.96/reconstructed_all_fp_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/stringtie_DefaultParam/stringtie_fp_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/cufflinks_DefaultParam/cufflinks_fp_log.txt')

files = ' '.join(files)

cmd = 'python waterfall.py -%d %s'%(choice, files)

pdb.set_trace()

run_cmd(cmd)

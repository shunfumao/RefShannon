from util import run_cmd
import pdb

choice = 0 # 0 - sens 1 - fp

files = []
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/stringtie_f_0_c_0.001/stringtie_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/cufflinks_F_0.001/cufflinks_log.txt')

files = ' '.join(files)

cmd = 'python waterfall.py -%d %s'%(choice, files)

pdb.set_trace()

run_cmd(cmd)

from util import run_cmd
import pdb

dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/'

files = []
files.append('%s/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all_log.txt'%dataDir)
files.append('%s/Per_Log_fpLog/cufflinks_log.txt'%dataDir)
files.append('/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_c_0.001/stringtie_log.txt')

files = ' '.join(files)

cmd = 'python waterfall.py -0 %s'%files

pdb.set_trace()

run_cmd(cmd)

from util import run_cmd
import pdb

files = []
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/_res_0726b_stringtie_f_0_c_0.001/stringtie_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/Per_Log_fpLog/cufflinks_log.txt')


files = ' '.join(files)

cmd = 'python waterfall.py -0 %s'%files

pdb.set_trace()

run_cmd(cmd)

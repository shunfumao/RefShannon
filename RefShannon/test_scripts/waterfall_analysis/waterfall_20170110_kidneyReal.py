from util import run_cmd
import pdb

files = []
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/stringtie/stringtie_log.txt')
files.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/cufflinks/cufflinks_log.txt')

files = ' '.join(files)

cmd = 'python waterfall.py -0 %s'%files

pdb.set_trace()

run_cmd(cmd)

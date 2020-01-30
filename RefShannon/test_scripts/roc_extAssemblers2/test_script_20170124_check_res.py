from util import *
import pdb, subprocess
from run_parallel_cmds import *

#check res of refShannon/stringtie/TransCOmb/CLASS2

#sim
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/TransCombTest/SnyderSimChr15/'
rootResFile = '%s/res.txt'%(dataDir)

cases = []
cases.append(['star',  'RefShannon', 'reconstructed_res.txt'])
cases.append(['star',  'stringtie_f_0_c_0.001', 'stringtie_res.txt'])
cases.append(['star',  'CLASS2_DefaultParam', 'CLASS2_res.txt'])
cases.append(['star',  'CLASS2_F_0', 'CLASS2_res.txt'])

#cases.append(['tophat2',  'RefShannon', 'reconstructed_res.txt'])
cases.append(['tophat2',  'stringtie_f_0_c_0.001', 'stringtie_res.txt'])
cases.append(['tophat2',  'CLASS2_DefaultParam', 'CLASS2_res.txt'])
cases.append(['tophat2',  'CLASS2_F_0', 'CLASS2_res.txt'])
cases.append(['tophat2',  'TransComb_DefaultParam', 'TransComb_res.txt'])
cases.append(['tophat2',  'TransComb_f_1', 'TransComb_res.txt'])

#extract res from gen_res()
def check_res():
    results = []
    with open(rootResFile, 'w') as f:
        for case in cases:
            resFile = '%s/%s/%s/%s'%(dataDir, case[0], case[1], case[2])
            try:
                #pdb.set_trace()
                res = subprocess.check_output('tail %s'%resFile, shell=True)
                res = res.strip().split()[-4:]
                results.append([resFile, int(res[0]), float(res[3])])#sens, fp
                f.write('%s\n\tsens\t%d\tfp\t%f\n'%(resFile, int(res[0]), float(res[3])))
            except:
                print('except at %s'%resFile)
    #print(results)
    for f, s, fp in results:
        print('%s\n\t%d\t%f\t'%(f, s, fp))
    return

if __name__ == '__main__':
    check_res()
    #check_trec_fp()


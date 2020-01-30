#re-check the ROC curve using strict sensitivity criteria
from util import *

#DataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/'
#DataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/'
DataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/'

Cases = []
Cases.append(['stringtie_f_0_c_0.001', 'stringtie', 'StringTie - Max Sensitivity'])
Cases.append(['stringtie_DefaultParam', 'stringtie', 'StringTie - Default'])
Cases.append(['cufflinks_F_0.001', 'cufflinks', 'Cufflinks - Max Sensitivity'])
Cases.append(['cufflinks_DefaultParam', 'cufflinks', 'Cufflinks - Default'])
Cases.append(['refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0','reconstructed_all','RefShannon - 0'])
Cases.append(['refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.25','reconstructed_all','RefShannon - 0.25'])
Cases.append(['refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.4','reconstructed_all','RefShannon - 0.4'])
Cases.append(['refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.96','reconstructed_all','RefShannon - 0.96'])

for case in Cases:
    print('-----%s-----'%(case[2]))
    logFile = '%s/%s/%s_log.txt'%(DataDir, case[0], case[1])
    cmd = 'python filter_FP_batch.py --calcSensWrapper --logFile %s'%(logFile)
    run_cmd(cmd)    





import os
'''
scripts to analyze real data sensitivity based on grouping oracle ref transcripts by abundance or isoform multiplicity
'''

#oracleFa = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/oracle/snyder_oracle_et25_ct100_cut0.fa' ## snyder
oracleFa = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/ww/ww_oracle_et25_ct100_cut0.fa' ## ww
#oracleFa = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/kidney/oracle/kidney_oracle_et25_ct100_cut0.fa'

#numIsoFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/geneIsoNum.txt'
numIsoFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/ww/geneIsoNum.txt'

# snyder, hg19
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/Per_Log_fpLog/reconstructed_log.txt'
# snyder, hg19, nomerge type II
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/refShannon_test/nomergetypeII/0/reconstructed_all_log.txt'
# snyder, hg19, merge II tune cov approx 3
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/refShannon_test_mergeII_partial_tune_cov_approx3/mergeII_partial_tune_cov_approx3/0/reconstructed_all_log.txt'

# kidney, hg19, no merge type II
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/refShannon_test/nomergetypeII/0/reconstructed_all_log.txt'

# ww, hg19
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/Per_Log_fpLog/reconstructed_log.txt'
# ww, hg19, nomerge typeII
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/refShannon_test/nomergetypeII/0/reconstructed_all_log.txt'
# ww, hg19, merge II tune cov approx 3
#recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/refShannon_test_mergeII_partial_tune_cov_approx3/mergeII_partial_tune_cov_approx3/0/reconstructed_all_log.txt'
# ww, hg19, merge II approx 4, less P
recLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all_log.txt'

#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_c_0.001/stringtie_log.txt'
#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/Per_Log_fpLog/cufflinks_log.txt'

#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/_res_0726b_stringtie_f_0_c_0.001/stringtie_log.txt'
#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/Per_Log_fpLog/cufflinks_log.txt'

cmpLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/Per_Log_fpLog/reconstructed_log.txt' # ww hg9 (before mergeII_approx4_SF_lessP)
#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/Per_Log_fpLog/reconstructed_log.txt'
#cmpLog = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/refShannon/reconstructed_log.txt' #kidney hg19

#expFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/exp.isoforms.results'
expFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/ww/exp_kal_out/abundance.tsv'
#expFile = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/kidney_snyderRef/exp.isoforms.results'

#expFormat = 'RSEM'
expFormat = 'KAL'

#L=100 #kidney
#L = 101 #snyder
L = 50
S = 0
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/sens_refS_stringtie_L%d_S%d'%(L,S)
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/sens_refS_cufflinks_L%d_S%d'%(L,S)
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/sens_refShannonNoMergeII_refShannon_L%d_S%d'%(L,S)
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/sens_refShannonMergeIItuneCovApprox3_refShannon_L%d_S%d'%(L,S)

#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/sens_refShannonNoMergeII_refShannon_L%d_S%d'%(L,S)
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/sens_refShannonMergeIItuneCovApprox3_refShannon_L%d_S%d'%(L,S)
#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/sens_refShannonMergeIItuneCovApprox4LessP_refShannon_L%d_S%d'%(L,S)

#outFileStem = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/sens_refShannonNoMergeII_refShannon_L%d_S%d'%(L,S)

ablist = '0,1e10'
#ablist = '0.0,25.1,46.4,67.6,99.4,136.4,191.5,283.7,470.4,947.3,1e10'
#ablist = '0,50,100,500,1e10'#25,25,30+,20- percentile
ablist = '0,46.38,99.41,191.51,470.44,1e10'#20 percentile per group

#snyder
#ablist = '0,0.7,2.2,4.2,6.9,10.9,17.8,30.0,56.8,142.5,1e10'
#ablist = '0,0.49,2.9,6.35,11.32,18.44,29.73,49.08,87.39,191.34,1e10' #kidney percentile

#IM_list=[[1,5],[6,7],[8,9],[10,11],[12,13],[14,15],[16,18],[19,22],[23,28],[29,77]]
IM_list = [[1,100]]
#IM_list = [[1,1],[2,3],[4,8]]
#IM_list = [[1,2],[3,4],[5,7],[8,9],[10,11],[12,13],[14,16],[17,19],[20,25],[26,77]]

for iL, iH in IM_list:

    print('---------- iL=%d, iH=%d ----------'%(iL, iH))

    cmd = 'python analyze_realData_sensitivity.py '+\
          '--performancePlotExpress '+\
          '--oracleFa %s '%oracleFa+\
          '--numIsoFile %s '%numIsoFile+\
          '--recLog %s '%recLog+\
          '--cmpLog %s '%cmpLog+\
          '--expFile %s '%expFile+\
          '--expFormat %s '%expFormat+\
          '-L %d -S %d --isoLow %d --isoHigh %d '%(L,S,iL,iH)+\
          '--o %s_isoL%d_isoH%d.txt '%(outFileStem, iL, iH)+\
          '--ablist %s'%(ablist)

    os.system(cmd)

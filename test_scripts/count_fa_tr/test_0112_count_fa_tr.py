
def cnt_fasta(filename):
    minTr=200 #200
    res=0
    seq = []
    next_name = '_'
    with open(filename) as f:
        for line in f:
            if line[0]=='>':
                seq=''.join(seq)
                if next_name != '_' and len(seq)>=minTr:
                    res+=1
                seq=[]
                next_name=line.split()[0][1:]
            else:
                seq.append(line.strip())
    seq=''.join(seq)
    if next_name != '_' and len(seq)>=minTr:
        res+=1
    print(filename)
    print(res)
    return res

filenames = []

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/cufflinks/cufflinks.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/cufflinks/cufflinks.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/cufflinks/cufflinks.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/cufflinks_F_0.001/cufflinks.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/cufflinks_DefaultParam/cufflinks.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/cufflinks_F_0.001/cufflinks.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/cufflinks_DefaultParam/cufflinks.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/cufflinks_F_0.001/cufflinks.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/cufflinks_DefaultParam/cufflinks.fasta')


'''
filenames.append('/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/_res_0726b_stringtie_f_0_c_0.001/stringtie.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_c_0.001/stringtie.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/stringtie/stringtie.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/stringtie_f_0_c_0.001/stringtie.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/stringtie_DefaultParam/stringtie.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/stringtie_f_0_c_0.001/stringtie.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/stringtie_DefaultParam/stringtie.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/stringtie_f_0_c_0.001/stringtie.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/stringtie_DefaultParam/stringtie.fasta')
'''

'''
filenames = []
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/ww_1002a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidney_0916a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.25/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.4/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSim_1122a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.96/reconstructed_all.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.25/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.4/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSim_1018a/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.96/reconstructed_all.fasta')

filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.25/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.4/reconstructed_all.fasta')
filenames.append('/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/kidneySim/refShannon_test_mergeII_approx4_SF_lessP/mergeII_approx4_SF_lessP/0.96/reconstructed_all.fasta')
'''

for filename in filenames:
    cnt_fasta(filename)

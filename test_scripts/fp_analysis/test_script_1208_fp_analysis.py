
from util import run_cmd
from filter_FP_batch import calcFP
import sys, pdb

## functions
def load_kal_ab(kal_ab, transcripts):

    with open(kal_ab, 'r') as f:
        f.readline()
        for line in f:
            tokens = line.strip().split()
            transcript_id = tokens[0]
            length = int(tokens[1])#0
            eff_len = float(tokens[2])#1
            #est_cnts = float(tokens[3])#2
            cov = float(tokens[3]) / max(1,(length-RL)) * RL#2
            tpm = float(tokens[4])#3
            transcripts[transcript_id]=[length, eff_len, cov, tpm]
    return

def label_fp(fp_log, transcripts):
    tr_labels = {}
    calcFP(fp_log, mtl=200, tr_labels=tr_labels)
    for k,v in tr_labels.items():
        transcripts[k] += [v]#4 -1: skip 0: non-fp 1: fp
    return

def cal_cov_ratio(dep_file, transcripts):

    #pdb.set_trace()
    tr_hits = {}
    for line in open(dep_file, 'r'):
        fields = line.strip().split()
        tr_name = fields[0]
        #gpos = int(fields[1])
        cov = int(fields[2])
        if cov>0: tr_hits[tr_name] = tr_hits.get(tr_name, 0)+1

    #pdb.set_trace()
    for tr_name, hits_len in tr_hits.items():
        tr_len = transcripts[tr_name][0]
        transcripts[tr_name] += [float(hits_len)/tr_len]#5 cov ratio

    for tr_name, v in transcripts.items():
        if len(v)<6:
            #pdb.set_trace()
            v += [0.0]

    return

#4
#check hist-feature curve of fp and non-fp
#to see which feature (e.g. cov ratio) to filter that may bring a decreasing fp_ratio-cutoff  curve
def check_fp():

    import numpy as np 
    import pdb
    import matplotlib.pyplot as plt

    itms = []

    dmpRes = 'dmp/FP_analysis_wwSimChr15_0.4/dmp_stat_FP_vs_nonFP.txt'

    with open(dmpRes, 'r') as f:

        f.readline()

        for line in f:
            tokens = line.strip().split()
            itms.append(tokens)


    choice = 6

    if choice == 3:
        xlab = 'abundance'
        bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 1, 2, 5, 10, 20, 50, 100,200,500,1000,2000,5000,1e4, 2e4, 5e4, 1e10]
    elif choice == 1:
        xlab = 'len'
        bins = [0, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6]
    elif choice == 2:
        xlab = 'eff_len'
        bins = [2e1, 3e1, 4e1, 5e1, 6e1, 7e1, 8e1, 9e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]
    elif choice == 4:
        xlab = 'tpm'
        bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 2e1, 3e1, 4e1, 5e1, 6e1, 7e1, 8e1, 9e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]
    elif choice == 6:
        xlab = 'cov ratio'
        bins = [0, 1e-3, 1e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1]

    fp0 = [float(a[choice]) for a in itms if int(a[5])==0] #non fp
    fp1 = [float(a[choice]) for a in itms if int(a[5])==1] #fp

    print('fp ratio: %f'%(float(len(fp1))/(len(fp0)+len(fp1))))

    #pdb.set_trace()

    #'''
    #modified fp ratio versus cutoff threshold:
    weights = fp0 + fp1; weights.sort(key=lambda x: x)
    cuts = [0, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 9.5e-1, 9.9e-1]
    fp_ratios = []
    fp0_ft_perc = [] #percentage of fp0 that are filtered
    fp1_ft_perc = []
    for cut in cuts:
        threshold = weights[int(len(weights)*cut)]
        fp0_ft = [float(a[choice]) for a in itms if int(a[5])==0 and float(a[choice])>=threshold] #non fp
        fp1_ft = [float(a[choice]) for a in itms if int(a[5])==1 and float(a[choice])>=threshold] #fp
        try:
            fpr = float(len(fp1_ft))/(len(fp0_ft)+len(fp1_ft))
            fp_ratios.append(fpr)

            fp0_ft_perc.append(1-float(len(fp0_ft))/len(fp0))
            fp1_ft_perc.append(1-float(len(fp1_ft))/len(fp1))

            #pdb.set_trace()
            print('cut %f, threshold %f, fp_ft ratio: %f'%(cut, threshold, fpr))
        except:
            continue

    fig, ax = plt.subplots()
    plt.plot(cuts[0:len(fp_ratios)], fp_ratios, marker='o', label='fp ratio')
    plt.plot(cuts[0:len(fp0_ft_perc)], fp0_ft_perc, marker='o', label='fp0_ft_perc')
    plt.plot(cuts[0:len(fp1_ft_perc)], fp1_ft_perc, marker='o', label='fp1_ft_perc')
    ax.set_xlabel('cut ratio (feature: %s)'%xlab)
    ax.set_ylabel('fp ratio')
    ax.set_yscale('log')
    plt.legend()
    plt.show()
    return
    #'''


    #pdb.set_trace()

    [hist0, bins0] = np.histogram(fp0, bins=bins)
    [hist1, bins1] = np.histogram(fp1, bins=bins)
    fig, ax = plt.subplots()
    plt.plot(bins0[1:], hist0, marker='o', label='non fp')
    plt.plot(bins1[1:], hist1, marker='o', color='red', label='fp')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(xlab)
    ax.set_ylabel('frequency')
    plt.legend()
    plt.show()

    #pdb.set_trace()

    return

'''
usage

python test..py mode
'''

## global config
#Tref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/reference.fasta'
Tref = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSimChr15/reference.fasta'

#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon_test/0.96/' #minFP
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon/algo_output/' #maxSens
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon_test/0.3/'

dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/wwSimChr15/refShannon_test/0.4/'

prefix = 'reconstructed'
nJobs = 20

Trec = '%s/%s.fasta'%(dataDir, prefix)

#reads = ['/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads_1.fq',
#         '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads_2.fq']


reads = ['/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/ww_snyderRef/sim_reads.fa']
isFa = True #fa or fq

#RL = 101
RL = 50

#dmpRes = 'dmp/dmp_stat_FP_vs_nonFP_maxSens.txt'
#dmpRes = 'dmp/dmp_stat_FP_vs_nonFP_minFP.txt'
#dmpRes = 'dmp/wwSimChr15_dmp_stat_FP_vs_nonFP.txt'

if __name__ == '__main__':

    #pdb.set_trace()

    if '1' in sys.argv:

        #if '-d' in sys.argv:
        #    dataDir = sys.argv[sys.argv.index('-d')+1]

        #generate kal abundance

        KAL_prefix = '%s/%s'%(dataDir, prefix)

        if len(reads)==2:

            cmd = 'python analyze_simData_ROC.py --kallisto ' + \
                  '-o %s '%KAL_prefix + \
                  '--rec %s '%Trec + \
                  '-r1 %s '%reads[0] + \
                  '-r2 %s '%reads[1] + \
                  '-t %d'%nJobs

        else:

            cmd = 'python analyze_simData_ROC.py --kallisto ' + \
                  '-o %s '%KAL_prefix + \
                  '--rec %s '%Trec + \
                  '-r1 %s '%reads[0] + \
                  '-t %d'%nJobs

        run_cmd(cmd)

    elif '2' in sys.argv:

        #pdb.set_trace()

        #if '-d' in sys.argv:
        #    dataDir = sys.argv[sys.argv.index('-d')+1]
        #Build hisat file
        run_cmd('hisat2-build %s %s/%s.hisat'%(Trec, dataDir, prefix))

        #Align
        if len(reads)==2:
            run_cmd('hisat2 -p %d --no-spliced-alignment --no-discordant -q --fr -x %s/%s.hisat -1 %s -2 %s -S %s/%s.sam'%(nJobs, dataDir, prefix, reads[0], reads[1], dataDir, prefix))
        else:
            if isFa==True:
                run_cmd('hisat2 -p %d --no-spliced-alignment -f -x %s/%s.hisat -U %s -S %s/%s.sam'%(nJobs, dataDir, prefix, reads[0], dataDir, prefix))

        
        #Process SAM / BAM file to get depth information
        #run_cmd('samtools view -@ %d -bS -f 0x2 %s/%s.sam > %s/%s.bam'%(nJobs, dataDir, prefix, dataDir, prefix))
        #try samtools view -b -@ 20 reconstructed.sam > reconstructed.bam 
        #or samtools sort -@ 20 -o reconstructed.sorted.bam reconstructed.sam directly
        run_cmd('samtools sort -@ %d -o %s/%s.sorted.bam %s/%s.sam'%(nJobs, dataDir, prefix, dataDir, prefix))
        run_cmd('samtools depth %s/%s.sorted.bam > %s/%s.depth'%(dataDir, prefix, dataDir, prefix))

    elif '3' in sys.argv:

        #pdb.set_trace()
        #if '-d' in sys.argv:
        #    dataDir = sys.argv[sys.argv.index('-d')+1]

        transcripts = {} #key - tr_id val - [param0, param1, ...]
        kal_ab = '%s/%s_kal_out/abundance.tsv'%(dataDir, prefix)
        load_kal_ab(kal_ab, transcripts)
        print('load_kal_ab done')
        #pdb.set_trace()

        fp_log = '%s/%s_fp_log.txt'%(dataDir, prefix)
        label_fp(fp_log, transcripts)
        print('label_fp done')

        dep_file = '%s/%s.depth'%(dataDir, prefix)
        cal_cov_ratio(dep_file, transcripts)
        print('cal_cov_ratio done')
        #pdb.set_trace()

        #dmp
        #dmpRes to be checked by dmp/check_fp.py
        dmpRes = '%s/dmp_stat_FP_vs_nonFP.txt'%dataDir
        with open(dmpRes, 'w') as f_dmpRes:

            st = '#tr_id\tlen\teff_len\tcov\ttpm\tisFP\tcovRatio\t'
            f_dmpRes.write(st+'\n')

            for k,v in transcripts.items():
                st = '%s\t'%k
                st += '%d\t'%v[0]
                st += '%f\t'%v[1]
                st += '%f\t'%v[2]
                st += '%f\t'%v[3]
                st += '%d\t'%v[4]
                st += '%f\t'%v[5]
                f_dmpRes.write(st+'\n')

        print('dmpRes done -- to exit')
        #pdb.set_trace()

        #ab_fp = [itm[1][2] for itm in transcripts.items() if itm[1][4]==1]
        #ab_non_fp = [itm[1][2] for itm in transcripts.items() if itm[1][4]==0]

    elif '4' in sys.argv:

        #plot hist-feature curves for fp/non-fp transcripts
        #run on local pc

        check_fp()

    elif '5' in sys.argv:

        #further filter & check perf

        from filter_FP.filter_FP import write_filtered_tr

        #if '-d' in sys.argv:
        #    dataDir = sys.argv[sys.argv.index('-d')+1]

        #pdb.set_trace()

        dep_file = '%s/%s.depth'%(dataDir, prefix)
        #threshold = '0.95'
        threshold = '0.93' #'0.90'
        Trec_filterCovRatioDir = '%s/filteredCovRatio/%s/'%(dataDir, threshold)
        run_cmd('mkdir -p %s'%Trec_filterCovRatioDir)

        Trec_filterCovRatio = '%s/%s.fasta'%(Trec_filterCovRatioDir, prefix)
        Trec_filterCovRatioLog = '%s/filteredCovRatio/%s/%s.fasta.log'%(dataDir, threshold, prefix)

        #gen Trec
        write_filtered_tr(dep_file, Trec, Trec_filterCovRatio, Trec_filterCovRatioLog, threshold=float(threshold))

        #pdb.set_trace()

        #eval
        cmd = 'python filter_FP_batch.py --eval1Job ' + \
              '-t %s '%Tref + \
              '-r %s '%Trec_filterCovRatio + \
              '-O %s'%(Trec_filterCovRatioDir)
        run_cmd(cmd)

    elif 'batch' in sys.argv:

        steps = [1,2,3,5]
            
        for step in steps:
            cmd = 'python test_script_1208_fp_analysis.py %d'%(step)
            run_cmd(cmd)

        
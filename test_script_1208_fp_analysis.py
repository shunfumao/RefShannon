
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

def check_fp():

    import numpy as np 
    import pdb
    import matplotlib.pyplot as plt

    itms = []

    with open(dmpRes, 'r') as f:

        f.readline()

        for line in f:
            tokens = line.strip().split()
            itms.append(tokens)

    #choice = 3 # cov
    #xlab = 'abundance'
    #bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 1, 2, 5, 10, 20, 50, 100,200,500,1000,2000,5000,1e4, 2e4, 5e4, 1e10]

    #choice = 1
    #xlab = 'len'
    #bins = [2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6]

    #choice = 2
    #xlab = 'eff_len'
    #bins = [2e1, 3e1, 4e1, 5e1, 6e1, 7e1, 8e1, 9e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]

    #choice = 4
    #xlab = 'tpm'
    #bins = [0, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 2e1, 3e1, 4e1, 5e1, 6e1, 7e1, 8e1, 9e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]

    choice = 6
    xlab = 'cov ratio'
    bins = [0, 1e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1]

    fp0 = [float(a[choice]) for a in itms if int(a[5])==0] #non fp
    fp1 = [float(a[choice]) for a in itms if int(a[5])==1] #fp
    #pdb.set_trace()
    print('fp ratio: %f'%(float(len(fp1))/(len(fp0)+len(fp1))))

    '''
    #modified fp ratio:
    fp0_ft = [float(a[choice]) for a in itms if int(a[5])==0 and float(a[choice])>2e-2] #non fp
    fp1_ft = [float(a[choice]) for a in itms if int(a[5])==1 and float(a[choice])>2e-2] #fp
    #pdb.set_trace()
    print('fp_ft ratio: %f'%(float(len(fp1_ft))/(len(fp0_ft)+len(fp1_ft))))
    '''

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

    pdb.set_trace()

    return

'''
usage
'''

## global config
#dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon_test/0.96/'
dataDir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyderSimChr15/refShannon/algo_output/'
prefix = 'reconstructed'
nJobs = 20

Trec = '%s/%s.fasta'%(dataDir, prefix)

reads = ['/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads_1.fq',
         '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads_2.fq']
RL = 101

dmpRes = 'dmp/dmp_stat_FP_vs_nonFP_maxSens.txt'

if __name__ == '__main__':

    #pdb.set_trace()

    if '1' in sys.argv:

        #generate kal abundance

        KAL_prefix = '%s/%s'%(dataDir, prefix)

        cmd = 'python analyze_simData_ROC.py --kallisto ' + \
              '-o %s '%KAL_prefix + \
              '--rec %s '%Trec + \
              '-r1 %s '%reads[0] + \
              '-r2 %s '%reads[1] + \
              '-t %d'%nJobs

        run_cmd(cmd)

    elif '2' in sys.argv:

        #Build hisat file
        run_cmd('hisat2-build %s %s/%s.hisat'%(Trec, dataDir, prefix))

        #Align
        run_cmd('hisat2 -p %d --no-spliced-alignment --no-discordant -q --fr -x %s/%s.hisat -1 %s -2 %s -S %s/%s.sam'%(nJobs, dataDir, prefix, reads[0], reads[1], dataDir, prefix))
        
        #Process SAM / BAM file to get depth information
        run_cmd('samtools view -@ %d -bS -f 0x2 %s/%s.sam > %s/%s.bam'%(nJobs, dataDir, prefix, dataDir, prefix))
        run_cmd('samtools sort -@ %d -o %s/%s.sorted.bam %s/%s.bam'%(nJobs, dataDir, prefix, dataDir, prefix))
        run_cmd('samtools depth %s/%s.sorted.bam > %s/%s.depth'%(dataDir, prefix, dataDir, prefix))

    elif '3' in sys.argv:

        #pdb.set_trace()

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

        
import pdb
import sys
from util import reverse_complement, from_fasta, from_fasta_lens
import subprocess

def run_cmd(cmd):
    print(cmd)
    subprocess.call(cmd, shell=True)

def find_comm_diff(res_A, res_B):

    dic = {} #key: tr name/id, val: [m1, r1, m2, r2, t, trec1, trec2]

    files = [res_A, res_B]
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                if line[0]=='#':
                    continue

                tokens = line.split()
                tr_id = tokens[0]
                dic[tr_id] = [0,0,0,0,0,0,0]
    #pdb.set_trace()

    for fi in range(2):
        with open(files[fi], 'r') as f:
            for line in f:
                if line[0]=='#': #last two lines
                    continue

                tokens = line.split()
                tr_id = tokens[0]
                m = int(tokens[2])
                t = int(tokens[3])
                r = int(tokens[4])

                dic[tr_id][fi*2]=m
                dic[tr_id][fi*2+1]=r
                dic[tr_id][4]=t
                dic[tr_id][5+fi]=tokens[1]

    return dic.items()

def dump(fp, vec, type='float'):
    f = open(fp, 'w')
    for v in vec:
        if type=='float':
            f.write('%f\n'%v)
        elif type=='str':
            f.write('%s\n'%v)
        elif type=='str int':
            f.write('%s\t%d\n'%(v[0], v[1]))
        elif type=='str int str':
            f.write('%s\t%d\t%s\n'%(v[0], v[1], v[2]))
        elif type=='str int str str':
            f.write('%s\t%d\t%s\t%s\n'%(v[0], v[1], v[2], v[3]))

    f.close()

def load(fp):
    x = []
    with open(fp, 'r') as f:
        for line in f:
            x.append(float(line.strip()))
    return x

def load_format(fp, format='str int', col=0):
    x = []
    with open(fp, 'r') as f:
        if format=='str int':
            x = [line.split()[col] for line in f]
        elif format=='float':
            x = [float(line.split()[0]) for line in f]
    #convert
    if format.split()[col]=='int':
        x = [int(a) for a in x]
    return x

def main_1():

    #dmp_dir = 'dmp/'

    #res1 = dmp_dir + '/snyder_0807b_log.txt'
    #res2 = dmp_dir + '/snyder_0831a_log.txt'

    #res1 = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_chr15/stringtie_log_0908.txt'
    #res2 = '/data1/shunfu1/ref_shannon_modi/data/snyder_0831a/chr15_0807b/chrs/chr15/_res_SF_condense_modi_zero_overlap_all_seg_refChr15/log.txt'

    #res1 = '/data1/shunfu1/ref_shannon_modi/data/ww_0910a/reconstr_log.txt'
    #res2 = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/_res_0726b_stringtie_f_0_c_0.001/stringtie_log.txt'

    res1 = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2/RefShannon/reconstructed_log.txt'
    res2 = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/star/RefShannon/reconstructed_log.txt'

    #list of [tr_id, [m1, r1, m2, r2, t, trec1, trec2]]
    itms = find_comm_diff(res1, res2)
    pdb.set_trace()

    tr_res1 = [[itm[0], itm[1][4], itm[1][5]] for itm in itms if float(itm[1][0])/itm[1][4]>=0.9]
    tr_res2 = [[itm[0], itm[1][4], itm[1][6]] for itm in itms if float(itm[1][2])/itm[1][4]>=0.9]
    pdb.set_trace()

    tr_res1.sort(key=lambda x:x[1])
    dump('dmp/res1.txt', tr_res1, 'str int')

    tr_res2.sort(key=lambda x:x[1])
    dump('dmp/res2.txt', tr_res2, 'str int')

    tr_res1_only = [[itm[0], itm[1][4], itm[1][5]] for itm in itms if float(itm[1][0])/itm[1][4]>=0.9 and float(itm[1][2])/itm[1][4]<0.9]
    dump('dmp/res1_only.txt', tr_res1_only, 'str int str')

    tr_res2_only = [[itm[0], itm[1][4], itm[1][5], itm[1][6]] for itm in itms if float(itm[1][2])/itm[1][4]>=0.9 and float(itm[1][0])/itm[1][4]<0.9]
    dump('dmp/res2_only.txt', tr_res2_only, 'str int str str')

    '''
    tr_stringtie = [float(itm[1][2])/itm[1][3] for itm in itms if float(itm[1][2])/itm[1][4]>=0.9 \
                                                                and float(itm[1][0])/itm[1][4]>=0.9]
    dump('dmp/stringtie_match_over_rec.txt', tr_stringtie, 'float')

    tr_shannon = [float(itm[1][0])/itm[1][1] for itm in itms if float(itm[1][0])/itm[1][4]>=0.9 \
                                                                and float(itm[1][2])/itm[1][4]>=0.9]
    dump('dmp/shannon_match_over_rec.txt', tr_shannon, 'float')
    '''

    '''
    tr_of_Shannon_only = [[itm[0], itm[1][4]] for itm in itms if float(itm[1][0])/itm[1][4]>=0.9 \
                                                              and float(itm[1][2])/itm[1][4]<0.9]    
    tr_of_Shannon_only.sort(key=lambda x:x[1])
    dump('dmp/tr_of_Shannon_only_snyder.txt', tr_of_Shannon_only, 'str int')

    tr_of_Stringtie_only = [[itm[0], itm[1][4]] for itm in itms if float(itm[1][2])/itm[1][4]>=0.9 \
                                                                and float(itm[1][0])/itm[1][4]<0.9]
    tr_of_Stringtie_only.sort(key=lambda x:x[1])
    dump('dmp/tr_of_Stringtie_only_snyder.txt', tr_of_Stringtie_only, 'str int')
    '''
    
    '''
    tr_of_Shannon_strict = [[itm[0], itm[1][4]] for itm in itms if float(itm[1][0])/itm[1][4]>=0.9 and float(itm[1][0])/itm[1][1]>=0.9]
    tr_of_Stringtie_strict = [[itm[0], itm[1][4]] for itm in itms if float(itm[1][2])/itm[1][4]>=0.9 and float(itm[1][2])/itm[1][3]>=0.9]

    #to_remove = [itm for itm in itms if itm[1][1]==0 or itm[1][3]==0 or itm[1][4]==0]
    #pdb.set_trace()
    tr_of_only_Stringtie_strict = [[itm[0], itm[1][4]] for itm in itms 
                          if itm[1][1]!=0 and itm[1][3]!=0 and itm[1][4]!=0 and
                          float(itm[1][0])/itm[1][4]<0.5  and float(itm[1][0])/itm[1][1]<0.5 and 
                          float(itm[1][2])/itm[1][4]>=0.9 and float(itm[1][2])/itm[1][3]>=0.9]

    tr_of_only_Stringtie_strict.sort(key=lambda x:x[1])
    pdb.set_trace()
    '''
    return

def main_2():

    import numpy as np
    import matplotlib.pyplot as plt

    dmp_dir = 'dmp/'

    fig, ax = plt.subplots()

    f_names = ['rA_of_trB.txt',
               'rB_of_trA.txt',
               'rA_of_noRec.txt',
               'rB_of_noRec.txt']

    for f_name in f_names:
        print(f_name)
        x = load(dmp_dir+f_name)
        hist, bins = np.histogram(x, bins=50)
        plt.plot(bins[1:], hist, marker='o', label=f_name)
    
    plt.legend(loc='upper center')
    ax.set_xlabel('bins')
    ax.set_ylabel('hist')
    #ax.set_title('relationship b/w Y and t, %d-iterations per t'%num_realization)
    #ax.bar(1, np.mean(y_hist, axis=0), 1, yerr=np.std(y_hist, axis=0), ecolor='black')
    #ax.set_ylabel('Y')
    plt.show()


    return



def gen_read_from_sam(sam_file, chrom, tr_stt, tr_end, out_dir, nLines):

    reads = {} # name, [seqA, seqB]

    i = 0
    j = 0
    T = nLines/100

    with open(sam_file) as f:
        for line in f:

            i += 1
            if i>T:
                i=0
                j+=1
                print('1st iter: %d perc processed'%j)

            if line[0] == '@': continue

            fields = line.split()

            name, flags, g_name, start, cigar, seq = (
                fields[0], int(fields[1]), fields[2],
                int(fields[3])-1, fields[5], fields[9])

            if g_name != chrom: continue
            if start<tr_stt or start>tr_end: continue
            if (flags>>8)&1==1: continue #2nd ali

            if name not in reads:
                reads[name] = ['', '']

            rev = bool((flags>>4)&1)
            if rev==0:
                read_seq = seq
            elif rev==1:
                read_seq = reverse_complement(seq)

            if (flags >> 6) & 1:#1st seg, A
                reads[name][0] = read_seq
            elif (flags >> 7) & 1:#2nd seq, B
                reads[name][1] = read_seq

            #pdb.set_trace()

    #pdb.set_trace()
    
    i=0
    j=0

    #try to include read who is not in reads while its pair is in reads
    with open(sam_file) as f:
        for line in f:

            i += 1
            if i>T:
                i=0
                j+=1
                print('2nd iter: %d perc processed'%j)

            if line[0] == '@': continue

            fields = line.split()

            name, flags, g_name, start, cigar, seq = (
                fields[0], int(fields[1]), fields[2],
                int(fields[3])-1, fields[5], fields[9])

            if g_name != chrom: continue
            #if start<tr_stt or start>tr_end: continue
            if (flags>>8)&1==1: continue #2nd ali

            if name not in reads: continue

            rev = bool((flags>>4)&1)
            if rev==0:
                read_seq = seq
            elif rev==1:
                read_seq = reverse_complement(seq)

            if (flags >> 6) & 1 and reads[name][0]=='':#1st seg, A
                #pdb.set_trace()
                reads[name][0] = read_seq
            elif (flags >> 7) & 1 and reads[name][1]=='':#2nd seq, B
                #pdb.set_trace()
                reads[name][1] = read_seq

            #pdb.set_trace()

    #pdb.set_trace()

    f_1 = open(out_dir + '/reads_1.fa', 'w')
    f_2 = open(out_dir + '/reads_2.fa', 'w')

    for name, seqs in reads.items():
        f_1.write('>%s\n%s\n'%(name, seqs[0]))
        f_2.write('>%s\n%s\n'%(name, seqs[1]))

    f_1.close()
    f_2.close()

    return

def main_3():

    '''
    case = 'case1'
    tr_id = 'ENST00000583635.1' #recovered by stringtie not refShannon    
    chrom = 'chr1'
    tr_stt = 59762658
    tr_end = 59781978
    '''

    #tr_id = 'ENST00000413017.2' #recovered by stringtie not refShannon    
    #chrom = 'chr21'
    #tr_stt = 34988593-1
    #tr_end = 35014028-1

    '''tr_id = 'ENST00000591854.1' #recovered by stringtie not refShannon    
    chrom = 'chr18'
    tr_stt = 55297534-1
    tr_end = 55334868-1'''

    #cp genome data
    #genome_file = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/genome/%s.fa'

    #extract reads from .sam
    #sam_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_all_sam/hits.sorted.sam'
    #sam_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/chr21/hits.sam'
    
    '''
    case = 'case_b_1'
    tr_id = 'ENST00000565181.1'
    chrom = 'chr15'
    tr_stt = 72765185
    tr_end = 72767509
    sam_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/%s/hits.sam'%chrom
    '''

    '''
    case = 'debug_int_chr7'
    tr_id = 'NA'
    chrom = 'chr7'
    tr_stt = 26766400
    tr_end = 26766600
    sam_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/%s/hits.sam'%chrom
    '''

    #case = 'debug_int2_modi_all_chrs'
    #tr_id 


    nLines = sum([1 for line in open(sam_file, 'r')])
    print('%s has %d lines'%(sam_file, nLines))
    out_dir = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/%s_%s'%(tr_id, case)
    cmd = 'mkdir -p %s'%out_dir
    run_cmd(cmd)

    gen_read_from_sam(sam_file, chrom, tr_stt, tr_end, out_dir, nLines)
    print('reads written to \n%s'%out_dir)

    return

def main_4():
    #from util import *

    #file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/' + 'reference.fasta'
    #fo = 'dmp/tr_len_snyder.txt'

    #file = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/reference.fasta'
    #fo = 'dmp/tr_len_ww.txt'

    file = '/data1/shunfu1/ref_shannon_modi/data/sim_0728b/reference.fasta'
    fo = 'dmp/tr_len_sim.txt'

    lens = from_fasta_lens(file)
    lens.sort(key=lambda x:x[1])

    dump(fo, lens, type='str int')

    pdb.set_trace()

    return

def main_5():

    import numpy as np
    import matplotlib.pyplot as plt

    #parameters
    dmp_dir = 'dmp/'
    '''f_names = [
               'tr_len_snyder.txt',
               'tr_of_Shannon_snyder.txt',
               'tr_of_Stringtie_snyder.txt'
               #'tr_len_ww.txt',
               #'tr_of_Shannon_ww.txt',
               #'tr_of_Stringtie_ww.txt'
               #'tr_len_sim.txt',
               #'tr_of_Shannon_sim.txt',
               #'tr_of_Stringtie_sim.txt'
               ]
    fmt = 'str int'
    col = 1
    nbins = 2000'''

    f_names = [
               'shannon_match_over_rec.txt',
               'stringtie_match_over_rec.txt'
               ]
    fmt = 'float'
    col = 0
    nbins = 100
    x_lab = 'match_len / rec_len'
    y_lab = 'hist'

    fig, ax = plt.subplots()
    for i in range(len(f_names)):
        f_name = f_names[len(f_names)-1-i]
        #f_name = f_names[i]
        print(f_name)
        x = load_format(dmp_dir+f_name, fmt, col)
        pdb.set_trace()
        hist, bins = np.histogram(x, bins=nbins)
        #plt.plot(bins[1:], hist, marker='o', label=f_name)
        plt.plot(bins[1:], hist, label=f_name)
    
    plt.legend(loc='upper right')
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    #ax.set_title('relationship b/w Y and t, %d-iterations per t'%num_realization)
    #ax.bar(1, np.mean(y_hist, axis=0), 1, yerr=np.std(y_hist, axis=0), ecolor='black')
    #ax.set_ylabel('Y')
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    plt.show()

    pdb.set_trace()
    return

def main_6():

    from tester import analyzer_blat_noExp

    dmp_dir = 'dmp/'

    choice = 'snyder'

    if choice == 'snyder':

        perf = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_chr15/stringtie_per.txt'
        logf = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_chr15/stringtie_log_0908.txt'
        analyzer_blat_noExp(perf,logf,'',0)

        '''
        loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_refShannon_0807b_chr15_only/'
        perf = loc + '/reconstr_per.txt'    
        logf = dmp_dir + '/snyder_0807b_log.txt'
        analyzer_blat_noExp(perf,logf,'',0)

        loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0821a/snyder_0729a_chr15/'
        perf = loc + '/reconstr_per.txt'
        logf = dmp_dir + '/snyder_int2_log.txt'
        analyzer_blat_noExp(perf,logf,'',0)
        '''

    elif choice == 'ww':
        loc = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/'

        perf_shannon = loc + '_res_0726c/reconstr_per.txt'
        logf_shannon = dmp_dir + 'refShannon_ww_log.txt'
        analyzer_blat_noExp(perf_shannon,logf_shannon,'',0)

        perf_stringtie = loc + '_res_0726b_stringtie_f_0_c_0.001/stringtie_per.txt'
        logf_stringtie = dmp_dir + 'stringtie_ww_log.txt'
        analyzer_blat_noExp(perf_stringtie,logf_stringtie,'',0)

    elif choice == 'sim':
        loc = '/data1/shunfu1/ref_shannon_modi/data/sim_0728b/'

        perf_shannon = loc + '_res_refShannon_0728b/reconstr_per.txt'
        logf_shannon = dmp_dir + 'refShannon_sim_log.txt'
        analyzer_blat_noExp(perf_shannon,logf_shannon,'',0)

        perf_stringtie = loc + '_res_stringtie_0728b_f_0_c_0.001/stringtie_per.txt'
        logf_stringtie = dmp_dir + 'stringtie_sim_log.txt'
        analyzer_blat_noExp(perf_stringtie,logf_stringtie,'',0)


    return

def load_snps(fp):
    snps = {}
    with open(fp,'r') as f:
        for line in f:
            fields = line.split()
            if 'snp' in fields:
                snps[int(fields[1])]=fields[2]
    return snps

def main_7():

    case = 'b_1'

    if case==1:
        tr_id = 'ENST00000583635.1'
        chrom = 'chr1'
        exons = [[59762658 ,59762822], [59781264,59781408],  [59781818,59781978]]#1 based, ref gtf
        genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000583635.1_case1/genome/chr1.fa'
        k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000583635.1_case1/chrs/chr1/algo_input/k1mer.dict'

    elif case==2:
        tr_id = 'ENST00000413017.2'#'ENST00000583635.1'#'ENST00000583635.1' #recovered by stringtie not refShannon    
        chrom = 'chr21' #'chr1' #case 2
        exons = [[34988593,        34988909],
                 [34989012  ,      34989056],
                 [34994302  ,      34994374],
                 [34996989   ,     34997066],
                 [35003792   ,     35003863],
                 [35013987    ,    35014028]]#1-based
        genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/genome/chr21.fa'
        #case2
        #old
        #k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/chr21/algo_input/k1mer.dict_org'
        #k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/chr21/algo_input/k1mer.dict'
        #k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs_old_sam_full/chr21/algo_input/k1mer.dict_org'
        #old modi
        #k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs_old_sam_full/chr21/algo_input/modi_k1mer.dict'
        #old --> re-run, ext ec off, minweight 3
        k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs/chr21/algo_input/k1mer.dict'
        #new
        #k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs_new_sam/chr21/algo_input/k1mer.dict_org'
    
    elif case==3:
        tr_id = 'ENST00000591854.1'
        chrom = 'corrected' #'corrected' #'chr18'
        exons = [ [55297534,        55297750],
                  [55308791,        55309284],
                  [55334474,        55334868]]
        #genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000591854.1_case3/genome/chr18.fa' #chrs/chr18/%s.fa'%chrom #genome/chr18.fa'
        genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000591854.1_case3/chrs/chr18/%s.fa'%chrom
        k1mer_file =  '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000591854.1_case3/chrs/chr18/algo_input/k1mer.dict'

    elif case=='b_1':
        tr_id = 'ENST00000565181.1'
        chrom = 'chr15' #'corrected' #'chr18'
        exons = [[72765185, 72767509]]

        genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/ENST00000565181.1_case_b_1/genome/chr15.fa'
        #0807b/large
        #k1mer_file =  '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/chr15/algo_input/k1mer.dict'
        #0807b/small
        k1mer_file =  '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/ENST00000565181.1_case_b_1/chrs_0807b/chr15/algo_input/k1mer.dict'
        #0821a/large
        #k1mer_file =  '/data1/shunfu1/ref_shannon_modi/data/snyder_0821a/snyder_0729a_chr15/chrs/chr15/algo_input/k1mer.dict'
        #0821a/small
        #k1mer_file =  '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/ENST00000565181.1_case_b_1/chrs_0821a/chr15/algo_input/k1mer.dict'



    snps = {} #load_snps('dmp/dmp_error_correction.txt')
    print('snps loaded')

    print('k1mer file:\n%s\n'%k1mer_file)

    genome = ''
    for name, seq in from_fasta(genome_file):
        if name==chrom:
            genome = seq
            break
    print('len of genome: %d'%len(genome))

    seqs = [genome[e[0]-1:e[1]] for e in exons]#exon 1-based
    tr_seq = ''.join(seqs)
    print('len of tr_seq: %d\n%s'%(len(tr_seq), tr_seq))

    K1=25
    dic={} #key - k1mer val - [pos,...,]
    for i in range(len(tr_seq)-K1+1):
        k1mer = tr_seq[i:i+K1]
        if k1mer in dic:
            dic[k1mer].append(i)
        else:
            dic[k1mer]=[i]
    print('# of k1mers:%d'%len(dic))

    dic2={} #key - k1mer val - cnt
    for k1mer in dic:
        dic2[k1mer]=0
    #pdb.set_trace()

    with open(k1mer_file) as f:
        for line in f:
            k1mer, cnt = line.split()[:2]
            if k1mer in dic2:
                dic2[k1mer]=int(cnt)
    #pdb.set_trace()

    tr_seq2 = list(tr_seq)
    cnt_0=0
    for k1mer, cnt in dic2.items():
        if cnt==0:
            #pdb.set_trace()
            cnt_0+=1
            if len(dic[k1mer])>1:
                pdb.set_trace()
            pos = dic[k1mer][0]
            genome_pos = exons[0][0]-1+pos
            if genome_pos in snps:
                tr_seq2[pos]='*'
            else:
                tr_seq2[pos]='-'
    tr_seq2 = ''.join(tr_seq2)
    #print('tr covered by available k1mers:%s'%tr_seq2)
    print('# of k1mers with 0 counts: %d\n'%cnt_0)

    acc_len = 0
    for i in range(len(exons)):
        exon = exons[i]
        exon_l = exon[1]-exon[0]+1
        exon_st = tr_seq2[acc_len:acc_len+exon_l]
        acc_len += exon_l
        print('exon %d: [%d, %d) (len=%d, acc_len=%d)\n%s\n'%(i, exon[0]-1, exon[1], exon_l, acc_len, exon_st))

    pdb.set_trace()

    return

def main_8():

    file = 'dmp/tr_of_Stringtie_only_snyder.txt'
    print('file:%s'%file)

    tr_ids = load_format(file, format='str int', col=0)
    print('# of tr:%d'%(tr_ids))



    return

def get_alignment_info(reads_dic, sam_file, chrom):
    #reads_dic: key - read name ; val - list of [[pos of A, +/- of A, cigar], [pos of B, +/- of B, cigar]]

    i = 0
    j = 0
    nLines = sum([1 for line in open(sam_file)])
    T = nLines/100

    with open(sam_file) as f:
        for line in f:

            i += 1
            if i>T:
                i=0
                j+=1
                print('%d perc processed'%j)

            if line[0] == '@': continue

            fields = line.split()

            name, flags, g_name, start, cigar, seq = (
                fields[0], int(fields[1]), fields[2],
                int(fields[3])-1, fields[5], fields[9])

            if g_name != chrom: continue
            #if start<tr_stt or start>tr_end: continue
            #if (flags>>8)&1==1: continue #2nd ali

            if name not in reads_dic: continue

            #pdb.set_trace()

            rev = bool((flags>>4)&1)
            sign = '+'
            if rev==1:
                sign = '-'
            
            if (flags >> 6) & 1:#1st seg, A
                itm = 0
            elif (flags >> 7) & 1:#2nd seq, B
                itm = 1
            else:
                print('unexpected: %d'%flags)

            reads_dic[name][-1][itm]=[start, sign, cigar]
            if reads_dic[name][-1][0] is not None and reads_dic[name][-1][1] is not None:
                reads_dic[name].append([None, None])

    return reads_dic

def dmp_reads_dic(reads_dic, fp):
    with open(fp, 'w') as f:
        for name, alignments in reads_dic.items():
            st = name + '\n'
            for alignment in alignments:
                if alignment[0] is not None:
                    st += '%d,%s,%s'% \
                         (alignment[0][0],
                          alignment[0][1],
                          alignment[0][2])
                if alignment[1] is not None:
                    st += '//%d,%s,%s'% \
                         (alignment[1][0],
                          alignment[1][1],
                          alignment[1][2])
                st += '\n'
            st += '\n'
            f.write(st)

def main_9():

    #loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000583635.1_case1/'
    loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/'
    reads_files = [loc+'reads_1.fa',
                   loc+'reads_2.fa']
    chrom = 'chr21'
    old_sam = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/%s/hits.sam'%chrom
    new_sam = loc + 'hits.sam'

    dmp_reads_old = 'dmp/case2_reads_old.txt'
    dmp_reads_new = 'dmp/case2_reads_new.txt'

    reads_old = {} #key - read name val - list of [[pos of A, +/- of A, cigar], [pos of B, +/- of B, cigar]]
    reads_new = {}
    with open(reads_files[0]) as f:
        for line in f:
            if line[0]=='>':
                reads_old[line.split()[0][1:]]=[[None, None]]
                reads_new[line.split()[0][1:]]=[[None, None]]
    #pdb.set_trace()

    #old sam
    reads_old = get_alignment_info(reads_old, old_sam, chrom)
    #pdb.set_trace()

    #new sam
    reads_new = get_alignment_info(reads_new, new_sam, chrom)
    #pdb.set_trace()

    #export reads_old and new
    dmp_reads_dic(reads_old, dmp_reads_old)
    dmp_reads_dic(reads_new, dmp_reads_new)
    #pdb.set_trace()

    return

def main_10():

    #loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000583635.1_case1/'
    loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/'
    reads_files = [loc+'reads_1.fa',
                   loc+'reads_2.fa']
    chrom = 'chr21'
    old_sam = '/data1/shunfu1/ref_shannon_modi/data/snyder_0729a/chrs_0807b/%s/hits.sam'%chrom
    new_sam = loc + 'hits.sam'

    filtered_old_sam = loc + 'hits_old.sam'

    reads = []
    with open(reads_files[0]) as f:
        for line in f:
            if line[0]=='>':
                reads.append(line.split()[0][1:])
    #pdb.set_trace()

    i = 0
    j = 0
    nLines = sum([1 for line in open(old_sam)])
    T = nLines/100

    fo = open(filtered_old_sam, 'w')
    with open(old_sam) as f:
        for line in f:

            i += 1
            if i>T:
                i=0
                j+=1
                print('%d perc processed'%j)

            if line[0] == '@': continue

            fields = line.split()

            name, flags, g_name, start, cigar, seq = (
                fields[0], int(fields[1]), fields[2],
                int(fields[3])-1, fields[5], fields[9])

            if g_name != chrom: continue
            if name not in reads: continue

            fo.write(line)
    fo.close()
    print('%s written'%filtered_old_sam)

    return

def main_11():

    k1mer_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs_old_sam_full/chr21/algo_input/k1mer.dict_org'
    k1mer_file_modi = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/chrs_old_sam_full/chr21/algo_input/modi_k1mer.dict'

    from mb_sf_denovo.ext_corr import run_correction

    run_correction(k1mer_file, k1mer_file_modi, 3, 75, True)
    print('%s written'%k1mer_file_modi)

    #pdb.set_trace()
    
    

    return

def main_12():

    from util import reverse_complement

    f_c0 = 'dmp/k1mer_c_off_L_1'
    f_c1 = 'dmp/k1mer_c_on_L_1'

    d_c0 = {}
    d_c1 = {}

    d_c0_merge = {}

    with open(f_c0, 'r') as f:
        for line in f:
            fields = line.split()
            d_c0[fields[0]]=int(fields[1])
    list_c0 = d_c0.items()
    list_c0.sort(key=lambda x:x[0])
    print('%s loaded, # of k1mers=%d'%(f_c0, len(d_c0)))
    #pdb.set_trace()

    with open(f_c1, 'r') as f:
        for line in f:
            fields = line.split()
            d_c1[fields[0]]=int(fields[1])
    list_c1 = d_c1.items()
    list_c1.sort(key=lambda x:x[0])
    print('%s loaded, # of k1mers=%d'%(f_c1, len(d_c1)))
    #pdb.set_trace()

    for k1mer, cnt in list_c0:

        rc = reverse_complement(k1mer)

        if k1mer in d_c0_merge:
            d_c0_merge[k1mer] += cnt
        elif rc in d_c0_merge:
            d_c0_merge[rc] += cnt
        else:
            d_c0_merge[k1mer] = cnt
    print('merge k1mer and its rev comp into d_c0_merge')
    print('# of k1mers=%d'%(len(d_c0_merge)))
    #pdb.set_trace()

    #list_c1 = d_c1.items()
    #list_c1.sort(key=lambda x:x[0])
    dump('dmp/c1_sorted', list_c1, type='str int')

    list_c0_merge = d_c0_merge.items()
    list_c0_merge.sort(key=lambda x:x[0])
    dump('dmp/c0_merge_sorted', list_c0_merge, type='str int')
    pdb.set_trace()
    return

def main_13():

    tr_id = 'ENST00000413017.2'#'ENST00000583635.1'#'ENST00000583635.1' #recovered by stringtie not refShannon    
    chrom = 'chr21' #'chr1'
    #case 2
    exons = [[34988593,        34988909],
             [34989012  ,      34989056],
             [34994302  ,      34994374],
             [34996989   ,     34997066],
             [35003792   ,     35003863],
             [35013987    ,    35014028]]#1-based
    genome_file = '/data1/shunfu1/ref_shannon_modi/data/snyder_0809a/ENST00000413017.2_case2/genome/chr21.fa'

    snps = load_snps('dmp/dmp_error_correction.txt')
    print('snps loaded')
    pdb.set_trace()

    k1mer_file = 'dmp/k1mer_c_on_L_1'

    c1_dic = {}
    with open(k1mer_file, 'r') as f:
        for line in f:
            fields = line.split()
            c1_dic[fields[0]]=int(fields[1])

    print('k1mer file:\n%s\n'%k1mer_file)

    genome = ''
    for name, seq in from_fasta(genome_file):
        if name==chrom:
            genome = seq
            break
    print('len of genome: %d'%len(genome))

    seqs = [genome[e[0]-1:e[1]] for e in exons]#1-based
    tr_seq = ''.join(seqs)
    print('len of tr_seq: %d'%len(tr_seq))

    K1=25
    dic={} #key - k1mer val - [pos,...,]
    for i in range(len(tr_seq)-K1+1):
        k1mer = tr_seq[i:i+K1]
        if k1mer in dic:
            dic[k1mer].append(i)
        else:
            dic[k1mer]=[i]
    print('# of k1mers:%d'%len(dic))

    dic2={} #key - k1mer val - cnt
    for k1mer in dic:
        dic2[k1mer]=0
    #pdb.set_trace()

    with open(k1mer_file) as f:
        for line in f:
            k1mer, cnt = line.split()[:2]
            if k1mer in dic2:
                dic2[k1mer]=int(cnt)
    #pdb.set_trace()

    k1mer_cnt0 = set()
    tr_seq2 = list(tr_seq)
    cnt_0=0
    for k1mer, cnt in dic2.items():
        if cnt==0:
            cnt_0+=1
            pos = dic[k1mer][0]
            genome_pos = exons[0][0]-1+pos
            if genome_pos in snps:
                tr_seq2[pos]='*'
            elif reverse_complement(k1mer) in c1_dic:
                tr_seq2[pos]='+'
                k1mer_cnt0.add(k1mer)
            else:
                tr_seq2[pos]='-'
    tr_seq2 = ''.join(tr_seq2)
    #print('tr covered by available k1mers:%s'%tr_seq2)
    print('# of k1mers with 0 counts: %d\n'%cnt_0)

    acc_len = 0
    for i in range(len(exons)):
        exon = exons[i]
        exon_l = exon[1]-exon[0]+1
        exon_st = tr_seq2[acc_len:acc_len+exon_l]
        acc_len += exon_l
        print('exon %d: [%d, %d) (len=%d, acc_len=%d)\n%s\n'%(i, exon[0]-1, exon[1], exon_l, acc_len, exon_st))

    pdb.set_trace()

    return

if __name__ == '__main__':

    args = sys.argv

    choice = int(args[1])

    if choice==1: #stats
        main_1()
    elif choice==2: #histogram
        main_2()
    elif choice==3: #find genome regions and reads for certain tr
        main_3()
    elif choice==4:
        main_4() #histogram of tr len of reference
    elif choice==5:
        main_5() #draw histogram of tr len of reference
    elif choice==6: #re-gen log files with tr_len included
        main_6()
    elif choice==7: #check actual kmer weights of an expected tr
        main_7()
    elif choice==8: #check average read coverage of tr recovered only by stringtie
        main_8()
    elif choice==9: #compare newly alignments of extracted reads with original alignments of same reads
        main_9()
    elif choice==10: #filter old sam based on extracted reads
        main_10()
    elif choice==11: #ext ec modi
        main_11()
    elif choice==12: #convert k1mer.dict_org from jellyfish -C off to -C on res
        main_12()
    elif choice==13:
        #based on main_7

        #when we find the target isoform has a 0-cnt k1mer in k1mer.dict_org 
        #(jellyfish -C), we check if its reverse complementary is in k1mer.dict_org,
        # and get its count

        main_13()

from util import *
import numpy as np

def check_sam_file(sam_file, target=''):

    flags_stat = {}
    cigar_stat = {'S':[], 'M':[], 'N':[], 'D':[], 'I':[]}

    num_lines = sum([1 for l in open(sam_file)]); print('%d lines in %s'%(num_lines, sam_file))

    i = 0; j=0; T=num_lines/100;

    with open(sam_file, 'r') as f:

        for line in f:
            i+=1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of sam processed'%j); sys.stdout.flush()

            if line[0]=='@': continue

            fields = line.split()
            name, flags, g_name, start, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]

            if flags in flags_stat:
                flags_stat[flags]+=1
            else:
                flags_stat[flags]=1

            if target != '' and g_name != target:
                pdb.set_trace()
                continue

            if ((flags >> 2) & 1):
                pdb.set_trace()
                continue #unmapped

            [valid_cigar, _, modi_blocks, splice_sites, cigar_stat0] = cigar2blocks6(start, cigar)
            if valid_cigar==0:
                pdb.set_trace()
                continue #not M,N,S,I,D

            keys = ['S', 'M', 'N', 'D', 'I']
            for key in keys:
                if key in cigar_stat0:
                    cigar_stat[key]+=cigar_stat0[key]

            #for sp_l, sp_r in splice_sites:
            #   pdb.set_trace()
            #   splice_starts[sp_l] += 1
            #   splice_ends[sp_r] += 1

            #for start, length in modi_blocks:
            #    for wp in xrange(start, start + length):
            #        #if wp == 30925692: pdb.set_trace()
            #        weights[wp] += 1

    print('')
    #pdb.set_trace()
    return flags_stat, cigar_stat

def PrintFlagsStat(flags_stat):

    flags_stat_of = '%s/accepted_hits_flags_stat.txt'%out_dir

    with open(flags_stat_of, 'w') as f:
    
        for k, v in flags_stat.items():
            bit_stream = '{0:08b}'.format(k)
            st = '%d\t%s\t%d\n'%(k, bit_stream, v)
            f.write(st)
            print(st)
    print('%s written'%flags_stat_of)
    return

def PrintCigarStat(cigar_stat):

    for c, n_list in cigar_stat.items():
        if len(n_list)>0:
            if c=='M':
                bins = [i for i in xrange(1,110,10)]
            elif c=='N':
                bins = [1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e10]
            elif c=='S':
                bins = [i for i in xrange(1,110,10)]
            elif c=='I':
                bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]
            elif c=='D':
                bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]
            [hist, bins] = np.histogram(n_list, bins=bins)
            print(c)
            print(hist)
            print(bins)
    return

sam_file = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2_unsorted/accepted_hits.sam'
#sam_file = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2/accepted_hits.sam'
#sam_file = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/star/hits.sam'
out_dir = parent_dir(sam_file)

flags_stat, cigar_stat = check_sam_file(sam_file, 'chr15')
#pdb.set_trace()

'''
{385: 26809, 387: 177, 179: 22, 129: 10317, 137: 596352, 371: 68, 337: 183723, 115: 22, 401: 181619, 147: 628523, 433: 26176, 409: 74660, 153: 655133, 393: 79183, 161: 1337104, 163: 627148, 65: 10317, 369: 26176, 355: 103179, 419: 100506, 177: 11021, 435: 68, 417: 183723, 329: 79694, 321: 26809, 323: 177, 353: 181619, 73: 596713, 403: 103179, 81: 1337104, 339: 100506, 89: 653644, 97: 1337745, 67: 45, 99: 628523, 131: 45, 145: 1337745, 113: 11021, 83: 627148, 345: 74982}
'''
PrintFlagsStat(flags_stat)
#pdb.set_trace()

'''
tophat2 unsorted:

['I', 'S', 'M', 'D', 'N']
[101354, 0, 16553921, 13114, 4480728]

I - 100071 1s and 1283 2s
M - hist array([              902739,  944910, 1099676,  882678,  899044,  870036,  860099,  822793,  775344, 8496602])
    bins array([      1.,     11.,     21.,     31.,     41.,     51.,     61.,     71.,     81.,     91.,    101.])
D - 12638 1s and 476 2s
N - 
    hist array([             4222185,      0447,       45455,     7885,      627,        3736,       1760,      176845,        545,     1243])
    bins array([  5.2e1,   4.92882e4,   9.852e4,    1.4776e5,   1.96e5,   2.46e5,    2.9546e5,   3.4470e5,   3.93941e5,    4.431e5,   4.92e5])

'''
#PrintCigarStat(cigar_stat)
#pdb.set_trace()

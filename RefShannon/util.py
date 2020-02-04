from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import subprocess, re, time, sys, os, shutil, pdb

def run_cmd(cmd, shell=True):
    #print(cmd)
    subprocess.call(cmd, shell=shell)

class Clock:
    def __init__(self):
        self.last = time.time()
    def time(self):
        cur = time.time()
        elapsed = cur - self.last
        self.last = cur
        return elapsed
    def asctime(self):
        return str(time.asctime())

class Counter(dict):
    def __missing__(self, key):
        return float()
    def max(self):
        max_key, max_value = None, -1
        for k, v in self.items():
            if v >= max_value:
                max_key, max_value = k, v
        return max_key
    def total(self):
        return sum(v for _, v in self.items())

class Alignment:
    def __init__(self, g_name, blocks, rev, seq, weight, splice_sites=[]):
        self.g_name = g_name
        self.blocks = blocks
        self.weight = weight
        #self.seq = seq
        self.reversed = rev
        self.splice_sites = splice_sites
    def start(self):
        return self.blocks[0][0]
    def end(self):
        start, length = self.blocks[-1]
        return start + length - 1

def from_fasta(filename):
    sequences = {} #[]
    seq = []
    next_name = '_'
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                #sequences.append((next_name, ''.join(seq)))
                sequences[next_name]=''.join(seq)
                seq = []
                next_name = line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
    #sequences.append((next_name, ''.join(seq)))
    sequences[next_name]=''.join(seq)
    del sequences['_']
    return sequences #[1:]

def combineSeperateLines(fasta_src, fasta_dst): #post process stringtie outputs
    seq = []
    next_name = '_'
    with open(fasta_src, 'r') as src, \
         open(fasta_dst, 'w') as dst:
        for line in src:
            if line[0] == '>':
                seq_str = ''.join(seq)
                if next_name != '_':
                    dst.write('>%s\n%s\n'%(next_name, seq_str))
                seq = []
                next_name = line.strip()[1:] #line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
        seq_str = ''.join(seq)
        if next_name != '_':
            dst.write('>%s\n%s\n'%(next_name, seq_str))
    return

def enforce_unique_tr_name(target_fasta):

    tmp_fasta = target_fasta + '.tmp'

    seq = []
    next_name_line = '_'
    names = {}

    with open(target_fasta, 'r') as src, \
         open(tmp_fasta, 'w') as dst:

        for line in src:
            if line[0] == '>':
                seq_str = ''.join(seq)
                if next_name_line != '_':
                    dst.write('>%s\n%s\n'%(next_name_line, seq_str))
                seq = []
                next_name_line = line.strip()[1:] #line.split()[0][1:]
                
                part1 = next_name_line.split()[0]
                part2 = next_name_line.split()[1:]
                if part1 in names:
                    names[part1]+=1
                else:
                    names[part1]=1

                next_name_line = '%s.%d'%(part1, names[part1])+' '+''.join(part2)
            else:
                seq.append(line.strip().upper())
        seq_str = ''.join(seq)
        if next_name_line != '_':
            dst.write('>%s\n%s\n'%(next_name_line, seq_str))

    cmd = 'mv %s %s'%(tmp_fasta, target_fasta)
    run_cmd(cmd)               
    
    return

def to_fasta(filename, seq, seq_name):
    seqs_to_fasta(filename, [(seq_name, seq)])

def seqs_to_fasta(filename, seqs):
    LINE_LENGTH = 50
    with open(filename, 'w') as f:
        for seq_name, seq in seqs:
            f.write(">{}\n".format(seq_name))
            for i in xrange(0, len(seq), LINE_LENGTH):
                f.write("{}\n".format(seq[i:i+LINE_LENGTH]))

def cigar2blocks5(ref_i, cigar):

    num_pattern = re.compile('[\d]+')
    chr_pattern = re.compile('[^\d]+')
    s = num_pattern.findall(cigar)

    splice_n = [int(s[i]) for i in range(len(s))]
    splice_c = chr_pattern.findall(cigar)

    modi_seq = ''
    blocks = []
    splice_sites = [] # (l,r) if read has 4M50N94M, l is left pos of read (inclusive) before 50N
                      #                             r is right pos of read (inclusive) after 50N

    if splice_c[0]=='S':
        #L = splice_n[0]
        splice_c = splice_c[1:]
        splice_n = splice_n[1:]
        #read = read[L:]

    if splice_c[-1]=='S':
        #L = splice_n[-1]
        splice_c = splice_c[:-1]
        splice_n = splice_n[:-1]
        #read = read[:-L]

    for i in range(len(splice_c)):
        c = splice_c[i]
        n = splice_n[i]
        if c=='M':
            #read_part = read[read_i:read_i+n]
            #modi_seq += read_part
            #blocks.append([ref_i, ref_i+n-1]) #blocks.append((start, length))
            blocks.append((ref_i, n))
            #ref_part = ref[ref_i:ref_i+n]
            #blocks.append([ref_i, ref_i+n-1, read_part, ref_part])
            #read_i += n
            ref_i += n
        elif c=='N':
            splice_left = ref_i-1
            splice_right = ref_i+n
            splice_sites.append((splice_left, splice_right)) #inclusive
            #pdb.set_trace()
            #read_i = read_i
            ref_i += n
        elif c=='D':
            #read_i = read_i
            ref_i += n
        elif c=='I':
            #read_i += n
            ref_i = ref_i
        else:
            return [0, '', [], []]
    return [1, modi_seq, blocks, splice_sites]

#based on cigar2blocks5, offer cigar pattern stats
def cigar2blocks6(ref_i, cigar):

    num_pattern = re.compile('[\d]+')
    chr_pattern = re.compile('[^\d]+')
    s = num_pattern.findall(cigar)

    splice_n = [int(s[i]) for i in range(len(s))]
    splice_c = chr_pattern.findall(cigar)

    modi_seq = ''
    blocks = []
    splice_sites = [] # (l,r) if read has 4M50N94M, l is left pos of read (inclusive) before 50N
                      #                             r is right pos of read (inclusive) after 50N

    cigar_stat = {} #key - S/M/N/D/I val - list of val

    if splice_c[0]=='S':
        #L = splice_n[0]
        splice_c = splice_c[1:]
        splice_n = splice_n[1:]
        #read = read[L:]

        cigar_stat['S']=[splice_n[0]]

    if splice_c[-1]=='S':
        #L = splice_n[-1]
        splice_c = splice_c[:-1]
        splice_n = splice_n[:-1]
        #read = read[:-L]

        if 'S' in cigar_stat:
            cigar_stat['S'].append(splice_n[-1])
        else:
            cigar_stat['S']=[splice_n[-1]]

    for i in range(len(splice_c)):
        c = splice_c[i]
        n = splice_n[i]
        if c=='M':
            #read_part = read[read_i:read_i+n]
            #modi_seq += read_part
            #blocks.append([ref_i, ref_i+n-1]) #blocks.append((start, length))
            blocks.append((ref_i, n))
            #ref_part = ref[ref_i:ref_i+n]
            #blocks.append([ref_i, ref_i+n-1, read_part, ref_part])
            #read_i += n
            ref_i += n

            if 'M' in cigar_stat:
                cigar_stat['M'].append(n)
            else:
                cigar_stat['M']=[n]
        elif c=='N':
            splice_left = ref_i-1
            splice_right = ref_i+n
            splice_sites.append((splice_left, splice_right)) #inclusive
            #pdb.set_trace()
            #read_i = read_i
            ref_i += n

            if 'N' in cigar_stat:
                cigar_stat['N'].append(n)
            else:
                cigar_stat['N']=[n]

        elif c=='D':
            #read_i = read_i
            ref_i += n

            if 'D' in cigar_stat:
                cigar_stat['D'].append(n)
            else:
                cigar_stat['D']=[n]

        elif c=='I':
            #read_i += n
            ref_i = ref_i

            if 'I' in cigar_stat:
                cigar_stat['I'].append(n)
            else:
                cigar_stat['I']=[n]

        else:
            return [0, '', [], [], cigar_stat]
    return [1, modi_seq, blocks, splice_sites, cigar_stat]

def parent_dir(dir_path):
    dir_path = dir_path.split('/')
    dir_path = [itm for itm in dir_path if itm != '']
    return '/'+'/'.join(dir_path[:-1])+'/'


def do_copy_chrs(args):

    chrs_dir_src = args[args.index('-I')+1]
    chrs_dir_dst = args[args.index('-O')+1]

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir_src)

    cnt = 0
    for target in target_list:
        if 'chr' in target and os.path.exists('%s/%s/hits.sam'%(chrs_dir_src, target))==True:
            cnt += 1
            print('%d %s/%s'%(cnt, chrs_dir_src, target))
            run_cmd('mkdir -p %s/%s/'%(chrs_dir_dst, target))
            shutil.copyfile('%s/%s/hits.sam'%(chrs_dir_src, target), \
                            '%s/%s/hits.sam'%(chrs_dir_dst, target))

    return

def do_splitMultiFasta(args):

    fasta_file = args[args.index('-i')+1]
    out_dir = args[args.index('-O')+1]
    run_cmd('mkdir -p %s'%out_dir)

    cnt = 0
    for name, seq in from_fasta(fasta_file).items():
        cnt += 1
        to_fasta('%s/%s.fa'%(out_dir, name), seq, name)
        print('%5d written: %s/%s.fa'%(cnt, out_dir, name))

    return

def do_splitRead(args):

    read_file = args[args.index('-r')+1]
    out_file = args[args.index('-o')+1]
    #stt = int(args[args.index('-s')+1])
    L = int(args[args.index('-l')+1])

    num_lines = sum([1 for line in open(read_file)])
    print('%d lines at%s'%(num_lines, read_file))
    i=0; j=0; T=num_lines/100;

    with open(read_file, 'r') as inf, \
         open(out_file, 'w') as of:

         for line in inf:
            i += 1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed'%j); sys.stdout.flush()

            #if i>10: break #test purpose

            if line.strip()=='': continue

            if line[0]=='@':
                st = line.split()[0]
                of.write(st+'\n')
            elif line[0]=='+':
                of.write('+\n')
            else:
                if L>0:
                    st = line.strip()[0:0+L]
                else:
                    st = line.strip()[L:]
                    #pdb.set_trace()
                of.write(st+'\n')
    print('%s written'%(out_file))
    return

'''
extract reads from sam file
'''
def reverse_complement(s):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return ''.join(comp[c] for c in s[::-1])
'''
def gen_read_from_sam(sam_file, chrom, out_dir):

    reads = {} # name, [seqA, seqB]

    i = 0
    j = 0
    nLines = sum([1 for l in open(sam_file)]); print('%d lines in %s'%(nLines, sam_file))
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
            #if start<tr_stt or start>tr_end: continue
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
'''

#DEBUG from 0922a
#get lens of tr of fasta file
def from_fasta_lens(filename):
    sequences = []
    seq = []
    next_name = '_'
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                #sequences.append((next_name, ''.join(seq)))
                sequences.append((next_name, len(''.join(seq))))
                seq = []
                next_name = line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
    #sequences.append((next_name, ''.join(seq)))
    sequences.append((next_name, len(''.join(seq))))
    return sequences[1:]

'''
usage:

#copy chrs_dir_src/chr_i/hits.sam to chrs_dir_dst/chr_i/hits.sam

python util.py --copy -I chrs_dir_src -O chrs_dir_dst [-chrs chr_a[,chr_b,...]]

#split multi fasta file to seperate ones (e.g. hg19.fa to chr1.fa etc)

python util.py --splitMultiFasta -i fasta_file -O out_dir

#split read file into one with shorter lengths (e.g. take first 50 bases (len=50) or last 50 bases (len=-50))

python util.py --splitRead -r read_file -o read_file_out -l len

'''
if __name__ == '__main__':
    args = sys.argv
    if '--copy' in args:
        do_copy_chrs(args)
    if '--splitMultiFasta' in args:
        do_splitMultiFasta(args)
    if '--splitRead' in args:
        do_splitRead(args)
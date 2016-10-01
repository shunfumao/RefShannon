import subprocess, re, time, sys, os, shutil

def run_cmd(cmd):
    #print(cmd)
    subprocess.call(cmd, shell=True)

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

def combineSeperateLines(fasta_src, fasta_dst):
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
                next_name = line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
        seq_str = ''.join(seq)
        if next_name != '_':
            dst.write('>%s\n%s\n'%(next_name, seq_str))
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
        L = splice_n[0]
        splice_c = splice_c[1:]
        splice_n = splice_n[1:]
        #read = read[L:]

    if splice_c[-1]=='S':
        L = splice_n[-1]
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

'''
usage:

#copy chrs_dir_src/chr_i/hits.sam to chrs_dir_dst/chr_i/hits.sam

python util.py --copy -I chrs_dir_src -O chrs_dir_dst [-chrs chr_a[,chr_b,...]]

#split multi fasta file to seperate ones (e.g. hg19.fa to chr1.fa etc)

python util.py --splitMultiFasta -i fasta_file -O out_dir

'''
if __name__ == '__main__':
    args = sys.argv
    if '--copy' in args:
        do_copy_chrs(args)
    if '--splitMultiFasta' in args:
        do_splitMultiFasta(args)
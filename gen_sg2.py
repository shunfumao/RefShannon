import sys, os, pdb, subprocess, copy
from util import *

'''
memory debug
from memory_profiler import profile
memory_fp = open('memory_profiler.log', 'w+')
'''

ALLOWED_HOLE = 10 # <= ALLOWED_HOLE, interpolate directly
MAX_HOLE = 200 # >= MAX_HOLE, no interpolate

'''

sam --> splice graph

usage:

python gen_sg2.py -i sam_file -g genome_file [-O out_dir] [-paired] [-target chr]


'''

def parse_args_gen_splice_graph2():

    args = sys.argv

    args_error = False 

    if '-i' in args:
        sam_file = args[args.index('-i')+1]
    else:
        args_error = True

    if '-g' in args:
        genome_file = args[args.index('-g')+1]
    else:
        args_error = True

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
        run_cmd('mkdir -p %s'%out_dir)
    else:
        out_dir = parent_dir(sam_file)+'/' #default is .sam folder

    if '-paired' in args:
        paired = 1
    else:
        paired = 0

    if '-target' in args:
        target = args[args.index('-target')+1]
    else:
        target = ''

    return [args_error, sam_file, genome_file, out_dir, paired, target]

def gen_splice_graph2():

    [args_error, sam_file, genome_file, subdir, paired, target] = parse_args_gen_splice_graph2()
    if args_error == True:
        print('gen_sg2 arguments error')
        return

    #parameters
    #'''
    #loc = args[1] #chrs loc
    #target = args[2]
    #paired = int(args[3])
    #genome_file = args[4]
    #'''

    #loc = '/data1/shunfu1/ref_shannon_modi/data/snyder_0922a/chrs/'
    #target = 'chr1_gl000191_random' #small data
    #target = 'chr15'
    #paired = 1
    #genome_file = '/data1/shunfu1/ref_shannon_modi/data/WingWongTest_K24_135M/genome/%s.fa'%target

    #subdir = out_dir #os.path.join(loc, target)
    #sam_file = os.path.join(subdir, 'hits.sam')

    output_dir = subdir + '/intermediate/'
    #pdb.set_trace()

    clock = Clock()

    #1. .sam --> regions_dic, genome2region_dic
    regions_dic, genome2region_dic = gen_regions_genome_mapping(sam_file, target)
    print('[%s] gen_regions_genome_mapping finished'%clock.asctime())
    #pdb.set_trace()

    #2. regions_dic, genome2region_dic, .sam --> region_pairs, known_paths
    region_pairs, known_paths = gen_edges_paths(regions_dic, genome2region_dic, sam_file, target, paired)
    print('[%s] gen_edges_paths finished'%clock.asctime())
    #pdb.set_trace()

    #3. regions --> seg
    region2seg, seg2region = gen_regions_segments_mapping(regions_dic, \
                                                          genome2region_dic, \
                                                          region_pairs.keys())
    print('[%s] gen_regions_segments_mapping finished'%clock.asctime())

    seg_res = check_segmentation(region2seg, seg2region)

    run_cmd('mkdir -p %s'%output_dir, shell=True)

    clock.time()
    if target=='':
        genome = from_fasta(genome_file)
        #del genome['_']
        genome = genome.items()[0][1]
    else:
        genome = from_fasta(genome_file)[target]
    print('{:.2f} to load genome\n\tgenome_file={}'.format(clock.time(), genome_file))
    #pdb.set_trace()
    
    export_nodes_edges_paths(output_dir, regions_dic, genome, region_pairs, known_paths, region2seg, seg2region)
    print('[%s] export_nodes_edges_paths finished'%clock.asctime())

    return

def calc_boundaries(weights, splice_starts, splice_ends):
    #key: genome pos (0-based) val:=1 (inclusive)/ -1 (inclusive)
    boundaries = {}

    weights_sorted = weights.items()
    weights_sorted.sort(key=lambda x:x[0])

    left = weights_sorted[0][0]
    for i in xrange(1,len(weights_sorted)):
        pos1 = weights_sorted[i-1][0]
        pos2 = weights_sorted[i][0]
        if pos2>pos1+1:
            right = pos1 #+1
            if right==left: #may lead to unpaired boundary
                print('island boundary at %d'%left)
                left = pos2
            else:
                boundaries[left]=1
                boundaries[right]=-1
                left = pos2

    right = weights_sorted[i][0] #+1
    if right != left:
        boundaries[left]=1
        boundaries[right]=-1
    #pdb.set_trace()

    '''
    debug unpaired boundary:
    use 0807b chr7 sam
    (Pdb) wt = weights.items()
    (Pdb) wt.sort(key=lambda x:x[0])
    (Pdb) islands = [wt[i] for i in xrange(1,len(wt)-1) if wt[i-1][0]+1!=wt[i][0] and wt[i][0]+1!=wt[i+1][0]]
    (Pdb) islands
    [(26766494, 1.0)]
    '''
    
    #splice starts
    for pos, _ in splice_starts.items():
        c1 = pos in weights
        c2 = pos+1 in weights
        c3 = pos not in boundaries
        c4 = pos+1 not in boundaries

        if c1 and c2 and c3 and c4:
            boundaries[pos]=-1
            boundaries[pos+1]=1
        else:
            st = ('unexpected at calc boundaries:\n')
            st += 'splice starts at pos=%d\n'%pos
            st += 'pos in weights: %s\n'%(c1)
            st += 'pos+1 in weights: %s\n'%(c2)
            st += 'pos not in boundaries: %s\n'%(c3)
            st += 'pos+1 not in boundaries: %s\n\n'%(c4)
            #logger.write(st)

    #pdb.set_trace()

    #splice ends
    for pos, _ in splice_ends.items():
        c1 = pos-1 in weights
        c2 = pos in weights
        c3 = pos-1 not in boundaries
        c4 = pos not in boundaries
        if c1 and c2 and c3 and c4:
           boundaries[pos-1]=-1
           boundaries[pos]=+1
        else:
            st = ('unexpected at calc boundaries\n')            
            st += 'splice ends at pos=%d\n'%pos
            st += 'pos-1 in weights: %s\n'%(c1)
            st += 'pos in weights: %s\n'%(c2)
            st += 'pos-1 not in boundaries: %s\n'%(c3)
            st += 'pos not in boundaries: %s\n\n'%(c4)
            #logger.write(st)

    #pdb.set_trace()

    boundaries_sorted = boundaries.items()
    boundaries_sorted.sort(key=lambda x:x[0])

    #pdb.set_trace()
    return boundaries_sorted #, weights_sorted

def gen_regions_modi(weights, boundaries, splice_starts, splice_ends): # dic of key - rid, val - [stt, stp, cov, # sp in, # sp out]

    #pdb.set_trace()

    regions_dic_0 = {}

    n_region = len(boundaries)/2

    for i in range(n_region):
        r_stt = boundaries[i*2][0]
        r_stp = boundaries[i*2+1][0]

        if not (boundaries[i*2][1]==1 and boundaries[i*2+1][1]==-1):
            print('gen regions modi - unpaired')
            pdb.set_trace()

        cov = 0
        for j in xrange(r_stt, r_stp+1):
            cov += weights[j]
        cov = float(cov)/(r_stp-r_stt+1)

        num_sp_in = 0
        if r_stt in splice_ends:
            num_sp_in = splice_ends[r_stt]

        num_sp_out = 0
        if r_stp in splice_starts:
            num_sp_out = splice_starts[r_stp]

        regions_dic_0[i] = [r_stt, r_stp, cov, num_sp_in, num_sp_out]

    #pdb.set_trace()
    return regions_dic_0 

def merge_regions_modi2(regions_dic_0):

    #pdb.set_trace()
    # dic of key - rid, val - [stt, stp, cov, # sp in, # sp out]

    regions_dic = {}
    msgs = []
    ind = 0 #new rid

    cnt_interpolated = 0

    stat_no_merge_has_splice = 0
    stat_no_merge_no_splice_large_gap = 0
    stat_merge_small_hole = 0
    stat_merge_close_hole_high_cov = 0
    stat_no_merge_misc = 0


    sorted_regions = regions_dic_0.items()
    sorted_regions.sort(key=lambda x:x[0]) #[rid0, [stt, stp, cov, # sp in, # sp out]]

    r1 = sorted_regions[0][1]
    for i in xrange(1, len(sorted_regions)):
        r2 = sorted_regions[i][1]

        if r1[4]>0 or r2[3]>0: #r1 has splice out or r2 has splice in, no merge
            regions_dic[ind]=copy.copy(r1[0:3])
            ind += 1
            r1 = r2

            stat_no_merge_has_splice += 1

        else: #r1 has no splice out and r2 has no splice in, consider merge (r1 and r2 has a gap; otherwise if they're connected, either r1 has a splice out or r2 has a splice in)
            if r1[1]+1==r2[0]:
                print('r1 and r2 connected with no splice out at r1 and no splice in at r2 -- unexpected')
                pdb.set_trace()

            elif r2[0]-r1[1]-1>=MAX_HOLE: #regions seperated far away, no merge            
               regions_dic[ind]=copy.copy(r1[0:3])
               ind+=1
               r1 = r2

               stat_no_merge_no_splice_large_gap += 1

            elif r2[0]-r1[1]-1<=ALLOWED_HOLE: #merge type I
                cnt_interpolated += (r2[0]-r1[1]-1)
                #print('r1=%s, r2=%s merged (cnt_interpolated=%d)'%(r1,r2, cnt_interpolated))
                msgs.append('r1=%s, r2(i=%d)=%s cnt_interpolated=%d (I)'%(r1,i,r2, cnt_interpolated))
                r_stt = r1[0]
                r_stp = r2[1]
                w1 = r1[2]*(r1[1]-r1[0]+1)
                w2 = r2[2]*(r2[1]-r2[0]+1)
                cov = float(w1+w2)/(r_stp-r_stt+1)
                num_sp_in = r1[3] #r1[4]==0
                num_sp_out = r2[4] #r2[3]==0
                r1 = [r_stt, r_stp, cov, num_sp_in, num_sp_out]

                stat_merge_small_hole += 1

            elif r2[0]-r1[1]-1<100 and r2[2]+r1[2]>=14: #merge type II
                #pdb.set_trace()
                cnt_interpolated += (r2[0]-r1[1]-1)
                #print('r1=%s, r2=%s merged (cnt_interpolated=%d)'%(r1,r2, cnt_interpolated))
                msgs.append('r1=%s, r2(i=%d)=%s cnt_interpolated=%d (II)'%(r1,i,r2, cnt_interpolated))
                r_stt = r1[0]
                r_stp = r2[1]
                w1 = r1[2]*(r1[1]-r1[0]+1)
                w2 = r2[2]*(r2[1]-r2[0]+1)
                cov = float(w1+w2)/(r_stp-r_stt+1)
                num_sp_in = r1[3] #r1[4]==0
                num_sp_out = r2[4] #r2[3]==0
                r1 = [r_stt, r_stp, cov, num_sp_in, num_sp_out]

                stat_merge_close_hole_high_cov += 1

            else: #no merge
                regions_dic[ind]=copy.copy(r1[0:3])
                ind+=1
                r1 = r2

                stat_no_merge_misc += 1

    regions_dic[ind]=copy.copy(r1[0:3])

    st = '\nsummary of merge_regions_modi2\n'
    st += 'stat_no_merge_has_splice=%d\n'%stat_no_merge_has_splice
    st += 'stat_no_merge_no_splice_large_gap=%d\n'%stat_no_merge_no_splice_large_gap
    st += 'stat_merge_small_hole (I) =%d\n'%stat_merge_small_hole
    st += 'stat_merge_close_hole_high_cov (II) =%d\n'%stat_merge_close_hole_high_cov
    st += 'stat_no_merge_misc=%d\n'%stat_no_merge_misc
    msgs.append(st)

    #pdb.set_trace()
    return regions_dic, msgs

def gen_genome2region(regions_dic):

    #pdb.set_trace()

    genome2region_dic = {} #key: genome pos, val: region index

    for region_idx, region in regions_dic.items():
        r_stt = region[0]
        r_stp = region[1]
        for i in xrange(r_stt, r_stp+1):
            genome2region_dic[i] = region_idx

    return genome2region_dic

def get_regions_from_read(genome2region_dic, alignment):
    regions = []
    for s, l in alignment.blocks:
        left = s
        right = s+l-1
        if left in genome2region_dic and right in genome2region_dic: #in case an island is not in genome2region_dic
            id1 = genome2region_dic[left]
            id2 = genome2region_dic[right]
            for i in xrange(id1, id2+1):
                regions.append(i)
    return regions

def gen_edges_paths_sam_SE(sam_file, target, genome2region_dic):

    #pdb.set_trace()

    region_pairs = {}
    known_paths = {} #set()

    num_lines = sum([1 for l in open(sam_file)]); print('%d lines in %s'%(num_lines, sam_file))
    i = 0; j=0; T=num_lines/100;

    with open(sam_file, 'r') as f:
    
        for line in f:
            i+=1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of sam processed'%j); sys.stdout.flush()

            if line[0]=='@': continue

            fields = line.split()
            name, flags, g_name, start, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]

            if target != '' and g_name != target: continue

            rev = bool((flags >> 4) & 1)

            if ((flags >> 2) & 1): continue #unmapped

            [valid_cigar, _, modi_blocks, splice_sites] = cigar2blocks5(start, cigar)
            if valid_cigar==0: continue #not M,N,S,I,D

            a = Alignment(g_name, modi_blocks, rev, '', 1.0, splice_sites)

            update_one_alignment(genome2region_dic, a, region_pairs, known_paths)

    print('')

    #known_paths = filter_known_paths(known_paths) #dic-->set
    return region_pairs, known_paths

def update_one_alignment(genome2region_dic, a, region_pairs, known_paths):

    regions = get_regions_from_read(genome2region_dic, a)
    regions = list(set(regions))
    regions.sort(key=lambda x:x)
    regions_tup = tuple(regions)

    if len(regions_tup)==2:
        if regions_tup in region_pairs:
            region_pairs[regions_tup]+=a.weight
        else:
            region_pairs[regions_tup]=a.weight

    elif len(regions_tup)>=3:
        
        pair = (regions_tup[0], regions_tup[1]) #e.g. (0,1)            
        if pair in region_pairs:
            region_pairs[pair]+=a.weight
        else:
            region_pairs[pair]=a.weight

        for i in xrange(2, len(regions_tup)):
            pair = (regions_tup[i-1], regions_tup[i]) #e.g. (1,2)            
            if pair in region_pairs:
                region_pairs[pair]+=a.weight
            else:
                region_pairs[pair]=a.weight

            path = (regions_tup[i-2], regions_tup[i-1], regions_tup[i]) #e.g. (0,1,2)
            #known_paths.add(path)
            add_known_paths(path, known_paths)
    return regions

def add_known_paths(path, known_paths):
    if path in known_paths:
        known_paths[path]+=1
    else:
        known_paths[path]=1

def filter_known_paths(known_paths): #dic-->set

    paths = known_paths.items()
    paths.sort(key=lambda x:-x[1])
    res = [p[0] for p in paths if p[1]>=10] #return known paths if confidence/cnt >=10
    print('%d/%d known paths kept'%(len(res), len(paths)))

    return res

def update_one_pair(genome2region_dic, rg, region_pairs, known_paths):

    #pdb.set_trace()
    #dmp_path = open('dmp/dmp_path.txt', 'a')
    #st = rg[0][0]+'\n'

    a = rg[0][1]
    a2 = rg[1][1]

    #from filter_mates
    if a.start() > a2.start(): a, a2 = a2, a
    if a.reversed or not a2.reversed: return #very few cases, assume it's in same region

    #from get_bridged_regions_PE_modi
    regions_1 = update_one_alignment(genome2region_dic, a, region_pairs, known_paths)    
    regions_2 = update_one_alignment(genome2region_dic, a2, region_pairs, known_paths)

    #if len(regions_1)>=3:
    #    st += '# blocks1='+str(a.blocks)+' regions1='+str(regions_1)+'\n'
    #if len(regions_2)>=3:
    #    st += '# blocks2='+str(a2.blocks)+' regions2='+str(regions_2)+'\n'

    #merge regions_1 and regions_2
    if regions_1[-1]+1==regions_2[0]:
        r_tup = (regions_1[-1], regions_2[0])
        if r_tup in region_pairs:
            region_pairs[r_tup]+=float(a.weight+a2.weight)/2
        else:
            region_pairs[r_tup]=float(a.weight+a2.weight)/2

        if len(regions_1)>=2:
            #known_paths.add((regions_1[-2], regions_1[-1], regions_2[0]))
            add_known_paths((regions_1[-2], regions_1[-1], regions_2[0]), known_paths)
            #st += '# interpolated regions_1 and 2 (1)\n'

        if len(regions_2)>=2:
            #known_paths.add((regions_1[-1], regions_2[0], regions_2[1]))
            add_known_paths((regions_1[-1], regions_2[0], regions_2[1]), known_paths)
            #st += '# interpolated regions_1 and 2 (2)\n'

    #st += '\n'

    #if '#' in st: dmp_path.write(st)
    #dmp_path.close()

    return


#@profile(stream=memory_fp)
def gen_edges_paths_sam_PE(sam_file, target, genome2region_dic):

    region_pairs = {}
    known_paths = {} #set()

    num_lines = sum([1 for l in open(sam_file)]); print('%d lines in %s'%(num_lines, sam_file))
    i = 0; j=0; T=num_lines/100;

    rg = [] #to store pair of reads
    with open(sam_file, 'r') as f:

        for line in f:
            i+=1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of sam processed'%j); sys.stdout.flush()

            if line[0]=='@': continue

            fields = line.split()
            name, flags, g_name, start, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]

            if target != '' and g_name != target: continue

            rev = bool((flags >> 4) & 1)

            if ((flags >> 2) & 1): continue #unmapped

            [valid_cigar, _, modi_blocks, splice_sites] = cigar2blocks5(start, cigar)
            if valid_cigar==0: continue #not M,N,S,I,D

            a = Alignment(g_name, modi_blocks, rev, '', 1.0, splice_sites)

            if rg == []:
                rg = [(name, a)] # list of (name, alignment)
            else:
                old_name = rg[-1][0]
                if old_name == name:
                    rg.append((name, a)) # a pair of reads prepared
                    update_one_pair(genome2region_dic, rg, region_pairs, known_paths)
                    rg = []
                    #pdb.set_trace()
                else:
                    #print('unexpected here')
                    #pdb.set_trace()
                    rg = [(name, a)]
    print('')

    #known_paths = filter_known_paths(known_paths) #dic-->set
    return region_pairs, known_paths

def gen_weights_splice(sam_file, target=''):

    weights, splice_starts, splice_ends = Counter(), Counter(), Counter()

    num_lines = sum([1 for l in open(sam_file)]); print('%d lines in %s'%(num_lines, sam_file))

    i = 0; j=0; T=num_lines/100;

    with open(sam_file, 'r') as f:

        for line in f:
            i+=1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of sam processed'%j); sys.stdout.flush()

            if line[0]=='@': continue

            fields = line.split()
            name, flags, g_name, start, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]

            if target != '' and g_name != target: continue

            if ((flags >> 2) & 1): continue #unmapped

            [valid_cigar, _, modi_blocks, splice_sites] = cigar2blocks5(start, cigar)
            if valid_cigar==0: continue #not M,N,S,I,D

            for sp_l, sp_r in splice_sites:
                splice_starts[sp_l] += 1
                splice_ends[sp_r] += 1

            for start, length in modi_blocks:
                for wp in xrange(start, start + length):
                    weights[wp] += 1

    print('')
    return weights, splice_starts, splice_ends


#@profile(stream=memory_fp)
def gen_regions_genome_mapping(sam_file, target=''):

    clock = Clock()

    weights, splice_starts, splice_ends = gen_weights_splice(sam_file, target)
    print('%.2f to gen_weights_splice'%clock.time())

    boundaries = calc_boundaries(weights, splice_starts, splice_ends)
    print('%.2f to calculate exon boundaries'%clock.time())

    regions_dic_0 = gen_regions_modi(weights, boundaries, splice_starts, splice_ends) # dic of key - rid, val - [stt, stp, cov, # sp in, # sp out]
    print('%.2f to generate regions (modi)'%clock.time())

    del weights; del splice_starts; del splice_ends

    regions_dic, msgs = merge_regions_modi2(regions_dic_0)
    print('%.2f to merge regions'%clock.time())

    genome2region_dic = gen_genome2region(regions_dic)
    print('%.2f to gen genome2region'%clock.time())

    return regions_dic, genome2region_dic

def gen_edges_paths(regions_dic, genome2region_dic, sam_file, target, paired):

    #region_pairs = {} #key - (r1, r2) if r1 & r2 are bridged
    #known_paths = set() #set of tuples, 1 tuple contains 3 elements

    if paired==0:
        return gen_edges_paths_sam_SE(sam_file, target, genome2region_dic)
    else:
        return gen_edges_paths_sam_PE(sam_file, target, genome2region_dic)

def gen_regions_segments_mapping(regions_dic, genome2region_dic, region_pairs):

    region2seg = {}
    seg2region = {}

    #init
    rid_list = []
    for r in regions_dic.items():
        rid = r[0]
        region2seg[rid] = rid
        seg2region[rid] = [rid]
        rid_list.append(rid)
    next_sid = max(rid_list)+1
    #pdb.set_trace()

    #merge connected regions with same segment id
    for rid1, rid2 in region_pairs:
        sid1 = region2seg[rid1]
        sid2 = region2seg[rid2]
        
        if sid1==sid2:
            continue

        rids_sid1 = seg2region[sid1]
        rids_sid2 = seg2region[sid2]

        rids_merge = list(set(rids_sid1 + rids_sid2))
        rids_merge.sort(key=lambda x:x)

        for rid in rids_merge:
            region2seg[rid] = next_sid

        seg2region[next_sid] = rids_merge

        del seg2region[sid1]
        del seg2region[sid2]

        next_sid += 1

    #pdb.set_trace()

    return region2seg, seg2region

def check_segmentation(region2seg, seg2region):
    res = True #True: segmentation looks good

    rids1 = set()
    sids1 = set()
    for rid, sid in region2seg.items():
        rids1.add(rid)
        sids1.add(sid)

    rids2 = set()
    sids2 = set()
    for sid, rids in seg2region.items():
        sids2.add(sid)
        for rid in rids:
            rids2.add(rid)

    if len(rids1)!=len(rids2):
        print('unexpected segmentation: len(rids1)==%d, len(rids2)==%d'%(len(rids1), len(rids2)))
        res = False

    if len(sids1)!=len(sids2):
        print('unexpected segmentation: len(sids1)==%d, len(sids2)==%d'%(len(sids1), len(sids2)))
        res = False

    print('check segmentation: %s'%res)

    return res

def str_region2(rid, stt_end_cov, genome):
    stt = stt_end_cov[0]
    stp = stt_end_cov[1]
    cov = stt_end_cov[2]

    st = ''
    st += '%d\t'%rid
    st += '%s\t'%genome[stt:stp+1]
    st += '%f\t'%cov
    st += '%d\t'%(stp-stt+1)
    #gtf purpose
    st += '%d\t'%stt #0-based, inclusive
    st += '%d\n'%stp #0-based, inclusive

    return st

#@profile(stream=fp)
def export_1seg(output_dir, sorted_rids, regions_dic, genome, region_pairs, paths_by_start, component):

    num_nodes = 0
    num_edges = 0
    num_paths = 0

    with open(output_dir+'/nodes'+str(component)+'.txt', 'w') as nodefile, \
         open(output_dir+'/edges'+str(component)+'.txt', 'w') as edgefile, \
         open(output_dir+'/paths'+str(component)+'.txt', 'w') as pathfile, \
         open(output_dir+'/node_hash'+str(component)+'.txt', 'w') as node_hash_file:

        node_hash_file.write('org_rid\thash_id\n')
        node_hash_dic = {} #key rid, val hash of rid
        r_hash = 0
        for rid in sorted_rids:
            node_hash_dic[rid] = r_hash
            node_hash_file.write('%d\t%d\n'%(rid, r_hash))
            r_hash += 1
            num_nodes += 1

        #pdb.set_trace()

        #nodefile.write('ID\tBases\tCopycount\tNormalization\n')
        nodefile.write('ID\tBases\tCopycount\tNormalization\tGenome_start\tGenome_stop\n')
        for rid in sorted_rids:
            #st = str_region(node_hash_dic[rid], regions_dic[rid], genome)
            st = str_region2(node_hash_dic[rid], regions_dic[rid], genome)
            nodefile.write(st)

        #pdb.set_trace()

        pathfile.write('ID1\tID2\tEtc.\n')
        for rid in sorted_rids:
            if rid not in paths_by_start:
                continue

            paths = paths_by_start[rid]
            for path in paths:
                path = [str(node_hash_dic[n]) for n in path]
                pathfile.write('\t'.join(path) + '\n')
                num_paths += 1

        #pdb.set_trace()

        edgefile.write('InID\tOutID\tWeight\tCopycount\tNormalization\n')
        for edge, cnt in region_pairs:
            st = ''
            st += '%d\t'%node_hash_dic[edge[0]]
            st += '%d\t'%node_hash_dic[edge[1]]
            st += '0\t' #weight: overlap
            st += '%d\t'%cnt
            st += '1\n' #normalization
            edgefile.write(st)
            num_edges += 1

    return [num_nodes, num_edges, num_paths]

def export_nodes_edges_paths(output_dir, regions_dic, genome, region_pairs, known_paths, region2seg, seg2region):
    '''
    regions_dic: key rid val [stt, stp, cov] 0-based, inclusive
    genome: str
    region_pairs: key (rid1, rid2) val cnt of reads passing the boundary of rid1 and rid2
    known_paths: set of tuple of rids which are passed by certain read
    region2seg: key rid val sid
    seg2region: key sid val list of rids (e.g. connected components)
    '''

    if output_dir is None:
        return

    clock = Clock()

    clock.time()

    paths_by_start = {} #key: rid, val: list of paths, each of which starts at rid
    for path in known_paths:
        if path[0] not in paths_by_start:
            paths_by_start[path[0]]=[]
        paths_by_start[path[0]].append(path)

    print('{:.2f} to get paths_by_start'.format(clock.time()))

    clock.time()

    with open(output_dir+'/single_nodes.txt', 'w') as single_file:
        #single_file.write('org_rid\tBases\tCopycount\tNormalization\n')
        single_file.write('org_rid\tBases\tCopycount\tNormalization\tGenome_start\tGenome_stop\n')
        for sid, rids in seg2region.items():
            if len(rids)==1:
                #st = str_region(rids[0], regions_dic[rids[0]], genome)
                st = str_region2(rids[0], regions_dic[rids[0]], genome)
                single_file.write(st)

    print('{:.2f} to export single nodes'.format(clock.time()))

    clock.time()

    sid2component = {}
    component = 0
    with open(output_dir+'/stats.txt', 'w') as stat_file:
        stat_file.write('org_sid\tcomp_id\tnum_nodes\n') #\tnum_edges\tnum_paths\n')
        for sid, rids in seg2region.items():
            num_nodes = len(rids)
            if num_nodes==1: continue
            stat_file.write('%d\t%d\t%d\n'%(sid, component, num_nodes))
            sid2component[sid]=component
            component += 1

    print('{:.2f} to export stats'.format(clock.time()))

    clock.time()
    seg2edges = {} # key sid val list of [edge, cnt] belonging to same seg

    for edge, cnt in region_pairs.items():

        rid0 = edge[0]
        sid0 = region2seg[rid0]
        rid1 = edge[1]
        sid1 = region2seg[rid1]
        
        if sid0 != sid1:
            print('unexpected edge: edge[0]=%d, edge[1]=%d'%(edge[0], edge[1]))
            pdb.set_trace()

        if sid0 in seg2edges:
            seg2edges[sid0].append([edge, cnt])
        else:
            seg2edges[sid0] = [[edge, cnt]]

    print('{:.2f} to build seg2edges'.format(clock.time()))
    #pdb.set_trace()

    clock.time()

    seg_list = seg2region.items()
    seg_list.sort(key=lambda x:x[0]) #sort by sid

    T = len(seg_list)/100
    i = 0
    j = 0
    #pdb.set_trace()

    for sid, rids in seg_list:

        i += 1
        if T==0:
            sys.stdout.write('\r'); sys.stdout.write('%.2f %% processed (export)'%(float(i)*100/len(seg_list))); sys.stdout.flush()
        else:
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (export)'%j); sys.stdout.flush()

        if len(rids)==1: continue

        [num_nodes, num_edges, num_paths] = export_1seg(output_dir, rids, regions_dic, genome, seg2edges[sid], paths_by_start, sid2component[sid])

    print('\n{:.2f} to exports >2 nodes'.format(clock.time()))

    return

if __name__ == '__main__':

    gen_splice_graph2()
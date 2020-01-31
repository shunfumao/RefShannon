#a wrapper to run sparse flow
#python sf.py chrs_loc target scheduler_file_index python_path shannon_dir

import sys
import pdb
import os
from RefShannon.util import run_cmd

ROOT = os.path.dirname(__file__)
pdb.set_trace()

def load_stat_file(stat_file):

    #'org_sid\tcomp_id\tnum_nodes\n

    #pdb.set_trace()

    stats = []
    with open(stat_file, 'r') as f:

        i=-1
        for line in f:
            i+=1
            if i>0: #skip first line
                tokens = line.split()
                org_sid = int(tokens[0])
                comp_id = int(tokens[1])
                num_nodes = int(tokens[2])
                stats.append([org_sid, comp_id, num_nodes])
                #num_edges = int(tokens[3])
                #num_paths = int(tokens[4])
                #stats.append([org_sid, comp_id, num_nodes, num_edges, num_paths])

    #pdb.set_trace()

    return stats

def build_scheduler_files(stats, N_jobs, loc):

    #stats 'org_sid\tcomp_id\tnum_nodes\n #\tnum_edges\tnum_paths\n'
    #pdb.set_trace()

    scheduler_indice = []
    scheduler_files = []
    for j in range(N_jobs):
        scheduler_indice.append(j)
        f = open(loc+'/scheduler%d.txt'%j, 'w')
        scheduler_files.append(f)

    stats.sort(key=lambda x:x[2]) #sort based on num of nodes

    next_scheduler = 0
    for itm in stats:
        comp_id = itm[1]
        scheduler_files[next_scheduler].write('%d\n'%comp_id)
        next_scheduler = (next_scheduler + 1)%N_jobs

    for j in range(N_jobs):
        scheduler_files[j].close()

    #pdb.set_trace()

    return scheduler_indice

#def run_sf(chrs_loc, target, scheduler_file_index):
def run_sf(sg_dir, res_dir, tr_name_tag, target_str, scheduler_file_index, F_val, outputFasta_str):

    component_list = []
    if scheduler_file_index==-1: #serial
        ncomp = 0
        while os.path.isfile(sg_dir + '/nodes'+str(ncomp)+'.txt'):
            component_list.append(ncomp)
            ncomp += 1
    else: #parallel
        with open(sg_dir + '/scheduler'+str(scheduler_file_index)+'.txt', 'r') as f:
            for line in f:
                component_list.append(int(line.strip()))

    for i in range(len(component_list)):
        ncomp = component_list[i]
        sys.stdout.write('\r%d/%d of components processed %s'%((i+1), len(component_list), target_str)); sys.stdout.flush()
        run_cmd('python %s/algorithm_SF.py '%ROOT + str(ncomp) + ' -I '+ sg_dir + ' -O ' + res_dir + \
            ' -tr_name_tag ' + tr_name_tag + ' ' + target_str + ' -F %f '%F_val + outputFasta_str)

    print('')

    return

'''
usage:

python sf.py -I sg_dir -O res_dir -tr_name_tag tr_name_tag [-target target] [-F F_val] -scheduler_index sid [--outputFasta]

# F_val: at local sparse flow decomposition, filter flows if fij/max(fij)<F_val to trade-off sens/fp
# outputFasta: sparse flow procedure will output both gtf and fasta files

'''

if __name__ == '__main__':

    #pdb.set_trace()

    #chrs_loc = sys.argv[1]
    #target = sys.argv[2]
    #scheduler_file_index = int(sys.argv[3])

    #run_sf(chrs_loc, target, scheduler_file_index)

    args = sys.argv

    sg_dir = args[args.index('-I')+1] #e.g. ../intermediate/
    res_dir = args[args.index('-O')+1] #e.g. ../algo_output/
    tr_name_tag = args[args.index('-tr_name_tag')+1] #e.g. refShannon_chr15

    if '-target' in args:
        target = args[args.index('-target')+1]
        target_str = '-target '+target
    else:
        target = ''
        target_str = ''

    if '-F' in args:
        F_val = float(args[args.index('-F')+1])
    else:
        F_val = 0.0

    if '--outputFasta' in args:
        outputFasta_str = '--outputFasta'
    else:
        outputFasta_str = ''

    scheduler_file_index = int(args[args.index('-scheduler_index')+1])

    run_sf(sg_dir, res_dir, tr_name_tag, target_str, scheduler_file_index, F_val, outputFasta_str)


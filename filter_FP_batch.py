import sys, pdb
from filter_FP.filter_FP import *

def process_oneThreshold(args):

    rec_fasta = args[args.index('-i')+1]
    fp_rec_file = args[args.index('-o')+1]
    depth_file = args[args.index('-d')+1]
    log_file = args[args.index('-l')+1]
    threshold = float(args[args.index('-t')+1])

    write_filtered_tr(depth_file, rec_fasta, fp_rec_file, log_file, threshold)
    print('%s and %s written (t=%f)'%(fp_rec_file, log_file, threshold))

    return

def process_multiThresholds(args):

    out_dir = args[args.index('-O')+1]
    run_cmd('mkdir -p %s'%out_dir)

    thr_str = args[args.index('-t')+1]
    thr_list = [float(itm) for itm in thr_str.split(',') if itm != '']

    fp_rec_files = []

    #pdb.set_trace()

    for thr in thr_list:

        args2 = []
        
        args2 += args[args.index('-i'):args.index('-i')+2]

        fp_rec_file = '%s/rec_t_%.2f.fasta'%(out_dir, thr)
        fp_rec_files.append((fp_rec_file, 'rec_t_%.2f'%thr))

        log_file = '%s/rec_t_%.2f.log'%(out_dir, thr)

        args2.append('-o'); args2.append(fp_rec_file)

        args2 += args[args.index('-d'):args.index('-d')+2]

        args2.append('-l'); args2.append(log_file)

        args2.append('-t'); args2.append(str(thr))

        process_oneThreshold(args2)

    if '--eval' in args:

        #pdb.set_trace()

        ref_file = args[args.index('-r')+1]

        out_dir2 = out_dir + '/eval/'
        run_cmd('mkdir -p %s'%out_dir2)      

        for fp_rec_file, name_tag in fp_rec_files:

            cmd = 'python refShannon.py --eval -i %s -r %s -O %s -n %s'%(fp_rec_file, ref_file, out_dir2, name_tag)
            run_cmd(cmd)
    return

'''
usage:

python filter_FP_batch.py --oneThreshold -i rec_fasta -o fp_rec_file -d depth_file -l log_fie -t threshold

###if a transcript is re-covered by over 90%, output it

python filter_FP_batch.py --multiThresholds -i rec_fasta -O out_dir -d depth_file -t [threshold1,[threshold2,...]] [--eval -r ref_file]

###filtered fp files are rec_t_thr.fasta (and .log) stored at out_dir,
###eval files rec_t_thr_log.txt (and _per.txt) stored at out_dir/eval/ 

'''
if __name__ == '__main__':

    args = sys.argv

    if '--oneThreshold' in args:
        process_oneThreshold(args)
    elif '--multiThresholds' in args:
        process_multiThresholds(args)
    else:
        print('unexpected mode')
        pdb.set_trace() 


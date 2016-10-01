import sys, pdb, os
from util import *
from sf import load_stat_file, build_scheduler_files
import tester

'''
reads --> sam

python refShannon.py --alignment -o sam_file -g genome_file -r read_file [-N num_Thread]
python refShannon.py --alignment -o sam_file -g genome_file -r1 read_file1 -r2 read_file2 [-N num_Thread]

sam --> chromosome-based sam

python refShannon.py --split -O chrs_dir -i sam_file

sam --> splice graph

python refShannon.py --sg -i sam_file -g genome_file [-O out_dir] [-paired] [-target chr]
python refShannon.py --sg -I chrs_dir -g multi_genome_file [-O out_dir] [-chrs chr_a[,chr_b,...]] [-paired] [-N nJobs]
python refShannon.py --sg -I chrs_dir -G genome_dir [-O out_dir] [-chrs chr_a[,chr_b,...]] [-paired] [-N nJobs]

splice graph --> transcripts

python refShannon.py --sf -I sg_dir [-O res_dir] [-N nJobs] [-target chr]
python refShannon.py --sf -Is chrs_dir [-O res_dir] [-chrs chr_a[,chr_b,...]] [-N nJobs]

transcripts --> performance (per.txt, log.txt)

python refShannon.py --eval -i fa_file -r ref_file [-O out_dir] [-n name_tag]
python refShannon.py --eval -I chrs_dir [-chrs chr_a[,chr_b,...]] -r ref_file [-O out_dir] [-n name_tag]

sam --> transcripts

python refShannon.py --batch -i sam_file -g genome_file [-O out_dir] [-paired] [-target chr] [-N nJobs] [-clear]
python refShannon.py --batch -I chrs_dir -g multi_genome_file [-O out_dir] [-chrs chr_a[,chr_b,...]] [-paired] [-N nJobs] [-clear]
python refShannon.py --batch -I chrs_dir -G genome_dir [-O out_dir] [-chrs chr_a[,chr_b,...]] [-paired] [-N nJobs] [-clear]

'''

def do_alignment(args):

    cmd = 'python aligner_star.py %s'%(' '.join(args[2:]))
    run_cmd(cmd)

    return

def do_split(args):

    cmd = 'python split.py %s'%(' '.join(args[2:]))
    run_cmd(cmd)

    return

def do_sg_i_g(args):

    cmd = 'python gen_sg2.py %s'%(' '.join(args)) #(' '.join(args[2:]))
    run_cmd(cmd)

    return

def do_sg_I_g(args):

    chrs_dir = args[args.index('-I')+1]
    multi_genome_file = args[args.index('-g')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir)

    if '-paired' in args:
        paired = '-paired'
    else:
        paired = ''

    if '-N' in args and len(target_list)>1:
        nJobs = int(args[args.index('-N')+1])
    else:
        nJobs = 1

    if nJobs==1:

        for target in target_list:

            sam_file = chrs_dir + '/' + target + '/hits.sam'
            target_out_dir = '%s/%s/'%(out_dir, target)
            cmd = 'python gen_sg2.py -i %s -g %s -O %s %s -target %s'%(sam_file, multi_genome_file, target_out_dir, paired, target)
            run_cmd(cmd)

    elif nJobs>1:

        pstring = ' '.join(target_list)
        cmd = 'time parallel --no-notice --jobs %d python gen_sg2.py -i %s/{}/hits.sam -g %s -O %s/{}/ %s -target {} ::: %s'%(nJobs, chrs_dir, multi_genome_file, out_dir, paired, pstring)
        run_cmd(cmd)

    return

def do_sg_I_G(args):

    chrs_dir = args[args.index('-I')+1]
    genome_dir = args[args.index('-G')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir)

    if '-paired' in args:
        paired = '-paired'
    else:
        paired = ''

    if '-N' in args and len(target_list)>1:
        nJobs = int(args[args.index('-N')+1])
    else:
        nJobs = 1

    if nJobs==1:

        for target in target_list:

            sam_file = chrs_dir + '/' + target + '/hits.sam'
            target_out_dir = '%s/%s/'%(out_dir, target)
            cmd = 'python gen_sg2.py -i %s -g %s/%s.fa -O %s %s -target %s'%(sam_file, genome_dir, target, target_out_dir, paired, target)
            run_cmd(cmd)

    elif nJobs>1:

        pstring = ' '.join(target_list)
        cmd = 'time parallel --no-notice --jobs %d python gen_sg2.py -i %s/{}/hits.sam -g %s/{}.fa -O %s/{}/ %s -target {} ::: %s'%(nJobs, chrs_dir, genome_dir, out_dir, paired, pstring)
        run_cmd(cmd)

    return

def do_sg(args):

    clock = Clock()
    stt_time = clock.asctime()

    if '-i' in args and '-g' in args:
        do_sg_i_g(args)
    elif '-I' in args and '-g' in args:
        do_sg_I_g(args)
    elif '-I' in args and '-G' in args:
        do_sg_I_G(args)
    else:
        print('unexpected args for refShannon.py --sg')

    print('[{}] to [{}] finish running splice graph generation'.format(stt_time, clock.asctime()))    

    return

def do_sf_I(args):

    sg_dir = args[args.index('-I')+1]

    if '-O' in args:
        res_dir = args[args.index('-O')+1]
    else:
        res_dir = parent_dir(sg_dir)+'/algo_output/'

    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-target' in args:
        target = args[args.index('-target')+1]
        target_str = ' -target '+target
        tr_name_tag = 'refShannon_%s'%target
        #pdb.set_trace()
    else:
        target = ''
        target_str = ''
        tr_name_tag = 'refShannon'

    stat_file = sg_dir + '/stats.txt'

    #single node
    run_cmd('python algorithm_SF.py -1 -I '+ sg_dir + ' -O ' + res_dir + \
            ' -tr_name_tag ' + tr_name_tag + target_str)

    if N_jobs == 1:        
        cmd = 'python sf.py -I %s -O %s -tr_name_tag %s %s -scheduler_index %d'%(sg_dir, res_dir, tr_name_tag, target_str, -1)
        run_cmd(cmd)

    elif N_jobs > 1:
        stats = load_stat_file(stat_file)
        scheduler_indice = build_scheduler_files(stats, N_jobs, sg_dir)
        pstring = ' '.join([str(i) for i in scheduler_indice])
        cmd = 'time parallel --no-notice --jobs %d python sf.py -I %s -O %s -tr_name_tag %s %s -scheduler_index {} ::: %s' \
              %(N_jobs, sg_dir, res_dir, tr_name_tag, target_str, pstring)
        run_cmd(cmd)
        
        cmd = 'rm %s/scheduler*.txt'%sg_dir
        run_cmd(cmd)

    #merge res
    wd = os.getcwd()
    os.chdir(res_dir) #in case cat arg too long
    reconstr_file = 'reconstructed.fasta'
    reconstr_gtf = 'reconstructed.gtf'
    os.system('cat reconstructed_comp_*.fasta > ' + reconstr_file) # '>>' to '>' for both Y on and Y off
    os.system('rm reconstructed_comp_*.fasta')
    os.system('cat reconstructed_comp_*.gtf > ' + reconstr_gtf) # '>>' to '>' for both Y on and Y off
    os.system('rm reconstructed_comp_*.gtf')
    os.chdir(wd)#change back cwd

    return

def do_sf_Is(args):

    chrs_dir = args[args.index('-Is')+1]

    if '-O' in args:
        res_dir = args[args.index('-O')+1]
    else:
        res_dir = chrs_dir # e.g. chrs/chr1/algo_output

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir)

    if '-N' in args:
        #pdb.set_trace()
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    for target in target_list:
        args_per_target =  '-I %s/%s/intermediate/ '%(chrs_dir, target)
        args_per_target += '-O %s/%s/algo_output/ '%(res_dir, target)
        args_per_target += '-N %d '%(N_jobs)
        args_per_target += '-target %s '%(target)
        args_per_target = args_per_target.split()
        #pdb.set_trace()
        do_sf_I(args_per_target)

    return

def do_sf(args):

    clock = Clock()
    stt_time = clock.asctime()

    if '-I' in args:
        do_sf_I(args)
    elif '-Is' in args:
        do_sf_Is(args)
    else:
        print('unexpected args for refShannon.py --sf')

    print('[{}] to [{}] finish running sparse flow decomposition'.format(stt_time, clock.asctime()))    

    return

def do_eval_i(args):

    fa_file = args[args.index('-i')+1] #reconstructed .fasta file

    if '-r' in args:
        ref_file = args[args.index('-r')+1]
    else:
        print('unexpected args for refShannon.py --eval -i -r')
        return

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = parent_dir(fa_file) #default: same as fasta folder

    if '-n' in args:
        name_tag = args[args.index('-n')+1]
    else:
        name_tag = 'reconstructed'

    per_file = out_dir + '/%s_per.txt'%name_tag
    log_file = out_dir + '/%s_log.txt'%name_tag

    run_cmd('python blat.py {} {} {}'.format(fa_file, ref_file, per_file))
    tester.analyzer_blat_noExp(per_file, log_file, '', 0)

    return

def do_eval_I(args):

    chrs_dir = args[args.index('-I')+1]

    if '-r' in args:
        ref_file = args[args.index('-r')+1]
    else:
        print('unexpected args for refShannon.py --eval -I -r')
        return

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:
        chrs_str = args[args.index('-chrs')+1]
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:
        target_list = os.listdir(chrs_dir)

    if '-n' in args:
        name_tag = args[args.index('-n')+1]
    else:
        name_tag = 'reconstructed'

    #merge
    target_file = out_dir + '/%s_all.fasta'%name_tag #e.g. reconstructed_all.fasta, stringtie_all.fasta
    temp_file = out_dir + '/%s_temp.fasta'%name_tag

    for target in target_list:        
        source_file = chrs_dir + '/' + target + '/algo_output/%s.fasta'%name_tag

        if os.path.exists(target_file)==False:
            run_cmd('touch %s'%target_file)

        if os.path.exists(source_file)==False:
            print('%s not found'%source_file)
            continue

        cmd = 'cat %s %s > %s'%(target_file, source_file, temp_file)
        run_cmd(cmd)

        cmd = 'mv %s %s'%(temp_file, target_file)
        run_cmd(cmd)

    args2 = '-i %s -r %s -O %s -n %s'%(target_file, ref_file, out_dir, name_tag)
    args2 = args2.split()
    do_eval_i(args2)

    return

def do_eval(args):

    clock = Clock()
    stt_time = clock.asctime()

    if '-i' in args:
        do_eval_i(args)
    elif '-I' in args:
        do_eval_I(args)
    else:
        print('unexpected args for refShannon.py --eval')

    print('[{}] to [{}] finish evaluation'.format(stt_time, clock.asctime()))    

    return

def do_batch_i_g(args):

    sam_file = args[args.index('-i')+1]
    genome_file = args[args.index('-g')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = parent_dir(sam_file)

    if '-paired' in args:
        paired_str = '-paired'
    else:
        paired_str = ''

    if '-target' in args:
        target = args[args.index('-target')+1]
        target_str = '-target %s'%target
    else:
        target = ''
        target_str = ''

    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-clear' in args:
        clear = True
    else:
        clear = False

    #gen sg
    sg_args = '-i %s -g %s -O %s %s %s'%(sam_file, genome_file, out_dir, paired_str, target_str)
    do_sg_i_g(sg_args.split())

    sg_dir = out_dir + '/intermediate/'
    res_dir = out_dir + '/algo_output/'
    sf_args = '-I %s -O %s -N %d %s'%(sg_dir, res_dir, N_jobs, target_str)
    do_sf_I(sf_args.split())

    if clear==True:
        run_cmd('rm -r %s'%sg_dir)

    return

def do_batch_I_g(args):

    chrs_dir = args[args.index('-I')+1]
    multi_genome_file = args[args.index('-g')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:        
        chrs_str = args[args.index('-chrs')+1]
        chrs_str_arg = '-chrs %s'%chrs_str
        target_list = [itm for itm in chrs_str.split(',') if itm != '']        
    else:        
        chrs_str = ''
        chrs_str_arg = ''
        target_list = os.listdir(chrs_dir)

    if '-paired' in args:
        paired_str = '-paired'
    else:
        paired_str = ''

    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-clear' in args:
        clear = True
    else:
        clear = False

    sg_args = '-I %s -g %s -O %s %s %s -N %d'%(chrs_dir, multi_genome_file, out_dir, chrs_str_arg, paired_str, N_jobs)
    do_sg_I_g(sg_args.split())

    sf_args = '-Is %s -O %s %s -N %d'%(out_dir, out_dir, chrs_str_arg, N_jobs)
    do_sf_Is(sf_args.split())

    if clear==True:
        for target in target_list:
            fld = out_dir + '/' + target + '/intermediate/'
            run_cmd('rm -r %s'%fld)

    return

def do_batch_I_G(args):

    chrs_dir = args[args.index('-I')+1]
    genome_dir = args[args.index('-G')+1]

    if '-O' in args:
        out_dir = args[args.index('-O')+1]
    else:
        out_dir = chrs_dir

    if '-chrs' in args:        
        chrs_str = args[args.index('-chrs')+1]
        chrs_str_arg = '-chrs %s'%chrs_str
        target_list = [itm for itm in chrs_str.split(',') if itm != '']
    else:        
        chrs_str = ''
        chrs_str_arg = ''
        target_list = os.listdir(chrs_dir)

    if '-paired' in args:
        paired_str = '-paired'
    else:
        paired_str = ''

    if '-N' in args:
        N_jobs = int(args[args.index('-N')+1])
    else:
        N_jobs = 1

    if '-clear' in args:
        clear = True
    else:
        clear = False

    sg_args = '-I %s -G %s -O %s %s %s -N %d'%(chrs_dir, genome_dir, out_dir, chrs_str_arg, paired_str, N_jobs)
    do_sg_I_G(sg_args.split())

    sf_args = '-Is %s -O %s %s -N %d'%(out_dir, out_dir, chrs_str_arg, N_jobs)
    do_sf_Is(sf_args.split())

    if clear==True:
        for target in target_list:
            fld = out_dir + '/' + target + '/intermediate/'
            run_cmd('rm -r %s'%fld)

    return

def do_batch(args):

    if '-i' in args and '-g' in args:
        do_batch_i_g(args)
    elif '-I' in args and '-g' in args:
        do_batch_I_g(args)
    elif '-I' in args and '-G' in args:
        do_batch_I_G(args)
    else:
        print('unexpected args for refShannon.py --batch')

    return

def main():

    args = sys.argv

    mode = args[1]

    if mode=='--alignment':
        do_alignment(args)

    elif mode=='--split':
        do_split(args)

    elif mode=='--sg':
        do_sg(args)

    elif mode=='--sf':
        do_sf(args)

    elif mode=='--eval':
        do_eval(args)

    elif mode=='--batch':
        do_batch(args)

    else:
        print('unknown mode: %s'%mode)

    return

if __name__ == '__main__':
    main()

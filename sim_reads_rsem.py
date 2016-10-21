import sys, pdb

'''
#### example rsem scripts:

(1) generate transcript reference

rsem-prepare-reference \
[--gtf /data1/shunfu1/ref_shannon_modi/data/snyder_0729a/reference/gencode.v15_PB_enhanced_HOP_Gm12878.ExonLinesOnly.gtf] \
--star \
-p 20 \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/genome/human/ \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref

(2) transcript abundance estimation

at [/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/], run:

rsem-calculate-expression \
--star \
--paired-end \
-p 20 \
/data1/lucille2/Full_Assembler/Snyder/PE_algo_input/reads_1.fastq \
/data1/lucille2/Full_Assembler/Snyder/PE_algo_input/ reads_2.fastq \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref \
exp

(3) learn statistics from (1)~(2) to generate simulated reads:

at [/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/], run:

rsem-simulate-reads \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref \
exp.stat/exp.model \
exp.isoforms.results \
theta0(first value of the third line of the file 'exp.stat/exp.theta' ==> 0.0429317239932184) \
N (the 4th number of the first line of the file 'sample_name.stat/sample_name.cnt' ==> 150M) \
output_name(path/to/sim_reads ==> /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads) \
--seed 0

example:
at [/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/], run:
rsem-simulate-reads /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/exp.stat/exp.model \
/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/exp.isoforms.results \
0.043 150000000 sim_reads --seed 0

#### Note:
(2) Time Used for EM.cpp : 3 h 04 m 55 s

for generated sim_reads, strandness is not specified, 
should be same as rsem-calculate-expression (e.g. --strandedness is none by default)

'''


def filterGTF(args):

    gtfFile = args[args.index('-i')+1]
    filtered_gtfFile = args[args.index('-o')+1]

    features_str = args[args.index('--keepFeatures')+1]
    features_list = [itm for itm in features_str.split(',') if itm != '']

    num_lines = sum([1 for l in open(gtfFile)]); print('%d lines in %s'%(num_lines, gtfFile))
    i = 0; j=0; T=num_lines/100;

    with open(gtfFile, 'r') as fin, \
         open(filtered_gtfFile, 'w') as fout:

        for line in fin:

            i+=1
            if i>T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% of gtf processed'%j); sys.stdout.flush()

            if line.split()[2] not in features_list:
                #pdb.set_trace()
                continue

            fout.write(line)

        print('')
        
    return

'''
Usage:

(1) python sim_reads_rsem.py --filterGTF -i path/to/gtf -o path/to/gtf_filtered --keepFeatures feature1[,feature2,...]

    originally used to keep only exon features in gtf file in WingWong reference transcript data, in order to run RSEM

'''
if __name__ == '__main__':

    args = sys.argv

    if '--filterGTF' in args:
        filterGTF(args)

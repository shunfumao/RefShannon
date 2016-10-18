'''
example rsem scripts:

(1) generate transcript reference

rsem-prepare-reference \
--gtf /data1/shunfu1/ref_shannon_modi/data/snyder_0729a/reference/gencode.v15_PB_enhanced_HOP_Gm12878.ExonLinesOnly.gtf \
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
theta0(first value of the third line of the file ‘exp.stat/exp.theta') \
N(the 4th number of the first line of the file 'sample_name.stat/sample_name.cnt') \
output_name(path/to/sim_reads) \
--seed 0
'''
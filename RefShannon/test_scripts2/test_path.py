"""
Modules
"""
path_test_aligner_hisat2 = {
  "wwSimChr15Hisat2": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/wwSimChr15/hisat2.sam",
    "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/ww_snyderRef/sim_reads.fa",
    # small file for test
    # "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/test_sgRefShannon/reads_1.fa",
    "-fasta",
    "-N", "20",
    "-sorted_bam"
  ],
  "snyderSimChr15": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15/hisat2.sam",
    "-r1", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_1.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_1.fq",
    "-r2", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_2.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_2.fq",
    "-N", "20",
    "-sorted_bam"
  ],
  "kidneySimChr15": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/kidneySimChr15/hisat2.sam",
    "-r1", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_1.fq",
    "-r2", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_2.fq",
    "-N", "20",
    "-sorted_bam"
  ], 
}

"""
Performance
"""
path_test_roc = {
  "alignment": "/data/shunfu/ref_shannon_modi/snyderSimChr15/hits.sorted.bam",
  # "alignment": "/data/shunfu/ref_shannon_modi/snyderSimChr15/tophat2.bam",
  "genomeFile": "/data/shunfu/ref_shannon_modi/genome/chr15.fa",

  "cases": [
    # ['stringtie', '', 'stringtie_DefaultParam'],
    # ['stringtie', '--maxSens -addHead', 'stringtie_f_0_c_0.001'],
    # ['cufflinks', '', 'cufflinks_DefaultParam'],
    # ['cufflinks', '--maxSens', 'cufflinks_F_0.001'],
    # ['strawberry', '', 'strawberry_DefaultParam'],
    # ['ryuto', '', 'ryuto_DefaultParam'],
    ['trinity', '', 'trinity_DefaultParam']
  ],

  "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15/", # + case[0]/case[2]
  # "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15_tophat2/", # + case[0]/case[2]

  "reference": "/data/shunfu/ref_shannon_modi/reference/ref_chr15.fasta",

  "nJobs": 20,
}
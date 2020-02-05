from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

"""
Modules
"""
path_test_aligner_hisat2 = {
  # "wwSimChr15Hisat2": [
  #   # ksreeram
  #   "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
  #   "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/wwSimChr15/hisat2.sam",
  #   "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/ww_snyderRef/sim_reads.fa",
  #   # small file for test
  #   # "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/test_sgRefShannon/reads_1.fa",
  #   "-fasta",
  #   "-N", "20",
  #   "-sorted_bam"
  # ],
  # "snyderSimChr15": [
  #   # ksreeram
  #   "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
  #   "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15/hisat2.sam",
  #   "-r1", 
  #   "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_1.fq",
  #   # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_1.fq",
  #   "-r2", 
  #   "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_2.fq",
  #   # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_2.fq",
  #   "-N", "20",
  #   "-sorted_bam"
  # ],
  # "kidneySimChr15": [
  #   # ksreeram
  #   "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
  #   "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/kidneySimChr15/hisat2.sam",
  #   "-r1", 
  #   "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_1.fq",
  #   "-r2", 
  #   "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_2.fq",
  #   "-N", "20",
  #   "-sorted_bam"
  # ], 
  "wwSimHg19Hisat2": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/wwSimHg19/hisat2.sam",
    "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/ww_snyderRef/sim_reads.fa",
    # small file for test
    # "-r", "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/test_sgRefShannon/reads_1.fa",
    "-fasta",
    "-N", "20",
    "-sorted_bam"
  ],
  "snyderSimHg19Hisat2": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimHg19/hisat2.sam",
    "-r1", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_1.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_1.fq",
    "-r2", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/sim_reads_2.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_2.fq",
    "-N", "20",
    "-sorted_bam"
  ],
  "kidneySimHg19Hisat2": [
    # ksreeram
    "-g", "/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa",
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/kidneySimHg19/hisat2.sam",
    "-r1", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_1.fq",
    "-r2", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_2.fq",
    "-N", "20",
    "-sorted_bam"
  ], 
}

path_test_split_bam = {
  # ksreeram, snyderSimChr15, star
  # 'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/hits.sorted.bam',
  # 'chroms':'all',
  # 'output_sam':True,
  # 'nJobs':20,
  # 'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/',
  #
  "wwSimHg19Hisat2": {
    # ksreeram, snyderSimHg19, hisat2
    'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSimHg19/hisat2.sorted.bam',
    'chroms':'chr15',
    'output_sam':True,
    'nJobs':20,
    'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/wwSimHg19_hisat2/',
  },
  # "snyderSimHg19Hisat2": {
  #   # ksreeram, snyderSimHg19, hisat2
  #   'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSimHg19/hisat2.sorted.bam',
  #   'chroms':'chr15',
  #   'output_sam':True,
  #   'nJobs':20,
  #   'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/',
  # },
  "kidneySimHg19Hisat2": {
    # ksreeram, snyderSimHg19, hisat2
    'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/kidneySimHg19/hisat2.sorted.bam',
    'chroms':'chr15',
    'output_sam':True,
    'nJobs':20,
    'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/kidneySimHg19_hisat2/',
  },
}

"""
Performance
"""
""" turing
path_test_refShannon_roc = {
  "snyderSimChr15": {
    "rootResDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15/refShannon/",
    "T": ['0'], # ['0', '0.25', '0.4', '0.96']
    "sam_file": "/data/shunfu/ref_shannon_modi/snyderSimChr15/hits.sam",
    "genomeFile": "/data/shunfu/ref_shannon_modi/genome/chr15.fa",
    "pairedStr": "-paired",
    "chrom": "chr15",
    "Tref": "/data/shunfu/ref_shannon_modi/reference/ref_chr15.fasta",
    "nJobs": 20,}
}
"""
# """ ksreeram
path_test_refShannon_roc = {
  "snyderSimChr15_hisat2": {
    "rootResDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSimChr15/hisat2/refShannon_roc/",
    "T": ['0', '0.25', '0.4', '0.96'],
    "sam_file": "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/chr15/hisat2.sorted.sam",
    "genomeFile": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
    "pairedStr": "-paired",
    "chrom": "chr15",
    "Tref": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/reference.fasta",
    "nJobs": 20,}
}
# """

""" turing
path_test_exAssembler_roc = {
  # "alignment": "/data/shunfu/ref_shannon_modi/snyderSimChr15/hits.sorted.bam",
  # "alignment": "/data/shunfu/ref_shannon_modi/snyderSimChr15/tophat2.bam",
  "alignment": "/data/shunfu/ref_shannon_modi/snyderSimChr15/hisat2.sorted.bam",
  "genomeFile": "/data/shunfu/ref_shannon_modi/genome/chr15.fa",

  "cases": [
    # ['stringtie', '', 'stringtie_DefaultParam'],
    # ['stringtie', '--maxSens -addHead', 'stringtie_f_0_c_0.001'],
    ['TransComb', '', 'TransComb_DefaultParam'],
    ['CLASS2', '', 'CLASS2_DefaultParam'],
    ['scallop', '', 'scallop_DefaultParam'],
    ['cufflinks', '', 'cufflinks_DefaultParam'],
    ['cufflinks', '--maxSens', 'cufflinks_F_0.001'],
    # ['strawberry', '', 'strawberry_DefaultParam'],
    # ['ryuto', '', 'ryuto_DefaultParam'],
    # ['trinity', '', 'trinity_DefaultParam']
  ],

  # "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15/", # + case[0]/case[2]
  # "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15_tophat2/", # + case[0]/case[2]
  "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15_hisat2/", # + case[0]/case[2]

  "reference": "/data/shunfu/ref_shannon_modi/reference/ref_chr15.fasta",

  "nJobs": 20,
}
"""
#""" ksreeram
path_test_exAssembler_roc = {
  # "alignment":  "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2/accepted_hits.bam",
  # "alignment":  "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/hits.sorted.bam",
  "alignment":  "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/chr15/hisat2.sorted.bam",
  "genomeFile": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",

  "cases": [
    # ['stringtie', '', 'stringtie_DefaultParam'],
    # ['stringtie', '--maxSens -addHead', 'stringtie_f_0_c_0.001'],
    # ['TransComb', '', 'TransComb_DefaultParam'],
    ['CLASS2', '--maxSens', 'CLASS2_F0'],
    # ['CLASS2', '', 'CLASS2_DefaultParam'],
    # ['scallop', '', 'scallop_DefaultParam'],
    # ['cufflinks', '', 'cufflinks_DefaultParam'],
    # ['cufflinks', '--maxSens', 'cufflinks_F_0.001'],
    # ['strawberry', '', 'strawberry_DefaultParam'],
    # ['ryuto', '', 'ryuto_DefaultParam'],
    # ['trinity', '', 'trinity_DefaultParam']
  ],

  # "resDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15/test_exAssembler_roc/", # + case[0]/case[2]
  # "resDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15_tophat2/test_exAssembler_roc/", # + case[0]/case[2]
  "resDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15_hisat2/test_exAssembler_roc/", # + case[0]/case[2]

  "reference": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/reference.fasta",

  "nJobs": 20,
}
#"""
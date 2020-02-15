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
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/wwSimHg19_2/hisat2.sam",
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
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimHg19_2/hisat2.sam",
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
    "-o", "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/kidneySimHg19_2/hisat2.sam",
    "-r1", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_1.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_1.fq",
    "-r2", 
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/kidney_snyderRef/sim_reads_2.fq",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/sim_star/snyder/lc_sim_reads_100k_2.fq",
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
  # 'ratio':None,
  # 'combine': False,
  # 'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/',
  #
  # "wwSimHg19Hisat2": {
  #   # ksreeram, snyderSimHg19, hisat2
  #   'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSimHg19/hisat2.sorted.bam',
  #   'chroms':'chr15',
  #   'output_sam':True,
  #   'nJobs':20,
  #   'ratio': None,
  #   'combine': False,
  #   'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/wwSimHg19_hisat2/',
  # },
  # "snyderSimHg19Star": {
  #   # ksreeram, snyderSimHg19, STAR
  #   # 'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/hits.sorted.bam',
  #   'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSim_1018a/hits.sorted.bam', # unsorted
  #   'chroms':'chr1,chr15',
  #   'output_sam':False,
  #   'nJobs':20,
  #   'ratio': 0.05,
  #   'combine': True,
  #   'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSim_Star_Chr1Chr3_Ratio0.05/',
  # },
  # "snyderSimHg19Hisat2": {
  #   # ksreeram, snyderSimHg19, hisat2
  #   # 'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSimHg19/hisat2.sorted.bam',
  #   'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimHg19_2/hisat2.bam', # unsorted
  #   'chroms':'chr15',
  #   'output_sam':True,
  #   'nJobs':20,
  #   'ratio': None,
  #   'combine': False,
  #   'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/',
  # },
  "snyderRealHg19Star": {
    'input_bam':'/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/hits.sorted.bam', # unsorted
    'chroms':'chr1,chr15',
    'output_sam':False,
    'nJobs':20,
    'ratio': 0.05,
    'combine': True,
    'outDir':'/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/',
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
  "snyderSimChr15_tophat2": {
    "rootResDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSimChr15/tophat2/refShannon_roc/",
    "T": ['0', '0.25', '0.4', '0.96'],
    "sam_file": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2/accepted_hits.sam",
    "genomeFile": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
    "pairedStr": "-paired",
    "chrom": "chr15",
    "Tref": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/reference.fasta",
    "nJobs": 20,}
  # "snyderSimChr15_hisat2": {
  #   "rootResDir": "/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSimChr15/hisat2/refShannon_roc/",
  #   "T": ['0', '0.25', '0.4', '0.96'],
  #   "sam_file": "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/chr15/hisat2.sam",
  #   "genomeFile": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",
  #   "pairedStr": "-paired",
  #   "chrom": "chr15",
  #   "Tref": "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/reference.fasta",
  #   "nJobs": 20,}
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

path_exAssembler_run = {
  'dummy_run': {
    'case': 
      ['stringtie', '', 'stringtie_DefaultParam'],
  #     # ['trinity', '', 'trinity_DefaultParam'],
    'alignment': 
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSim_Star_Chr15_Ratio0.05/hits.sorted.bam',
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSim_Star_Chr1Chr3_Ratio0.05//hits.sorted.bam',
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSim_1018a/hits.sorted.bam',
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/hits.sorted.bam', # dummy file
    'genomeFile':
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa',
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_run/snyderSim_Star_Chr15_Ratio0.05/',
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_run/trinity/',
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_run/snyderSim_Star_Hg19_All/trinity/',
  #     # '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/trinity/',
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/stringtie/',
    'nJobs': 20,
  },
  # 'ww_sim_star_all_ryuto': {
  #   'case': 
  #     ['ryuto', '', 'ryuto_DefaultParam'],
  #   'alignment': 
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/wwSim_1122a/hits.sorted.bam',
  #   'genomeFile':
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
  #   'resDir':
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSim_Star_Hg19_All/ryuto/',
  #   'nJobs': 20,
  # },
  # 'snyder_sim_star_all_ryuto': {
  #   'case': 
  #     ['ryuto', '', 'ryuto_DefaultParam'],
  #   'alignment': 
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSim_1018a/hits.sorted.bam',
  #   'genomeFile':
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
  #   'resDir':
  #     '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSim_Star_Hg19_All/ryuto/',
  #   'nJobs': 20,
  # },
  'kidney_sim_star_all_ryuto': {
    'case': 
      ['ryuto', '', 'ryuto_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/kidneySim/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/kineySim_Star_Hg19_All/ryuto/',
    'nJobs': 20,
  },
  'kidney_sim_star_all_trinity': {
    'case': 
      ['trinity', '', 'trinity_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/kidneySim/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/kineySim_Star_Hg19_All/trinity/',
    'nJobs': 20,
  },
  'ww_real_star_all_ryuto': {
    'case': 
      ['ryuto', '', 'ryuto_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/WingWongTest_K24_135M/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwReal_Star_Hg19_All/ryuto/',
    'nJobs': 20,
  },
  'ww_real_star_all_trinity': {
    'case': 
      ['trinity', '', 'trinity_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/WingWongTest_K24_135M/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwReal_Star_Hg19_All/trinity/',
    'nJobs': 20,
  },
  'snyder_real_star_all_ryuto': {
    'case': 
      ['ryuto', '', 'ryuto_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderReal_Star_Hg19_All/ryuto/',
    'nJobs': 20,
  },
  'snyder_real_star_all_trinity': {
    'case': 
      ['trinity', '', 'trinity_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderReal_Star_Hg19_All/trinity/',
    'nJobs': 20,
  },
  'kidney_real_star_all_ryuto': {
    'case': 
      ['ryuto', '', 'ryuto_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/kidney_0916a/_whole_sam/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/kidneyReal_Star_Hg19_All/ryuto/',
    'nJobs': 20,
  },
  'kidney_real_star_all_trinity': {
    'case': 
      ['trinity', '', 'trinity_DefaultParam'],
    'alignment': 
      '/data1/shunfu1/ref_shannon_modi/data/_copy/kidney_0916a/_whole_sam/hits.sorted.bam',
    'genomeFile':
      '/data1/shunfu1/ref_shannon_modi/data/_copy/snyder_0729a/genome/hg19.fa',
    'resDir':
      '/data1/shunfu1/ref_shannon_modi/data/_copy3/kidneyReal_Star_Hg19_All/trinity/',
    'nJobs': 20,
  },
}

path_test_gen_logs = {
  # 'case_dummy' : { # snyder, real, star, chr1chr15, sample rate 5%
  #     'Tref': 
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyder_0807b/reference.fasta',
  #     'Trec':
  #       # '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/trinity/trinity.fasta',
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/stringtie/stringtie.fasta',
  #     'resDir':
  #       # '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/trinity/',
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/stringtie/',
  #   },
  # 'case_ww_sim_all_trinity_default': { # ww, sim, star, all chroms, sample rate 100%
  #     'Tref': 
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyder_0807b/reference.fasta',
  #     'Trec':
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_run/wwSim_Star_Hg19_All/trinity/trinity.fasta',
  #     'resDir':
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSim_Star_Hg19_All/trinity/',
  # },
  'case_ww_sim_all_ryuto_default': { # ww, sim, star, all chroms, sample rate 100%
      'Tref': 
        '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyder_0807b/reference.fasta',
      'Trec':
        '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSim_Star_Hg19_All/ryuto/ryuto.fasta',
      'resDir':
        '/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSim_Star_Hg19_All/ryuto/',
  },
  # 'case_snyder_sim_all_trinity_default': { # snyder, sim, star, all chroms, sample rate 100%
  #     'Tref': 
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyder_0807b/reference.fasta',
  #     'Trec':
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_run/snyderSim_Star_Hg19_All/trinity/trinity.fasta',
  #     'resDir':
  #       '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSim_Star_Hg19_All/trinity/',
  # },
  'case_snyder_sim_all_ryuto_default': { # snyder, sim, star, all chroms, sample rate 100%
      'Tref': 
        '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyder_0807b/reference.fasta',
      'Trec':
        '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSim_Star_Hg19_All/ryuto/ryuto.fasta',
      'resDir':
        '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyderSim_Star_Hg19_All/ryuto/',
  },
}

path_test_gen_sens = {
  'case_dummy' : { # snyder, real, star, chr1chr15, sample rate 5%
      'IM_list': [[1,100]],
      'oracleFa': '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/reference/snyder/oracle/snyder_oracle_et25_ct100_cut0.fa',
      'numIsoFile': '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/reference/snyder/geneIsoNum.txt',
      'recLog': '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/stringtie/stringtie_log.txt',
      'cmpLog': '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/trinity/trinity_log.txt',
      'expFile': '/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/rsem/exp_star/snyder/exp.isoforms.results',
      'expFormat': 'RSEM',
      'L': 101,
      'S': 0,
      'outFileStem': '/data1/shunfu1/ref_shannon_modi/data/_copy3/snyder_real_star_chr1chr15_Ratio0.05/sens_stringtie_trinity',
      'ablist': '0.0,2.23,6.09,17.83,56.85,1e10',
    },
}

#""" ksreeram
path_test_exAssembler_roc = {
  "alignment":
    # ww sim chr15
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/wwSimChr15/hits.sorted.bam",
    # snyder sim chr15
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2/accepted_hits.bam",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/hits.sorted.bam",
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_split_bam/snyderSimHg19_hisat2/chr15/hisat2.sorted.bam",
    # kidney sim chr15
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/kidneySimChr15/hits.sorted.bam",

  "genomeFile":
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/genome/human/chr15.fa",

  "cases": [
    # ['stringtie', '', 'stringtie_DefaultParam'],
    # ['stringtie', '--maxSens -addHead', 'stringtie_f_0_c_0.001'],
    # ['TransComb', '', 'TransComb_DefaultParam'],
    ['CLASS2', '--maxSens', 'CLASS2_F0'],
    ['CLASS2', '', 'CLASS2_DefaultParam'],
    ['scallop', '', 'scallop_DefaultParam'],
    # ['cufflinks', '', 'cufflinks_DefaultParam'],
    # ['cufflinks', '--maxSens', 'cufflinks_F_0.001'],
    ['strawberry', '', 'strawberry_DefaultParam'],
    ['ryuto', '', 'ryuto_DefaultParam'],
    ['trinity', '', 'trinity_DefaultParam']
  ],

  "resDir": 
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15/test_exAssembler_roc/", # + case[0]/case[2]
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15_tophat2/test_exAssembler_roc/", # + case[0]/case[2]
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_aligner_hisat2/snyderSimChr15_hisat2/test_exAssembler_roc/", # + case[0]/case[2]
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/test_exAssembler_roc/cgmemtime/snyderSimChr15/", # + case[0]/case[2]
    # ww sim chr15
    # "/data1/shunfu1/ref_shannon_modi/data/_copy3/wwSim_Star_Chr15/",
    # kidney sim chr15
    "/data1/shunfu1/ref_shannon_modi/data/_copy3/kidneySim_Star_Chr15/",

  "reference":
    # ww 7703 refs
    # "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/wwSimChr15/reference.fasta",
    # snyder 7703 refs
    "/data1/shunfu1/ref_shannon_modi/data/_copy2/sgRefShannon/snyderSimChr15/reference.fasta",

  "nJobs": 20,
}
#"""
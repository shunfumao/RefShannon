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
    # ['strawberry', '--maxSens', 'strawberry_DefaultParam'],
    []
  ],

  "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15/", # + case[0]/case[2]
  # "resDir": "/data/shunfu/ref_shannon_modi/results/roc/snyderSimChr15_tophat2/", # + case[0]/case[2]

  "reference": "/data/shunfu/ref_shannon_modi/reference/ref_chr15.fasta",

  "nJobs": 20,
}
from util import *

ref_gtf = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2_gffcmp/gencode.v15_PB_enhanced_HOP_Gm12878.ExonLinesOnly.gtf'
Dir = '/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/MoreExtAssemblers/SnyderSimChr15/tophat2_gffcmp/'

targets = []
targets.append(['%s/reconstructed.gtf'%Dir, '%s/reconstructed/'%Dir, 'reconstructed'])
targets.append(['%s/stringtie.gtf'%Dir, '%s/stringtie/'%Dir, 'stringtie'])
targets.append(['%s/TransComb.gtf'%Dir, '%s/TransComb/'%Dir, 'TransComb'])

for target in targets:
	print('target:')
	print(target)
	print('')

	cmd = 'mkdir -p %s'%target[1]
	run_cmd(cmd)

	cmd = 'gffcompare -R -r %s -o %s/%s %s'%(ref_gtf, target[1], target[2], target[0])
	run_cmd(cmd)
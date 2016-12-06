import subprocess
import os
import pdb
import math
from util import run_cmd
from global_values import *
import run_parallel_cmds

def cut_file(in_name,out_name,line_start,line_end):
	run_cmd('awk \'NR > ' + str(line_end) + ' { exit } NR >= ' + str(line_start) +  '\' '+ in_name + ' > ' + out_name )

def parallel_blat(target_fasta,query_fasta,out_file,QUERY_SPLIT):
	'''Function takes in target,query and output file. parallelizes blat by running GNU parallel
	- Currently only parallelizes on query space
	- Also assumes that query fasta file takes two lines per sequence (not wrapped)'''
	target_length = float(subprocess.check_output('grep -c \'>\' ' + target_fasta,shell=True))
	query_length = float(subprocess.check_output('grep -c \'>\' ' + query_fasta,shell=True))
	run_cmd( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +query_fasta + ' > '+ query_fasta +'_nospace')
	#run_cmd( 'awk \'/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}\' ' +target_fasta + ' > '+ target_fasta)
	query_fasta = query_fasta + '_nospace'
	#TARGET_SPLIT = 1
	#QUERY_SPLIT = 4
	#Alernately
	#QUERY_SPLIT = min(int(math.ceil(float(query_length)/float(target_length))),50)
	#QUERY_SPLIT = max(int(math.ceil(float(query_length)/float(target_length))),500) 
	#QUERY_SPLIT = int(min(QUERY_SPLIT,query_length))
	#QUERY_SPLIT= 100
	#pdb.set_trace()
	print('Query length: ' +str(query_length) + ' Target length: ' + str(target_length) + ' Query Split: ' + str(QUERY_SPLIT))
	split_size = int(math.floor(float(query_length)/QUERY_SPLIT))
	'''if split_size % 2 !=0:
		split_size +=1'''
	'''if query_length <= float(target_length):
		print('Cannot parallelize on query. Running Vanilla Blat')
		run_cmd('blat -noHead ' + target_fasta + ' ' +  query_fasta + ' ' + out_file)	
		return'''

	for n in range(QUERY_SPLIT):
		if n==QUERY_SPLIT-1:
                	cut_file(query_fasta,query_fasta+'_'+str(n+1),2*(n)*split_size+1,2*query_length)
		else:
			cut_file(query_fasta,query_fasta+'_'+str(n+1),2*(n)*split_size+1,2*(n+1)*split_size)
	#pdb.set_trace()
	q_range = range(QUERY_SPLIT)
	#x = [int(i)+1 for i in q_range]
	#q_str = " ".join(map(str,x))
	run_cmd('rm ' + out_file + '_*' )
	#print('parallel --jobs %d blat -noHead '%MAX_PARALLEL_PROCESS + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str )
	#run_cmd('time parallel --jobs %d blat -noHead '%MAX_PARALLEL_PROCESS + target_fasta + ' ' + query_fasta + '_{} ' +out_file + '_{} ::: ' + q_str  )
	cmds = []
	for i in q_range:
		cmd = 'blat -noHead ' + target_fasta + ' ' + query_fasta + '_%d '%i +out_file + '_%d'%i
		cmds.append(cmd)
	run_parallel_cmds.run_cmds(cmds, MAX_PARALLEL_PROCESS)

	run_cmd('cat ' + out_file + '_* > ' + out_file)
	run_cmd('rm ' + out_file + '_*' )
	run_cmd('rm ' + query_fasta + '_*' )

def main():
    import sys
    if len(sys.argv) == 5:
	Query_split = int(sys.argv[4])
    else: 
	Query_split = 100
    parallel_blat(sys.argv[1],sys.argv[2],sys.argv[3],Query_split)#rec -- 1,ref -- 2


if __name__ == '__main__':
    main()


import os, sys, pdb

def run_cmd(s1):
    print(s1); 
    os.system(s1) # + '>> temp_out.txt')


def write_filtered_tr(depth_file, in_tr_file, out_tr_file, log_file, threshold=0.9):
    import sys, pdb
    tr_hits = {}
    THRESH = threshold
    for a in open(depth_file): # as depth_file:
        #a=depth_file.readlines()
        fields = a.strip().split()
        tr_hits[fields[0]] = tr_hits.get(fields[0],0) +1;

    tr_len = {}; tr_name = ''
    with open(out_tr_file,'w') as tr_file, open(log_file,'w') as lf: #output transcripts
       for line in open(in_tr_file): #in_tr_file
            fields = line.strip().split()
            if fields[0][0] == '>': tr_name = fields[0][1:]; continue
            tr_len[tr_name]=len(fields[0]); trh = tr_hits.get(tr_name,0)
            lf.write(tr_name + '\t' + str(trh) + '\t' + str(len(fields[0])) + '\n')
            if tr_hits.get(tr_name,0) >= len(fields[0]) * THRESH:
                tr_file.write('>' + tr_name + '\n')
                tr_file.write(fields[0] + '\n')

def filter_FP(rec_fasta, read_1, read_2, out_dir, flags = '-f --ff', useHisat2 = False):

    #pdb.set_trace()
    #take rec_fasta, read_1, read_2 files of type. 
    #Output in out_fasta. Temp dir = out_dir. 
    #Output fasta is in out_dir/reconstructed.fasta 
    #Flags is '-f --ff' for fasta and forward-forward
    #Flags is '-q --fr' for fastq and forward-reverse
    hisat_dir = '' #include dir/ if specifying directory

    if useHisat2==True:
        ver = '2'
    else:
        ver = ''
    
    #Build hisat file
    run_cmd(hisat_dir + 'hisat%s-build '%ver + rec_fasta + ' ' + out_dir + '/rec.hisat')

    #Align
    run_cmd(hisat_dir + 'hisat%s'%ver + ' -p 20 --no-spliced-alignment --no-discordant '+ flags + ' -x ' + out_dir + '/rec.hisat -1 ' + read_1 +  ' -2 ' + read_2 + ' -S ' + out_dir + '/rec.sam' )
    
    #Process SAM / BAM file to get depth information
    run_cmd('samtools view -@ 20 -bS -f 0x2 ' + out_dir + '/rec.sam > ' + out_dir + '/rec.bam')
    run_cmd('samtools sort -@ 20 ' + '-o ' + out_dir + '/rec_sort.bam ' +  out_dir + '/rec.bam ')
    run_cmd('samtools depth ' +  out_dir + '/rec_sort.bam > ' +  out_dir + '/rec.depth')

    #Use the BAM file along with our tool to get new fasta file
    depth_file = out_dir + '/rec.depth'; 
    in_tr_file = rec_fasta
    out_tr_file = out_dir + '/rec.fasta'; 
    log_file = out_dir + '/rec.log'
    write_filtered_tr(depth_file, in_tr_file, out_tr_file, log_file)
    run_cmd('cp ' + rec_fasta + ' ' + out_dir +'/reconstructed_org.fasta')
    run_cmd('mv ' + out_tr_file + ' ' + out_dir + '/reconstructed.fasta')

'''
usage:

python filter_FP.py rec_fasta read_1 read_2 out_dir -f/q --fr/rf/ff [-hisat2]

# note
if a transcript is re-covered by over 90%, output it

'''
if __name__ == '__main__':
    rec_fasta = sys.argv[1]
    read_1 = sys.argv[2]
    read_2 = sys.argv[3]
    out_dir = sys.argv[4]
    run_cmd('mkdir -p %s'%out_dir)
    flags = '%s %s'%(sys.argv[5], sys.argv[6]) #-f --fr etc

    useHisat2=False
    if '-hisat2' in sys.argv:
        useHisat2=True
    filter_FP(rec_fasta, read_1, read_2, out_dir, flags, useHisat2)



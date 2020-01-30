import pdb, sys
from RefShannon.util import *
import numpy

'''
Note:

#### Description

Originally from Find_Covered_Transcripts.py.

This function prepares oracle set (a subset of reference transcripts that are fully covered by reads) from reference transcripts and reads.

Input:

- coverage_file ("tr genome_pos cov") indicates per genome pos per tr, what is the read coverage. can be obtained by reads aligned onto transcripts and samtools depth

- fasta_file: reference transcripts

- out_file: in fasta format, storing the oracle set

parameters:

- end_threshold (et): e.g. 25, focus on tr[et:trLen-et] region or tr[0:trLen] if trLen is too small
- coverage_threshold (ct): e.g. <=100. At least ct/100 percentage of tr (excluding et) should be covered (ct/100<1 means there's intermediate gaps)
- consecutive_uncovered_threshold (cut): e.g. 0. len of <= cut gap allowed in tr.
- a transcript is considered to be in oracle set if tr (excluding et) is highly covered (>= ct/100) with gaps <= cut allowed  

#### coverage file generation
- coverage file format: "tr genome_pos cov"
- procedure
  utilize rsem result (rsem-calculate-expression estimate transcripts abundance from reads) --> (exp).transcript.bam
  samtools sort -o (exp).transcript.sorted.sam -@ (20) (exp).transcript.bam
  samtools depth (exp).transcript.sorted.sam > coverage.txt

'''

def Find_Covered_Transcripts(end_threshold, coverage_threshold, consecutive_uncovered_threshold, coverage_file_name, fasta_file_name, output_file_name):
    coverage_file = open(coverage_file_name, 'r')
    fasta_file = open(fasta_file_name, 'r')
    lines1 = coverage_file.readlines()
    lines2 = fasta_file.readlines()
    
    nLines=len(lines1); T=nLines/100; p=0; q=0;
    transcript_dict = {}
    for line in lines1:
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (lines1)'%q); sys.stdout.flush()

        tokens = line.split()
        transcript = tokens[0]
        if transcript in transcript_dict:
            if int(tokens[2]) >= 1:
                transcript_dict[transcript].append(int(tokens[1])-1)
        else:
            if int(tokens[2]) >= 1:
                transcript_dict[transcript] = [(int(tokens[1])-1)]   
    print('')         
    
    #pdb.set_trace()
    
    Transcript_Strings = {}
    Transcript_Lengths = {}
    Transcript_Fasta_Line = {}
    
    i = 0
    nLines=len(lines2); T=nLines/100; p=0; q=0;
    for line2 in lines2:
        p += 1
        if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (lines2)'%q); sys.stdout.flush()

        tokens = line2.split()
        #print(tokens)
        if len(line2.split()) == 0:
           break
        if tokens[0][0] == ">":
            ID = tokens[0][1:]
            Transcript_Fasta_Line[ID] = line2
            string = ""
            #pdb.set_trace()
            j = i + 1
            while j < len(lines2):
                cand = lines2[j]
                if len(cand.split()) == 0:
                    break
                if cand.split()[0][0] == ">":
                    break
                j += 1
                string = string + cand.split()[0]
            Transcript_Strings[ID] = string
            Transcript_Lengths[ID] = len(string)
            #print(ID)
        i = i + 1  
    print('')

    #pdb.set_trace()    
    Reconstructable_Transcripts = []
    if float(coverage_threshold) == 100.0:
        nLines=len(transcript_dict); T=nLines/100; p=0; q=0;
        for transcript in transcript_dict:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (tr filtering)'%q); sys.stdout.flush()
            reconstructable_transcript = True
            #pdb.set_trace()
            if end_threshold >= float(Transcript_Lengths[transcript])/2.0:
                for ind in range(0, Transcript_Lengths[transcript]):
                    if ind not in transcript_dict[transcript]:
                        reconstructable_transcript = False   
                        break
            else:
                for ind in range(end_threshold, Transcript_Lengths[transcript] - end_threshold):
                    if ind not in transcript_dict[transcript]:
                        reconstructable_transcript = False
                        break
            if reconstructable_transcript == True:
                Reconstructable_Transcripts.append(transcript)
    else:
        nLines=len(transcript_dict); T=nLines/100; p=0; q=0;
        for transcript in transcript_dict:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (tr filtering)'%q); sys.stdout.flush()
            reconstructable_transcript = True
            #pdb.set_trace()
            num_covered_bases = 0
            if end_threshold >= float(Transcript_Lengths[transcript])/2.0:
                consecutive_uncovered = 0
                for ind in range(0, Transcript_Lengths[transcript]):
                    if consecutive_uncovered >= consecutive_uncovered_threshold:
                        reconstructable_transcript = False
                        break
                    if ind in transcript_dict[transcript]:
                        num_covered_bases += 1
                        consecutive_uncovered = 0
                    else:
                        consecutive_uncovered += 1
                if float(num_covered_bases)/float(Transcript_Lengths[transcript]) < float(coverage_threshold)/100.0:
                    reconstructable_transcript = False                
            else:
                consecutive_uncovered = 0
                for ind in range(end_threshold, Transcript_Lengths[transcript] - end_threshold):
                    if consecutive_uncovered >= consecutive_uncovered_threshold:
                        reconstructable_transcript = False
                        break
                    if ind in transcript_dict[transcript]:
                        num_covered_bases += 1
                        consecutive_uncovered = 0
                    else:
                        consecutive_uncovered += 1
                if float(num_covered_bases)/float(Transcript_Lengths[transcript] - end_threshold - end_threshold) < float(coverage_threshold)/100.0:
                    reconstructable_transcript = False   
            if reconstructable_transcript == True:
                Reconstructable_Transcripts.append(transcript)        
    #pdb.set_trace()
    print('')
    
    File = open(output_file_name, 'w')
    i = 0
    #File.write
    for transcript in Reconstructable_Transcripts:
        File.write(Transcript_Fasta_Line[transcript] )
        File.write(Transcript_Strings[transcript] + "\n")
        i += 1
    File.close()

    coverage_file.close()
    fasta_file.close()

#Find_Covered_Transcripts(100, 80, 100, "reads_pacbio.depth", "reference_PacBio.fasta", "reads_pacbio_oracle_set.fasta")
#Find_Covered_Transcripts(100, 100, 0, "./WingWongTest_algo_input/Bowtie/reads_GSE_depth.txt", "./WingWongTest_algo_output/reference_GSE_unique.fasta", "./WingWongTest_algo_output/GSE_oracle_set.fasta")
#Find_Covered_Transcripts(50,100,0,"/data/sreeramk/Full_Assembler/Snyder_algo_input/depth_chr15.txt","/data/sreeramk/Full_Assembler/Snyder_algo_input/reference_chr15.fasta","/data/sreeramk/Full_Assembler/Snyder_algo_input/oracle_set_chr15.fasta")         

def findCoveredTr(args):

    end_threshold = int(args[args.index('--et')+1])
    coverage_threshold = int(args[args.index('--ct')+1])
    consecutive_uncovered_threshold = int(args[args.index('--cut')+1])
    coverage_file_name = args[args.index('--trCov')+1]
    fasta_file_name = args[args.index('--trFa')+1]
    output_file_name = args[args.index('-o')+1]

    Find_Covered_Transcripts(end_threshold, \
                             coverage_threshold, \
                             consecutive_uncovered_threshold, \
                             coverage_file_name, \
                             fasta_file_name, \
                             output_file_name)
    return
       
def findTrGeneNum(args):

    infile = args[args.index('-i')+1] #path/to/rsem_exp_isoforms_results
    numLines = sum([1 for line in open(infile, 'r')]); T=numLines/100

    outfile = args[args.index('-o')+1] #path/to/numIsoFile

    out_dir = parent_dir(outfile)
    run_cmd('mkdir -p %s'%out_dir)

    geneIsoCnt = {} # key - gene_id val - num of isoforms
    with open(infile, 'r') as f:
        p=0;q=0;
        f.readline()
        #pdb.set_trace()
        for line in f:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (build geneIsoCnt)'%q); sys.stdout.flush()
            tokens = line.split()
            tid = tokens[0]
            gid = tokens[1]
            if gid not in geneIsoCnt:
                geneIsoCnt[gid]=0
            geneIsoCnt[gid]+=1

        print('')

    #pdb.set_trace()

    with open(infile, 'r') as fi, \
         open(outfile, 'w') as fo:

        p=0;q=0;
        fi.readline()
        #pdb.set_trace()
        for line in fi:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (export)'%q); sys.stdout.flush()
            tokens = line.split()
            tid = tokens[0]
            gid = tokens[1]
            st = '%s\t%s\t%d\n'%(tid, gid, geneIsoCnt[gid])
            fo.write(st)
        
        print('')            

    return

def findTrGeneNum_f(args):

    infile = args[args.index('-f')+1] #path/to/refFa
    numLines = sum([1 for line in open(infile, 'r')]); T=numLines/100

    outfile = args[args.index('-o')+1] #path/to/numIsoFile

    out_dir = parent_dir(outfile)
    run_cmd('mkdir -p %s'%out_dir)

    geneIsoCnt = {} # key - gene_id val - num of isoforms
    with open(infile, 'r') as f:
        p=0;q=0;
        #f.readline()
        pdb.set_trace()
        for line in f:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (build geneIsoCnt)'%q); sys.stdout.flush()
            
            if line[0]!='>': continue

            tokens = line.split()

            tid = tokens[0][1:] #skip '>'
            gid = tokens[1][5:] #skip 'GENE='
            if gid not in geneIsoCnt:
                geneIsoCnt[gid]=0
            geneIsoCnt[gid]+=1

        print('')

    #pdb.set_trace()

    with open(infile, 'r') as fi, \
         open(outfile, 'w') as fo:

        p=0;q=0;
        #fi.readline()
        pdb.set_trace()
        for line in fi:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (export)'%q); sys.stdout.flush()

            if line[0]!='>': continue

            tokens = line.split()
            tid = tokens[0][1:]
            gid = tokens[1][5:]
            st = '%s\t%s\t%d\n'%(tid, gid, geneIsoCnt[gid])
            fo.write(st)
        
        print('')            

    return

def performancePlotExpress(args): # a Wrapper

    oracle_fasta = args[args.index('--oracleFa')+1]
    num_isoform_file = args[args.index('--numIsoFile')+1]
    reconstr_log = args[args.index('--recLog')+1]
    stringtie_log = args[args.index('--cmpLog')+1]
    plot_file = args[args.index('--o')+1]
    exp_file = args[args.index('--expFile')+1]
    exp_format = args[args.index('--expFormat')+1]
    L = int(args[args.index('-L')+1])
    S = int(args[args.index('-S')+1])
    iso_low = int(args[args.index('--isoLow')+1])
    iso_high = int(args[args.index('--isoHigh')+1])

    if '--ablist' in args:
        ab_str = args[args.index('--ablist')+1]
        ab_list = [float(itm) for itm in ab_str.split(',') if itm != '']
    else:
        ab_list = []

    performance_plot_express_modi(oracle_fasta, \
                                  num_isoform_file, \
                                  reconstr_log, \
                                  stringtie_log, \
                                  plot_file, \
                                  exp_file, \
                                  exp_format, \
                                  L, \
                                  S, \
                                  iso_low, \
                                  iso_high, \
                                  use_oracle=True, \
                                  ab_list=ab_list)
    return

#from tester.performance_plot_express
'''
S: minTrLen
'''
def performance_plot_express_modi(oracle_fasta, \
                                  num_isoform_file, \
                                  reconstr_log, \
                                  stringtie_log, \
                                  plot_file, \
                                  exp_file, \
                                  exp_format, \
                                  L, \
                                  S, \
                                  iso_low, \
                                  iso_high, \
                                  use_oracle=True, \
                                  ab_list=[]):
    tr_dict_full = {}
    
    if exp_format=='KAL':
        tot_ab = 0
        #pdb.set_trace()
        for lines in open(exp_file):
            fields=lines.strip().split();
            tr_id = fields[0].upper()
            if not fields[1].isdigit():
                continue;
            if float(fields[1]) < S:
                continue;
            try:
                tr_len = float(fields[1]);
                tr_ab = float(fields[4])/1e6;
                tr_cov = float(fields[3]) / max(1,(tr_len-L)) * L # * 20 / 135
                tot_ab = tot_ab + tr_ab
                tr_dict_full[tr_id] = [tr_len, tr_cov, tr_ab, 0, 0, 0] #last 4 entries are num_iso, ours, stringtie
            except ValueError:
                pdb.set_trace()
                continue

    elif exp_format=='RSEM':
        tot_ab = 0
        #pdb.set_trace()
        for lines in open(exp_file):
            fields=lines.strip().split();
            tr_id = fields[0].upper()
            if not fields[2].isdigit():
                continue;
            if float(fields[2]) < S:
                continue;
            try:
                tr_len = float(fields[2]);
                tr_ab = float(fields[5])/1e6; # sum(TPM_i)=1e6
                tr_cov = float(fields[4]) / max(1,(tr_len-L)) * L # expected_cnt = sum_i(prob(read_i from this transcript))
                tot_ab = tot_ab + tr_ab
                tr_dict_full[tr_id] = [tr_len, tr_cov, tr_ab, 0, 0, 0] #last 4 entries are num_iso, ours, stringtie
            except ValueError:
                pdb.set_trace()
                continue

    '''Filter to files in oracle_fasta'''
    tr_os = {}
    #pdb.set_trace()
    for lines in open(oracle_fasta):
        fields=lines.strip().split()
        tr_id = fields[0][1:].upper()
        if fields[0][0]=='>':
            if tr_dict_full.get(tr_id):
                tr_os[tr_id] = tr_dict_full[tr_id]

    if not use_oracle:
       tr_os = tr_dict_full; #Use to bypass filtering

    #pdb.set_trace()
    '''Filter based on number of isoforms'''
    #No of isoforms is specified 
    iso_lower = iso_low;
    iso_upper = iso_high; #consider only transcripts that have iso_lower<=num_iso <= iso_upper
    tr_dict = {}
    tot_ab = 0
    tot_no = 0

    if ab_list==[]:
        ab_list = [0,1,10,25,50,100,1e10] #actually used to compare tr coverage
    no = [0] * len(ab_list)

    stat_num_iso_os = []
    stat_cov_os = []

    for lines in open(num_isoform_file):
        fields=lines.strip().split()
        tr_id = fields[0].upper()
        if tr_os.get(tr_id):
            num_iso = float(fields[2])
            stat_num_iso_os.append(num_iso)            
            if iso_lower<=num_iso and num_iso <= iso_upper:
                stat_cov_os.append(tr_os[tr_id][1])
                tr_dict[tr_id] = tr_os[tr_id]
                tr_dict[tr_id][3] = num_iso
                tot_no += 1
                tot_ab += tr_dict[tr_id][2];
                tr_cov = tr_dict[tr_id][1]
                for (i,ab) in enumerate(ab_list):
                    if i >= len(ab_list)-1:
                        continue;
                    if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                        no[i] +=1

    #pdb.set_trace()
    #h_numIso, b_numIso = numpy.histogram(stat_num_iso_os, bins=20)
    #print('h_numIso=%s, b_numIso=%s'%(str(h_numIso), str(b_numIso)))
    
    stat_num_iso_os.sort(key=lambda x:x)
    st = 'percentile\tnumber of isoforms (oracle set)'
    print(st)
    for cutoff in range(10): #0,1,...,9
        stt = int(cutoff * 0.1 * len(stat_num_iso_os))
        stp = int((cutoff+1) * 0.1 * len(stat_num_iso_os))
        st = '%d %% ~ %d %%:\t'%(cutoff*10, (cutoff+1)*10)
        sub_list = stat_num_iso_os[stt:min(stp, len(stat_num_iso_os))]
        st += '%d ~ %d'%(sub_list[0], sub_list[len(sub_list)-1])
        print(st)
    #pdb.set_trace()

    stat_cov_os.sort(key=lambda x:x)
    st = 'percentile\tabundance (isoforms %d - %d)'%(iso_lower, iso_upper)
    print(st)
    for cutoff in range(10):
        stt = int(cutoff * 0.1 * len(stat_cov_os))
        stp = int((cutoff+1) * 0.1 * len(stat_cov_os))
        st = '%d %% ~ %d %%:\t'%(cutoff*10, (cutoff+1)*10)
        sub_list = stat_cov_os[stt:min(stp, len(stat_cov_os))]
        try:
            st += '%.5f ~ %.5f'%(sub_list[0], sub_list[len(sub_list)-1])
            print(st)
        except:
            continue        

    #pdb.set_trace()
    

    #h_cov, b_cov = numpy.histogram(stat_cov_os, bins=20)
    #print('h_ab=%s, b_ab=%s'%(str(h_cov), str(b_cov)))


    '''Run through our reconstructed file'''
    tot = [0] * len(ab_list)
    frac = [0] * len(ab_list)
    length = [0] * len(ab_list)
    ab_rec = 0
    no_rec = 0

    for lines in open(reconstr_log):
        fields=lines.strip().split();
        tr_id = fields[0].upper()
        if (not tr_dict.get(tr_id)):
            continue

        tr_rec = float(fields[2]);
        tr_len = tr_dict[tr_id][0];
        tr_cov = tr_dict[tr_id][1];
        tr_ab = tr_dict[tr_id][2]

        for (i,ab) in enumerate(ab_list):
            if i >= len(ab_list)-1:
                continue;
        
            if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]:
                #no[i] +=1
                length[i] += tr_len
                tot[i] += min(tr_rec,tr_len)
                if tr_rec >= 0.9*tr_len:
                    tr_dict[tr_id][4]=1
                    frac[i]+=1
                    ab_rec += tr_ab
                    no_rec += 1

    #pdb.set_trace()
    '''Run through stringtie'''
    stringtie_tot = [0] * len(ab_list)
    stringtie_frac = [0] * len(ab_list)
    #stringtie_no = [0] * len(ab_list)
    #stringtie_length = [0] * len(ab_list)
    ab_stringtie = 0
    no_stringtie = 0

    for lines in open(stringtie_log):
        fields=lines.strip().split();
        tr_id = fields[0].upper()
        if (not tr_dict.get(tr_id)):
            continue
            #pdb.set_trace()
        tr_rec = float(fields[2]);
        tr_len = tr_dict[tr_id][0];
        tr_cov = tr_dict[tr_id][1];
        tr_ab = tr_dict[tr_id][2]
        for (i,ab) in enumerate(ab_list):
            #tr_cov = tr_ab*L / norm * N
            #pdb.set_trace() 
            if tr_cov>=ab_list[i] and tr_cov<ab_list[i+1]: #stringtie_no[i] +=1
                #stringtie_length[i] += tr_len
                stringtie_tot[i] += tr_rec 
                if tr_rec >= 0.9*tr_len:
                    tr_dict[tr_id][5]=1
                    stringtie_frac[i]+=1
                    ab_stringtie += tr_ab
                    no_stringtie += 1
 
    with open(plot_file,'w') as plotFile: 
        info_str = ''
        info_str += '# oracle tr: %s\n'%oracle_fasta
        info_str += '# num iso file: %s\n'%num_isoform_file
        info_str += '# log1: %s\n'%reconstr_log
        info_str += '# log2: %s\n'%stringtie_log
        info_str += '# exp file: %s\n'%exp_file
        info_str += '# \n'
        info_str += '# readLength = %d\n'%L
        info_str += '# minTrLength = %d\n'%S
        info_str += '# min iso number of each isoform\'s gene = %d\n'%iso_lower
        info_str += '# max iso number of each isoform\'s gene = %d\n#\n'%iso_high

        plotFile.write(info_str)

        #plotFile.write('#Abundance \t No. Trans \t Frac. reconstructed (Ours) \t Frac. reconstructed (Stringtie) \t Bases reconstructed (Ours) \t Bases reconstr (Stringtie) \t No bases\n')
        plotFile.write('#Abundance \t No. Trans (Oracle) \t No. reconstructed (RefShannon) \t No. reconstructed (Stringtie)\n')
        for (i,ab) in enumerate(ab_list):
            if i==len(ab_list)-1: continue
            #plotFile.write(str(ab)+'\t'+str(no[i])+'\t'+str(frac[i])+'\t'+str(stringtie_frac[i])+'\t'+str(tot[i])+'\t'+str(stringtie_tot[i])+'\t' + str(length[i])+'\n')
            plotFile.write('[%s,%s)'%(str(ab), str(ab_list[i+1]))+'\t'+str(no[i])+'\t'+str(frac[i])+'\t'+str(stringtie_frac[i])+'\n')

    '''
    with open(plot_file+'_all','w') as plotFile:  
        for (tr_id, tr_info) in tr_dict.items():
            plotFile.write(str(tr_id)+'\t'+str(tr_info[0])+'\t'+str(tr_info[1])+'\t'+str(tr_info[2])+'\t'+str(tr_info[3])+'\t'+str(tr_info[4])+'\t' + str(tr_info[5])+'\n')
    '''
    #print('(ab) Total,Ours,Stringtie='+str(tot_ab)+','+str(ab_rec)+','+str(ab_stringtie))
    #print('(no) Total,Ours,Stringtie='+str(tot_no)+','+str(no_rec)+','+str(no_stringtie))

    return

'''
Usage:

(1) python analyze_realData_sensitivity.py  --findCoveredTr \
                                   --et end_threshold \
                                   --ct coverage_threshold \
                                   --cut consecutive_uncovered_threshold \
                                   --trCov path/to/coverage_file_name \
                                   --trFa path/to/fasta_file_name \
                                   -o path/to/output_file_name

    See Find_Covered_Transcripts() for description

(2) python analyze_realData_sensitivity.py  --findTrGeneNum \
                                       [-i path/to/rsem_exp_isoforms_results or -f path/to/refFa] \
                                       -o path/to/numIsoFile

    We use rsem-calculate-expression's (sample).isoforms.results to get numIsoFile
    (sample).isoforms.results = "transcript_id gene_id length effective_length expected_count TPM FPKM IsoPct"
    numIsoFile = "transcript_id gene_id num_of_transcripts_of_same_gene". To be used by performance_plot_express()

    In case rsem res not availabe (e.g. due to errors in refGTF file), we try to build numIsoFile info from refFa.
    We assume transcript description pattern ">tid GENE=gid"

(3) python analyze_realData_sensitivity.py --performancePlotExpress \
                                    --oracleFa oracle_fasta \
                                    --numIsoFile num_isoform_file \
                                    --recLog reconstr_log \
                                    --cmpLog stringtie_log \
                                    --o plot_file \
                                    --expFile exp_file \
                                    --expFormat exp_format \
                                    -L L \
                                    -S S \
                                    --isoLow iso_low \
                                    --isoHigh iso_high \
                                    [--ablist abLvl_0,abLvl_1,abLvl_2,...,abLvl_k] # [0~1] [1~2] ... [k-1~k]

    missing refTrLen at either reconstr_log or stringtie_log don't affect the processing.

Examples:

(1) python analyze_realData_sensitivity.py --findCoveredTr --et 25 --ct 100 --cut 0 --trCov /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/coverage.txt --trFa /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref.transcripts.fa -o /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/oracle/snyder_oracle_et25_ct100_cut0.fa

(2) python analyze_realData_sensitivity.py --findTrGeneNum -i /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/exp.isoforms.results -o /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/geneIsoNum.txt

(3) python analyze_realData_sensitivity.py --performancePlotExpress --oracleFa /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/oracle/snyder_oracle_et25_ct100_cut0.fa --numIsoFile /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/reference/snyder/geneIsoNum.txt --recLog /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/reconstructed_log.txt --cmpLog /data1/shunfu1/ref_shannon_modi/data/snyder_0729a/_res_stringtie_f_0_c_0.001/stringtie_log.txt --o /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/snyder_0807b/sensitivity_analysis.txt --expFile /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/exp_star/snyder/exp.isoforms.results --expFormat RSEM -L 101 -S 0 --isoLow 1 --isoHigh 100

'''

if __name__ == '__main__':

    args = sys.argv

    if '--findCoveredTr' in args:
        findCoveredTr(args)
    elif '--findTrGeneNum' in args:
        if '-i' in args:
            findTrGeneNum(args)
        elif '-f' in args:
            findTrGeneNum_f(args)
    elif '--performancePlotExpress' in args:
        performancePlotExpress(args)

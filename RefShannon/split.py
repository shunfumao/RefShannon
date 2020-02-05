from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import pysam
import sys, os, pdb
from RefShannon.util import *
from RefShannon.dep_path import extAssembler_paths, tool_paths

'''
hits.sam --> chr1/hits.sam, chr2/hits.sam ...

usage:

python split.py -O chrs_dir -i sam_file

'''

def parse_args_split():

    args = sys.argv

    args_error = False 

    if '-O' in args:
        chrs_dir = args[args.index('-O')+1]
    else:
        args_error = True

    if '-i' in args:
        sam_file = args[args.index('-i')+1]
    else:
        args_error = True

    return [args_error, sam_file, chrs_dir]

def split():

    [args_error, filename, loc]=parse_args_split()
    if args_error==True:
        print('split arguments error')
        return

    out_files = {}
    with open(filename) as f:
        for line in f:
            if line[0] == '@': continue
            fields = line.split()
            flags, genome = int(fields[1]), fields[2]
            if (flags >> 2) & 1:
                genome = 'chrUnmapped'
            if genome not in out_files:
                run_cmd('mkdir -p %s/%s/'%(loc, genome)) #os.mkdir(loc+genome)
                out_files[genome] = open("{}/hits.sam".format(loc+genome), 'w', 0)
            out_files[genome].write(line)

    for _, out_file in out_files.items():
        out_file.close()

    return

def split_bam(args):
  """
  input:
    input_bam:
    chroms: e.g. 'chr1,chr15' or 'all'
    output_sam: whether put a sam copy
    nJobs: num of CPUs
  output:
    [outDir]/[chrom]/[input_bam_fn].bam
  """
  input_bam, chroms, output_sam, nJobs, outDir = args

  input_bam_fn = os.path.basename(input_bam)
  # pdb.set_trace()

  out_bams = {} # key: chrom val: [path, file handler]
  if chroms == 'all':
    chroms_set = set()
    with pysam.AlignmentFile(input_bam, 'rb') as fh_input_bam:
      for line in fh_input_bam:
        chroms_set.add(fh_input_bam.get_reference_name(line.reference_id))
    chrom_list = list(chroms_set)
  else:
    chrom_list = [x for x in chroms.split(',') if x!='']
  # pdb.set_trace()

  with pysam.AlignmentFile(input_bam, 'rb') as fh_input_bam:
    for chrom in chrom_list:
      dir_path = '%s/%s/'%(outDir, chrom)
      cmd = 'mkdir -p %s'%(dir_path)
      run_cmd(cmd)

      file_path = '%s/%s'%(dir_path, input_bam_fn)
      fh = pysam.AlignmentFile(file_path, 'wb', header=fh_input_bam.header)
      out_bams[chrom] = [file_path, fh]
    # pdb.set_trace()

    cnt = 0
    for line in fh_input_bam:
      cnt += 1
      if cnt % 10000 == 0:
        print ('%d lines processed'%cnt)

      if line.reference_id==-1 or line.reference_start==-1: #invalid RNAME (e.g. chr name) or POS
        continue
      chrom = fh_input_bam.get_reference_name(line.reference_id)
      if chrom in out_bams:
        out_bams[chrom][1].write(line)

    for chrom in chrom_list:
      print('%s written'%out_bams[chrom][0])
      out_bams[chrom][1].close()
    # pdb.set_trace()

  if output_sam is True:
    for chrom in chrom_list:
      bam_path = out_bams[chrom][0]
      bam_dir = os.path.dirname(bam_path)
      bam_fn = os.path.basename(bam_path)[:-3] + 'sam'
      sam_path = '%s/%s'%(bam_dir, bam_fn)
      cmd = '%s view -o %s -O SAM -@ %d %s'%(
        tool_paths['samtools'],
        sam_path,
        nJobs,
        bam_path)
      run_cmd(cmd)
      print('%s written'%sam_path)

  return    

if __name__ == '__main__':
    split()

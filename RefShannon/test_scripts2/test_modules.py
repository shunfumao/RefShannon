from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

"""
test module functions under RefShannon/*.py
"""

import os

from RefShannon.util import parent_dir
from RefShannon.aligner import aligner_hisat2
from RefShannon.split import split_bam 

from RefShannon.test_scripts2.test_path import \
  path_test_aligner_hisat2, \
  path_test_split_bam

PATH = os.path.dirname(__file__)
ROOT = parent_dir(PATH)

def test_aligner_hisat2():
  for case in path_test_aligner_hisat2.keys():
    args = path_test_aligner_hisat2[case]
    aligner_hisat2(args)
  return

def test_split_bam():
  for case in path_test_split_bam.keys():
    case_args = path_test_split_bam[case]
    args = (
      case_args['input_bam'],
      case_args['chroms'],
      case_args['output_sam'],
      case_args['nJobs'],
      case_args['outDir']
      )
    split_bam(args)
  return

if __name__ == "__main__":
  test_aligner_hisat2()
  # test_split_bam()
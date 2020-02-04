from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

"""
test module functions under RefShannon/*.py
"""

import os

from RefShannon.util import parent_dir
from RefShannon.aligner import aligner_hisat2
from RefShannon.test_scripts2.test_path import path_test_aligner_hisat2

PATH = os.path.dirname(__file__)
ROOT = parent_dir(PATH)

def test_aligner_hisat2():
  for case in path_test_aligner_hisat2.keys():
    args = path_test_aligner_hisat2[case]
    aligner_hisat2(args)
  return

if __name__ == "__main__":
  test_aligner_hisat2()
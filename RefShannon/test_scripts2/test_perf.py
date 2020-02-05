from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

"""
Performance Evaluation
"""

from RefShannon.test_scripts2.roc import \
  exAssembler_roc, refShannon_roc

from RefShannon.test_scripts2.test_path import \
  path_test_exAssembler_roc, path_test_refShannon_roc

def test_relative_path():
  import pdb
  pdb.set_trace()
  return

def test_refShannon_roc():
  args = (
      path_test_refShannon_roc["rootResDir"],
      path_test_refShannon_roc["T"],
      path_test_refShannon_roc["sam_file"],
      path_test_refShannon_roc["genomeFile"],
      path_test_refShannon_roc["pairedStr"],
      path_test_refShannon_roc["chrom"],
      path_test_refShannon_roc["Tref"],
      path_test_refShannon_roc["nJobs"]
    )
  refShannon_roc(args)
  return

def test_exAssembler_roc():
  args = (
    path_test_exAssembler_roc["alignment"],
    path_test_exAssembler_roc["genomeFile"],
    path_test_exAssembler_roc["cases"],
    path_test_exAssembler_roc["resDir"],
    path_test_exAssembler_roc["reference"],
    path_test_exAssembler_roc["nJobs"],
    )
  exAssembler_roc(args)
  return

if __name__ == "__main__":
  # test_relative_path()

  # test_refShannon_roc()
  test_exAssembler_roc()
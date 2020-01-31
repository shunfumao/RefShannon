from RefShannon.test_scripts2.roc import \
  exAssembler_roc

from RefShannon.test_scripts2.test_path import \
  path_test_roc

def test_exAssembler_roc():
  args = (
    path_test_roc["alignment"],
    path_test_roc["genomeFile"],
    path_test_roc["cases"],
    path_test_roc["resDir"],
    path_test_roc["reference"],
    path_test_roc["nJobs"],
    )
  exAssembler_roc(args)
  return

if __name__ == "__main__":
  test_exAssembler_roc()
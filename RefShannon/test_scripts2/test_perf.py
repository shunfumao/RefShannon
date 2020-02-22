from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

"""
Performance Evaluation
"""

import sys, pdb

from RefShannon.test_scripts2.roc_sens import \
  exAssembler_run, gen_logs, gen_sens, \
  exAssembler_roc, refShannon_roc

from RefShannon.test_scripts2.test_path import \
  path_exAssembler_run, path_test_gen_logs, path_test_gen_sens, \
  path_test_exAssembler_roc, path_test_refShannon_roc

def test_relative_path():
  import pdb
  pdb.set_trace()
  return

def test_refShannon_roc():
  for case in path_test_refShannon_roc.keys():
    case_args = path_test_refShannon_roc[case]
    args = (
        case_args["rootResDir"],
        case_args["T"],
        case_args["sam_file"],
        case_args["genomeFile"],
        case_args["pairedStr"],
        case_args["chrom"],
        case_args["Tref"],
        case_args["nJobs"]
      )
    refShannon_roc(args)
  return

def test_exAssembler_run(sys_args):
  print(sys_args)
  # pdb.set_trace()

  if len(sys_args) == 1:
    example_key_list = path_exAssembler_run.keys()
  elif len(sys_args) == 2:
    example_key_list = [sys_args[1]]

  for example_key in example_key_list:
    example_dic = path_exAssembler_run[example_key]
    args = (
      example_dic['case'],
      example_dic['alignment'],
      example_dic['genomeFile'],
      example_dic['resDir'],
      example_dic['nJobs'])
    exAssembler_run(args) 
  return

def test_gen_logs():
  for case_key in path_test_gen_logs.keys():
    print(case_key)
    case = path_test_gen_logs[case_key]
    case_args = (
      case['Tref'],
      case['Trec'],
      case['resDir']
      )
    gen_logs(case_args)
  return

def test_gen_sens():
  for case_key in path_test_gen_sens.keys():
    case = path_test_gen_sens[case_key]
    case_args = (
      case['IM_list'],
      case['oracleFa'],
      case['numIsoFile'],
      case['recLog'],
      case['cmpLog'],
      case['expFile'],
      case['expFormat'],
      case['L'],
      case['S'],
      case['outFileStem'],
      case['ablist']
      )
    gen_sens(case_args)
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

  """
  sens of exAssembler
  """
  # test_exAssembler_run(sys.argv)
  test_gen_logs()
  # test_gen_sens()

  """
  roc
  """
  # test_refShannon_roc()

  # test_exAssembler_roc()
  # test_exAssembler_run(sys.argv) # [cgmemtime] python test_perf.py [example_key]
  # test_gen_logs() # sim roc when fasta is available
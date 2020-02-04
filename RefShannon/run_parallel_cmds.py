from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import multiprocessing
import os
from RefShannon.util import run_cmd
from RefShannon.global_values import *
import pdb

def run_cmds(cmds,noJobs=MAX_PARALLEL_PROCESS):

    if noJobs==1:
        #pdb.set_trace()
        for cmd in cmds:
            run_cmd(cmd)
    else:
        #cmds is a tuple of strings
        #noJobs is a an integer with the no of Jobs to run in parallel
        p = multiprocessing.Pool(noJobs)
        p.map(run_cmd,cmds)
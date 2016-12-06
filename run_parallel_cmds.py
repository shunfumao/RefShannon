import multiprocessing
import os
from util import run_cmd
from global_values import *

def run_cmds(cmds,noJobs=MAX_PARALLEL_PROCESS):
    #cmds is a tuple of strings
    #noJobs is a an integer with the no of Jobs to run in parallel
    p = multiprocessing.Pool(noJobs)
    p.map(run_cmd,cmds)
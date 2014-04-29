#!/usr/bin/python2.6 -tt
# =======================================================
import os, sys, math

import run_event_tools
from multiprocessing import Pool
 
from ROOT import TFile, TChain, TFileCollection
from ROOT import gROOT, gSystem, gStyle
 
gROOT.Reset()
 
# =======================================================


def main():

  n = len(sys.argv)-1

  pool = Pool(n)

  inputLists = sys.argv[1:]
  triggerBits = [0,1,2]
  entries = -1
  
  processes = []
  for inputList in inputLists:
    processes.append(pool.apply_async(run_event_tools.EventTable,[inputList, triggerBits, entries]))
    
  for process in processes:
    process.get()
 
# _______________________________________________________  

  
if __name__ == '__main__':
  main()




# =======================================================

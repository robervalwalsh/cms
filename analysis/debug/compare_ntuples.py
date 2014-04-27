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

  pool = Pool(2)

  input1 = sys.argv[1]
  input2 = sys.argv[2]
  triggerBit = int(sys.argv[3])
  entries = -1
  if len(sys.argv) > 4:
    entries = int(sys.argv[4])
    
  output1 = "events_bit" + str(triggerBit) + input1
  output2 = "events_bit" + str(triggerBit) + input2
  
  process1 = pool.apply_async(run_event_tools.EventTable,[input1, output1, triggerBit, entries])
  process2 = pool.apply_async(run_event_tools.EventTable,[input2, output2, triggerBit, entries])
  
  eventTable1 = process1.get()
  eventTable2 = process2.get()
  
  run_event_tools.EventTableSummary(eventTable1,eventTable2)
 
# _______________________________________________________  

  
if __name__ == '__main__':
  main()




# =======================================================

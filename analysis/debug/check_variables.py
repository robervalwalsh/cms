#!/usr/bin/python2.6 -tt
# =======================================================
import os, sys, math

from ROOT import TFile, TChain, TFileCollection
from ROOT import gROOT, gSystem, gStyle
 
gROOT.Reset()
 
# =======================================================


def main():

  fileListFile = sys.argv[1]
  
  ## Open files with ntuples
  fileList = TFileCollection("dum","",fileListFile)
  chain = TChain("hbbanalysis/HBBTo4B")
  chain.AddFileInfoList(fileList.GetList())
  nentries = chain.GetEntries()

  print "Processing", nentries, "events"
  
  counter = 0
  
  for i in xrange(nentries):
    if i % 100000 == 0 and i > 0 : print i,"events processed"
    chain.GetEntry(i)
    if chain.trgAccept & (1<<0): counter += 1
    
  print "Total number of triggered events =", counter
    

 
# _______________________________________________________  

  
if __name__ == '__main__':
  main()




# =======================================================

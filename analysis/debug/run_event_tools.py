#!/usr/bin/python2.6 -tt
# =======================================================

import os, sys, math
import json

from ROOT import TFile, TChain, TFileCollection
from ROOT import gROOT, gSystem, gStyle
 
gROOT.Reset()
 
# =======================================================

# Given an list of root input files, the name of an output file
# a table consisting of run : [list of events] is created and
# dumped in the output file

def EventTable(fileListFile, outputFile, triggerBit, entries):

  ## Open files with ntuples
  fileList = TFileCollection("dum","",fileListFile)
  chain = TChain("hbbanalysis/HBBTo4B")
  chain.AddFileInfoList(fileList.GetList())
  nentries = entries
  if entries < 0 :
    nentries = chain.GetEntries()

  eventTable = {}
  
#  print "Processing", nentries, "events"
  
  for i in xrange(nentries):
  
#    if i % 100000 == 0 and i > 0 : print i,"events processed"
    
    chain.GetEntry(i)

    # check if trigger trigger was accepted. todo: pass the trigger bit as parameter
    isMyTrigger = False
    if chain.trgAccept & (1<<triggerBit): isMyTrigger = True
    
    if not isMyTrigger: continue 
    
    run = chain.Run
    event = chain.Event
    
    if not eventTable.get(run):
      eventTable[run] = [event]
    else:
#      if eventTable[run].count(event) == 0: # avoid duplicates for the moment
        eventTable[run] += [event]
    
  # end nentries loop
  if entries > 0 :
    OutputEventTable(eventTable,outputFile)
  
#  print eventTable
  
  return eventTable
  
# =======================================================

# Given the run-events table, save it to a file in the format
# run1 : event1, event2, ...
# run2 : event1, event2, ...
# ...

def OutputEventTable(eventTable, filename):
  f = open(filename, "w")
  
  for run in sorted(eventTable.keys()):
    f.write("%i : " % run)
    events = sorted(eventTable[run])
    formattedEvents = ""
    for event in events[:-1]:
      formattedEvents += str(event) + ", "
    formattedEvents += str(events[-1])
    f.write(formattedEvents+"\n")
    
  f.close()
  
# =======================================================

# Given two run-events table, check that the run lists are the same
# that there is no duplicated events, and search for partial intersects 
def EventTableSummary(eventTable1,eventTable2):
  total1 = 0
  total2 = 0
  totaldiff = 0
  
  runs1 = sorted(eventTable1.keys())
  runs2 = sorted(eventTable2.keys())
  runsDiff = list(set(runs1) - set(runs2)) + list(set(runs2) - set(runs1))
  if ( len(runsDiff) != 0 ):
    print "* warning * -- List of runs in table 1 and table 2 differ!!!"
    print "               The list of different runs is"
    print sorted(runsDiff)
    print "---------------------------------"
  else:
    print " Both tables have the same runs. It is a good start!"
    print "---------------------------------"
  
  
  for run in sorted(eventTable1.keys()):
    events1 = sorted(eventTable1[run])
    if len(events1) != len(set(events1)):
      print "* warning * -- In run", run, ",table1, there seems to be duplicated events!!!"
      print "---------------------------------"
    total1 += len(set(events1))
    if eventTable2.get(run):
      events2 = sorted(eventTable2[run])
      if len(events2) != len(set(events2)):
        print "* warning * -- In run", run, ", table2, there seems to be duplicated events!!!"
        print "---------------------------------"
      total2 += len(set(events2))
    else:
      print "* warning * -- Run", run, "not in table 2"
      print "---------------------------------"
      continue
    diff = list(set(events1)-set(events2)) + list(set(events2)-set(events1))
    print "---------------------------------"
    print "Run = ", run, " : number of events that differ between table 1 and table 2 =  ", len(diff)
    totaldiff += len(diff)
    
  print "---------------------------------"
  print "---------------------------------"
  print "Total events in table 1 = ", total1
  print "Total events in table 2 = ", total2
  print "Total different events  = ", totaldiff
    

# =======================================================

# Given two fiels with run-events table each, create the appropriate
# dictionaries to compare 

def MakeEventTable(filename):
  f = open(filename,"r")
  
  table = {}
  lines = f.readlines()
  for line in lines:
    run = int(line.split(":")[0])
    events = line.strip(str(run) +" : ").split(",") # list of events removing the run
    events = map (int , events)
    table[run] = sorted(events)
    
  return table


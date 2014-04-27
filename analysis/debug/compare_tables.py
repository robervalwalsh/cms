#!/usr/bin/python2.6 -tt
# =======================================================
import os, sys, math

import run_event_tools
 
# =======================================================


def main():
  
  input1 = sys.argv[1]
  input2 = sys.argv[2]
  
  eventTable1 = run_event_tools.MakeEventTable(input1)
  eventTable2 = run_event_tools.MakeEventTable(input2)
  
  run_event_tools.EventTableSummary(eventTable1,eventTable2)
 
    
# _______________________________________________________  

  
if __name__ == '__main__':
  main()




# =======================================================

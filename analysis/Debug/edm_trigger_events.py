#!/usr/bin/python2.6 -tt

import os
import json
import subprocess

# https://docs.python.org/2/library/optparse.html
from   optparse import OptionParser
import das_client

import ROOT
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle('Plain') # white background

# ___________________________________________________________________________________


class MyOptionParser:
    """
    My option parser
    """
    def __init__(self):
        usage  = "Usage: %prog [options]\n"
        usage += "For more help..."
        self.parser = OptionParser(usage=usage)
        self.parser.add_option("--run", action="store", type="int", default=0,
                               dest="run", help="specify the run")
        self.parser.add_option("--lumi", action="store", type="int", default=0,
                               dest="lumi", help="specify lumi section")
        self.parser.add_option("--trigger", action="store", type="string", default="HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV",
                               dest="trigger", help="specify the trigger")
        self.parser.add_option("--dataset", action="store", type="string", default="/MultiJet/Run2012A-22Jan2013-v1/AOD",
                               dest="dataset", help="specify the dataset")
    
    def get_opt(self):
        """
        Returns parse list of options
        """
        return self.parser.parse_args()

# ___________________________________________________________________________________

def get_list_of_files( run, lumi, dataset):
   
   query = 'file dataset=' + dataset
   if ( run > 0 ):  query += ' run='  + str(run)
   if ( lumi > 0 ): query += ' lumi=' + str(lumi)
   
   # Using DAS client to retrieve the filenames
   das_client_command = "das_client.py --limit=0 --query='" + query +"'"
   
   print "Running the command ", das_client_command, "..."
   filenames = os.popen(das_client_command).read().split('\n')
   
   return filenames[:-1] # the last entry is empty!?

# ___________________________________________________________________________________

def get_list_of_events(events,run,lumi,mytrigger):
   # Trigger information
   ''' TriggerResults is deprecated, but it is not clear how to use the ConfigProvider in Python'''
   handleTrigger  = Handle ('edm::TriggerResults')
   labelTrigger = ('TriggerResults','','HLT')
   
   # Loop over the events
   for event in events:
      evento = event.object()                                     # This is now an fwlite::Event object! Provides more information than python::Events.
      event.getByLabel (labelTrigger, handleTrigger)              # From the python::Events class instead of fwlite::Event; more straightforward
      if handleTrigger.isValid():
         triggerResults = handleTrigger.product()                 # Access the objects in the Handle
         triggerNames = evento.triggerNames(triggerResults)
         triggers = triggerNames.triggerNames()                   # List of triggers
         for trigger in triggers:
            if mytrigger in trigger:
               triggerIndex = triggerNames.triggerIndex(trigger)
               triggerAccepted = triggerResults.accept(triggerIndex)
               if triggerAccepted:
                  print evento.eventAuxiliary().run(), evento.eventAuxiliary().luminosityBlock(),evento.eventAuxiliary().event()
                     
   return 0


# ___________________________________________________________________________________

def main():
   
   # options from command line
   optmgr  = MyOptionParser()
   options, _ = optmgr.get_opt()  # tuple,  unpack the return value and assign it the variable named to the left; _ single value to unpack
   
   run       = options.run
   lumi      = options.lumi
   dataset   = options.dataset
   mytrigger = options.trigger
   
   filenames = get_list_of_files(run,lumi,dataset)
   if len(filenames) < 1:
      print "No file found!"
      return -1
   
   print "Running on", len(filenames), "files..."
   
   for filename in filenames:
      fullname = 'dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2' + filename
      print "   Reading file", filename, "..."
   
      ''' https://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_3_9/doc/html/d4/dab/classpython_1_1Events.html
      This is NOT a collection of fwlite::Event objects.
      The fwlite::Event object is obtained using the object() method (see below) '''
      events = Events(fullname)
      
      eventlist = get_list_of_events(events,run,lumi,mytrigger)


# ___________________________________________________________________________________


if __name__ == '__main__':
   main()


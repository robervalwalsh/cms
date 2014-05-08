#!/usr/bin/python2.6 -tt

import ROOT
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle('Plain') # white background


def main():
   
   ''' https://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_3_9/doc/html/d4/dab/classpython_1_1Events.html
   This is NOT a collection of fwlite::Event objects.
   The fwlite::Event object is obtained using the object() method (see below) '''
   events = Events('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2012D/BJetPlusX/AOD/22Jan2013-v1/10000/0016172B-A992-E211-8930-90B11C18E296.root')

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
            if "HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV" in trigger:
               triggerIndex = triggerNames.triggerIndex(trigger)
               triggerAccepted = triggerResults.accept(triggerIndex)
               if triggerAccepted:
                  print evento.eventAuxiliary().run(), evento.eventAuxiliary().luminosityBlock(),evento.eventAuxiliary().event()


''' ___________________________________________________________________________________ '''


if __name__ == '__main__':
   main()


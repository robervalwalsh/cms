// -*- C++ -*-
//
// Package:    KpipiTrueFilter
// Class:      KpipiTrueFilter
// 
/**\class KpipiTrueFilter KpipiTrueFilter.cc Analysis/CharmHadrons/src/KpipiTrueFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Roberval Walsh,01b/352,4849,
//         Created:  Fri Feb 12 18:13:21 CET 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


using namespace edm;
using namespace reco;

//
// class declaration
//

class KpipiTrueFilter : public edm::EDFilter {
   public:
      explicit KpipiTrueFilter(const edm::ParameterSet&);
      ~KpipiTrueFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
   InputTag src_;
      
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KpipiTrueFilter::KpipiTrueFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   src_      = iConfig.getParameter<InputTag> ("source");
}


KpipiTrueFilter::~KpipiTrueFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
KpipiTrueFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   Handle<View<Candidate> > genHandler;
   iEvent.getByLabel( src_ , genHandler );
   
   for ( View<Candidate>::const_iterator mcp = genHandler->begin(); mcp != genHandler->end(); ++mcp )
   {
      
      if ( abs(mcp -> pdgId()) != 411 && mcp -> numberOfDaughters() != 3 ) continue;
      int countK = 0;
      int countPi = 0;
      for ( Candidate::const_iterator dau = mcp -> begin(); dau != mcp -> end(); ++dau )
      {
          if ( abs(dau -> pdgId()) == 321 ) ++countK;
          if ( abs(dau -> pdgId()) == 211 ) ++countPi;
      }
      if ( countK == 1 && countPi == 2 )
      {
          return true;
      }
      
   }
   
   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
KpipiTrueFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KpipiTrueFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(KpipiTrueFilter);

// -*- C++ -*-
//
// Package:    KpipiGenAnalyzer
// Class:      KpipiGenAnalyzer
// 
/**\class KpipiGenAnalyzer KpipiGenAnalyzer.cc Analysis/CharmHadrons/src/KpipiGenAnalyzer.cc
 
 Description: <one line class summary>
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Roberval Walsh,01b/352,4849,
//         Created:  Tue Feb  9 10:07:31 CET 2010
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Common/interface/View.h"

#include "TROOT.h"
#include "TH1D.h"

using namespace edm;
using namespace reco;

//
// class declaration

class KpipiGenAnalyzer : public edm::EDAnalyzer
{
public:
   explicit KpipiGenAnalyzer(const edm::ParameterSet&);
   ~KpipiGenAnalyzer();
   
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   // ----------member data ---------------------------
   
   unsigned int nGenDplus;
   
   // Config parameters
   InputTag srcGen_;
   
   double ptMin_;
   std::vector<double> etaRange_;
   
   // Get the file service
   Service<TFileService> fs_;
   
   // Histogram handler
   std::map<std::string, TH1D *> Histogram_;
   //
   
};

KpipiGenAnalyzer::KpipiGenAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcGen_      = iConfig.getParameter<InputTag> ("GeneratedParticles");
   etaRange_    = iConfig.getParameter< std::vector<double> > ("EtaRange");
   ptMin_  = iConfig.getParameter<double> ("PtMin");
   
   // Histograms
   //   // D+
   Histogram_["Ncands"] = fs_ -> make<TH1D>( "Ncands", "", 10, 0., 10. );
   Histogram_["Mass"]   = fs_ -> make<TH1D>( "Mass", "", 40, 1.7, 2.1 );
   Histogram_["Pt"]     = fs_ -> make<TH1D>( "Pt", "", 300, 0., 30. );
   Histogram_["Eta"]    = fs_ -> make<TH1D>( "Eta", "", 200, -10, 10. );
   //   
   
   nGenDplus = 0;
}

KpipiGenAnalyzer::~KpipiGenAnalyzer()
{
   
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}

// ------------ method called to for each event  ------------
void KpipiGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Generated particles
   Handle<View<GenParticle> > genParticles;
   iEvent.getByLabel (srcGen_, genParticles);
   
   unsigned int counter = 0;
   
    for (View<GenParticle>::const_iterator mcp = genParticles->begin(); mcp != genParticles->end(); ++mcp)
   {
      if ( abs(mcp -> pdgId()) != 411 && mcp -> numberOfDaughters() != 3 ) continue;
      int countK = 0;
      int countPi = 0;
      for ( Candidate::const_iterator dau = mcp -> begin(); dau != mcp -> end(); ++dau )
      {
         if ( abs(dau -> pdgId()) == 321 ) ++countK;
         if ( abs(dau -> pdgId()) == 211 ) ++countPi;
      }
      if ( countK != 1 || countPi != 2 ) continue;
      Histogram_["Pt"] -> Fill(mcp -> pt());
      Histogram_["Eta"] -> Fill(mcp -> eta());
      Histogram_["Mass"] -> Fill(mcp -> mass());
      
      if ( mcp -> pt() < ptMin_ ) continue;
      if ( mcp -> eta() < etaRange_[0] || mcp -> eta() > etaRange_[1] ) continue;
      ++counter;
      ++nGenDplus;
      
   }
   if ( counter > 0 ) Histogram_["Ncands"] -> Fill(counter);
}


// ------------ method called once each job just before starting event loop  ------------
void KpipiGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpipiGenAnalyzer::endJob()
{
   std::cout << std::endl;
   std::cout << "=============================================" << std::endl;
   std::cout << "    *** endJob ***" <<  std::endl;
   
   std::cout << "Achei " << nGenDplus << " D+-" << std::endl;
   
   std::cout << std::endl;
   
}


// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(KpipiGenAnalyzer);

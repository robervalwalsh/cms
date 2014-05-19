// -*- C++ -*-
//
// Package:    DplusLeptonTrueAnalyzer
// Class:      DplusLeptonTrueAnalyzer
// 
/**\class DplusLeptonTrueAnalyzer DplusLeptonTrueAnalyzer.cc Analysis/CharmHadrons/src/DplusLeptonTrueAnalyzer.cc
 
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

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/CompositeRefCandidateT.h"

#include "DataFormats/Common/interface/View.h"

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace edm;
using namespace reco;

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > CovarianceMatrix;
typedef std::vector<CovarianceMatrix> CovarianceMatrixCollection;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;
typedef std::vector<Point> PointCollection;
typedef ROOT::Math::SVector<double, 3> Vector3D;


//
// class declaration

class DplusLeptonTrueAnalyzer : public edm::EDAnalyzer
{
public:
   explicit DplusLeptonTrueAnalyzer(const edm::ParameterSet&);
   ~DplusLeptonTrueAnalyzer();
   
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   const Candidate * Bparent ( const Candidate * );
   void FillHistograms( const Candidate * , const Candidate * );
   void CreateHistograms();
   
      
   // ----------member data ---------------------------
   
   // Config parameters
   InputTag srcGen_;
   double PtMin_;
   std::vector<double> EtaRange_;
   
   // Get the file service
   Service<TFileService> fs_;
   
   // Histogram handler
   std::map<std::string, TH1D *> Histogram_;
   
   
};

DplusLeptonTrueAnalyzer::DplusLeptonTrueAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcGen_              = iConfig.getParameter<InputTag> ("GeneratedParticles");
   PtMin_               = iConfig.getParameter<double> ("PtMin");
   EtaRange_            = iConfig.getParameter< std::vector<double> > ("EtaRange");

   //
   CreateHistograms();
   

}

DplusLeptonTrueAnalyzer::~DplusLeptonTrueAnalyzer()
{
   
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}

// ------------ method called to for each event  ------------
void DplusLeptonTrueAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // Generated Dpluses
   // Generated particles
   Handle<View<GenParticle> > genParticles;
   iEvent.getByLabel (srcGen_, genParticles);
   
   for (View<GenParticle>::const_iterator mcp = genParticles->begin(); mcp != genParticles->end(); ++mcp)
   {
      const Candidate * trueDplus = 0;
      if ( abs(mcp -> pdgId()) != 411 ) continue;
      if ( mcp -> pt()  < PtMin_ ) continue;
      if ( mcp -> eta() < EtaRange_[0] || mcp -> eta() > EtaRange_[1] ) continue;
      
      trueDplus = &*(mcp);
      Histogram_["DplusPt"]     ->  Fill(trueDplus -> pt());
      Histogram_["DplusEta"]    ->  Fill(trueDplus -> eta());

      const Candidate * mother = Bparent(trueDplus);
      if ( mother != 0 )
      {
         size_t nSisters = mother -> numberOfDaughters();
         for ( size_t i = 0; i < nSisters; ++i )
         {
            const Candidate * sis = mcp -> daughter(i);
            if ( sis == 0 ) continue;
            if ( abs(sis -> pdgId()) == 11 || abs(sis -> pdgId()) == 13 ) FillHistograms( trueDplus, sis );
         }
      }
   }
   
}



// ------------ method called once each job just before starting event loop  ------------
void DplusLeptonTrueAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DplusLeptonTrueAnalyzer::endJob()
{
}

const Candidate * DplusLeptonTrueAnalyzer::Bparent ( const Candidate * mcp )
{
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      char pdg[8];
      const char * b = "5";
      sprintf(pdg,"%d",abs(mom -> pdgId()));
      if ( pdg[0] == *b  ) return mom;
   }
   
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      Bparent(mom); 
   }
   
   return 0;
   
}


void DplusLeptonTrueAnalyzer::FillHistograms( const Candidate * dplus, const Candidate * lepton )
{
   Histogram_["DplusPtB"]     ->  Fill(dplus -> pt());
   Histogram_["DplusEtaB"]    ->  Fill(dplus -> eta());
   
   if ( abs(lepton -> pdgId()) == 11 ) 
   {
      Histogram_["ElectronPt"]     ->  Fill(lepton -> pt());
      Histogram_["ElectronEta"]    ->  Fill(lepton -> eta());
   }
   
   if ( abs(lepton -> pdgId()) == 13 ) 
   {
      Histogram_["MuonPt"]     ->  Fill(lepton -> pt());
      Histogram_["MuonEta"]    ->  Fill(lepton -> eta());
   }
   
   Histogram_["LeptonVsDplusPtRatio"] -> Fill(lepton -> pt()/dplus -> pt());
   
   GlobalVector dplusVector(dplus -> px(),dplus -> py(),dplus -> pz());
   GlobalVector leptonVector(lepton -> px(),lepton -> py(),lepton -> pz());
   double angle = acos(fabs(dplusVector.dot(leptonVector)/(dplusVector.mag()*leptonVector.mag())));
   Histogram_["LeptonVsDplusAngle"] -> Fill(angle);
}

void DplusLeptonTrueAnalyzer::CreateHistograms()
{
   // Histograms
   Histogram_["DplusPt"]      = fs_ -> make<TH1D>( "DPlusPt", "", 200, 0., 20. );
   Histogram_["DplusEta"]     = fs_ -> make<TH1D>( "DplusEta", "", 200, -3.5, 3.5 );
   Histogram_["DplusPtB"]     = fs_ -> make<TH1D>( "DPlusPtB", "", 200, 0., 20. );
   Histogram_["DplusEtaB"]    = fs_ -> make<TH1D>( "DplusEtaB", "", 200, -3.5, 3.5 );
   //   
   Histogram_["ElectronPt"]     = fs_ -> make<TH1D>( "ElectronPt", "", 200, 0., 5. );
   Histogram_["ElectronEta"]    = fs_ -> make<TH1D>( "ElectronEta", "", 200, -3.5, 3.5 );
   //
   Histogram_["MuonPt"]     = fs_ -> make<TH1D>( "MuonPt", "", 200, 0., 5. );
   Histogram_["MuonEta"]    = fs_ -> make<TH1D>( "MuonEta", "", 200, -3.5, 3.5 );
   
   Histogram_["LeptonVsDplusPtRatio"] = fs_ -> make<TH1D>( "LeptonVsDplusPtRatio", "", 200, 0., 2. );
   Histogram_["LeptonVsDplusAngle"]   = fs_ -> make<TH1D>( "LeptonVsDplusAngle", "", 3200, 0., 3.2 );
}

// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(DplusLeptonTrueAnalyzer);

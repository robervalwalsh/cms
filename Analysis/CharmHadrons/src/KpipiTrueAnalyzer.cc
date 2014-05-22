// -*- C++ -*-
//
// Package:    KpipiTrueAnalyzer
// Class:      KpipiTrueAnalyzer
// 
/**\class KpipiTrueAnalyzer KpipiTrueAnalyzer.cc Analysis/CharmHadrons/src/KpipiTrueAnalyzer.cc
 
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

class KpipiTrueAnalyzer : public edm::EDAnalyzer
{
public:
   explicit KpipiTrueAnalyzer(const edm::ParameterSet&);
   ~KpipiTrueAnalyzer();
   
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   double DistanceOfClosestApproach( const TrackRef & , const TrackRef &, const MagneticField * );
   bool IsGoodCandidate( const pat::PATObject<reco::VertexCompositeCandidate> & , const Vertex &, const MagneticField * );
   
   bool CheckParentage ( const GenParticle *, int pdg );
   bool CheckParentage ( const GenParticleRef &, int pdg );
   bool CheckParentage ( const Candidate *, int pdg );
   bool AreSiblings ( const Candidate * , const GenParticleRef & , const GenParticleRef & , const GenParticleRef & );
   
   VertexCompositeCandidate TrueKpipi ( const GenParticleRef & , const GenParticleRef & , const GenParticleRef & , const Candidate *);
   Vertex TruePrimaryVertex ( const Candidate * , bool & );
   
   bool FromBdecay ( const Candidate *  );
   bool FromCdecay ( const Candidate *  );
   void FillHistograms( const pat::PATObject<reco::VertexCompositeCandidate> &, const VertexCompositeCandidate & , const Vertex &, const Vertex &, const MagneticField * );
   
   void PrintPtEtaPhi ( const Candidate * );
   
   void CreateHistograms();
   void SetDummyCuts();
   void SetCuts(const edm::ParameterSet&);
   
   
      
   // ----------member data ---------------------------
   
   // Config parameters
   InputTag srcKpipi_;
   InputTag srcPVs_;
   InputTag srcGen_;
   InputTag srcPocas_;
   InputTag srcPocaErrors_;
   
   // Track selection
   double trkDxyMax_;
   double trkDzMax_;
   int    trkPixelHitsMin_;
   int    trkSiliconHitsMin_;
   double trkNormChi2Max_;
   std::vector<double> trkEtaRange_;
   double trkKaonPtMin_;
   double trkPionPtMin_;
   double trkDCAMax_;
   
   // D+ candidates selection
   double KpipiPtMin_;
   std::vector<double> KpipiEtaRange_;
   std::vector<double> KpipiMassWindow_;
   double KpipiDecayVertexProbMin_;
   double KpipiPVAngleMax_;
   
   double KpipiDecayLengthSigMin_;
   double KpipiDecayLengthMin_;
   double KpipiDecayLengthSig2DMin_;
   double KpipiDecayLength2DMin_;
   
   double KpipiProperDecayLengthMin_;
   
   double KpipiCosThStarKMin_;
   
   unsigned int quarkMother_;
   
   bool useCuts_;
   
   unsigned int counter;
   
   unsigned int nGenDplus;
   
   // Get the file service
   Service<TFileService> fs_;
   
   // Histogram handler
   std::map<std::string, TH1D *> Histogram_;
   std::map<std::string, TH2D *> Histogram2_;
   
   TTree * Tree_;
   double pt_rec_;
   double m_rec_;
   double lxy_rec_;
   double lxyz_rec_;
   double ct_rec_;
   
   double pt_true_;
   double m_true_;
   double lxy_true_;
   double lxyz_true_;
   double ct_true_;
   
   
   
   //
   
};

KpipiTrueAnalyzer::KpipiTrueAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcKpipi_            = iConfig.getParameter<InputTag> ("KpipiCandidates");
   srcPVs_              = iConfig.getParameter<InputTag> ("PrimaryVertex");
   srcPocas_            = iConfig.getParameter<InputTag> ("PointOfClosestApproach");
   srcPocaErrors_       = iConfig.getParameter<InputTag> ("PointOfClosestApproachError");
   srcGen_              = iConfig.getParameter<InputTag> ("GeneratedParticles");
   
   quarkMother_         = iConfig.getParameter<unsigned int> ("QuarkMother");
   
   //
   useCuts_ = iConfig.getParameter<bool> ("UseCuts");
   
   counter = 0;
   
   SetCuts( iConfig );

   CreateHistograms();

   
   nGenDplus = 0;
   
}

KpipiTrueAnalyzer::~KpipiTrueAnalyzer()
{
   
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}

// ------------ method called to for each event  ------------
void KpipiTrueAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   GlobalVector nullVector(0,0,0);
   
   ESHandle<MagneticField> bFieldHandler;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandler);
   
   // primary vertex
   Handle<VertexCollection> pvHandler;
   iEvent.getByLabel(srcPVs_, pvHandler);
   Vertex primaryVertex = pvHandler -> at(0);
   
   // Kpipi
   Handle<View<VertexCompositeCandidate> > kpipiHandler;
   iEvent.getByLabel (srcKpipi_, kpipiHandler);
   
   // Track candidates
   Handle<View<RecoChargedCandidate> > trackCandHandler;
   iEvent.getByLabel( "allTracks" , trackCandHandler );
   
   // Match RECO-TRUE (edm::Association)
   Handle<GenParticleMatch> trackMatchHandler;
   iEvent.getByLabel( "generalTracksGenParticlesMatch" , trackMatchHandler );
   
   unsigned int nMatchedDplus = 0;
   for ( unsigned int i = 0 ; i < kpipiHandler -> size() ; ++i )
   {
      // k-pi-pi candidate
      RefToBase<VertexCompositeCandidate> kpipiRef = kpipiHandler -> refAt(i);
      pat::PATObject<VertexCompositeCandidate> kpipi(kpipiRef);
      if ( useCuts_ )
         if ( ! IsGoodCandidate(kpipi, primaryVertex, &(*bFieldHandler))) continue;
      
      // track match
      int countTrackMatches = 0;
      GenParticleRef trueKaonRef;
      GenParticleRef truePi1Ref;
      GenParticleRef truePi2Ref;
      
      for ( unsigned int k = 0; k < 3; ++k )  // loop over dplus tracks
      {
         RecoChargedCandidate dplusDaughter = *(dynamic_cast<const RecoChargedCandidate *> (kpipi.daughter(k))) ;
         TrackRef dplusTrackRef = dplusDaughter.track();
         for ( unsigned int j = 0; j < trackCandHandler -> size(); ++j ) // loop over all tracks
         {
            RefToBase<RecoChargedCandidate> trackCandRef = trackCandHandler -> refAt(j);
            pat::PATObject<reco::RecoChargedCandidate> patTrackCand(trackCandRef);
            TrackRef trackRef = patTrackCand.track();
            if ( dplusTrackRef != trackRef ) continue;
            // get the match of the track
            GenParticleRef trueTrackRef = (*trackMatchHandler)[trackCandRef];
            patTrackCand.addGenParticleRef(trueTrackRef);
            if ( ! CheckParentage(trueTrackRef,kpipi.pdgId()) ) continue;
            
            if ( k == 0 ) trueKaonRef = trueTrackRef;
            if ( k == 1 ) truePi1Ref  = trueTrackRef;
            if ( k == 2 ) truePi2Ref  = trueTrackRef;
            
            ++countTrackMatches;
            
         }  // loop over all tracks
      }  // loop over dplus tracks
      if ( countTrackMatches > 3 ) std::cout << " *** KpipiTrueAnalyzer: More than 3 tracks match ***" << std::endl; 
      if ( countTrackMatches != 3 ) continue;
      
      // See if the daughters of a true D+ match the matched tracks
      const Candidate * trueDplus = 0;
      Handle<View<Candidate> > genHandler;
      iEvent.getByLabel( srcGen_ , genHandler );
      for ( View<Candidate>::const_iterator mcp = genHandler->begin(); mcp != genHandler->end(); ++mcp )
      {
         if ( mcp -> pdgId() != kpipi.pdgId() ) continue;
         size_t ndau =  mcp -> numberOfDaughters();
         if ( ndau != 3 ) continue;
         
         int cntKpipi = 0;
         for ( Candidate::const_iterator dau = mcp->begin(); dau != mcp->end(); ++dau )
            if ( (dau -> pdgId() == -(int)mcp -> charge()*321 || dau -> pdgId() ==  (int)mcp -> charge()*211) && dau -> status() == 1 ) ++cntKpipi;
         if ( cntKpipi != 3 ) continue;
         
         if ( ! AreSiblings( &*(mcp), trueKaonRef, truePi1Ref, truePi2Ref ) ) continue;
         
         trueDplus = &*(mcp);
         break;
      }
      
      if ( trueDplus == 0 ) continue;
      
      ++nMatchedDplus;
      
      if ( quarkMother_ == 4 || quarkMother_ == 5 )
      {
         if ( quarkMother_ == 4 && ! FromCdecay(trueDplus) ) continue;
         if ( quarkMother_ == 5 && ! FromBdecay(trueDplus) ) continue;
      }
      
      // true primary vertex
      bool foundVertex;
      Vertex truePrimaryVertex = TruePrimaryVertex(trueDplus, foundVertex);
      if ( ! foundVertex ) std::cout << " *** KpipiTrueAnalyzer::analyze *** True primary vertex not found" << std::endl;
      
      // Combine the true tracks
      VertexCompositeCandidate trueKpipi = TrueKpipi( trueKaonRef, truePi1Ref, truePi2Ref, trueDplus );
      
      // Fill the Histograms
      FillHistograms(kpipi, trueKpipi, primaryVertex, truePrimaryVertex,  &(*bFieldHandler));
      Tree_ -> Fill();
      
   } // kpipi candidate loop (i)
   Histogram_["NDplus"] -> Fill(nMatchedDplus);
   
// Generated Dpluses
   // Generated particles
   Handle<View<GenParticle> > genParticles;
   iEvent.getByLabel (srcGen_, genParticles);
   
   unsigned int countD = 0;
   
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
      Histogram_["PtGen"] -> Fill(mcp -> pt());
      Histogram_["EtaGen"] -> Fill(mcp -> eta());
      Histogram_["MassGen"] -> Fill(mcp -> mass());
      
      if ( mcp -> pt() < KpipiPtMin_ ) continue;
      if ( mcp -> eta() < KpipiEtaRange_[0] || mcp -> eta() > KpipiEtaRange_[1] ) continue;
      ++countD;
      ++nGenDplus;
      
   }
   if ( countD > 0 ) Histogram_["NGen"] -> Fill(countD);
   
}



// ------------ method called once each job just before starting event loop  ------------
void KpipiTrueAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpipiTrueAnalyzer::endJob()
{
   std::cout << std::endl;
   std::cout << "=============================================" << std::endl;
   std::cout << "    *** endJob ***" <<  std::endl;
   
//   std::cout << "Achei " << counter << " D+-" << std::endl;
   
   std::cout << std::endl;
   
}


double KpipiTrueAnalyzer::DistanceOfClosestApproach( const TrackRef & tRef1, const TrackRef & tRef2, const MagneticField * bField )
{
   double dca;
   double dummy = -999.;
   
   TransientTrack t1 (tRef1, &(*bField) );
   TransientTrack* tt1  = &t1;
   if ( ! tt1 -> impactPointTSCP().isValid() ) return dummy;
   TransientTrack t2 (tRef2, &(*bField) );
   TransientTrack* tt2  = &t2;
   if ( ! tt2 -> impactPointTSCP().isValid() ) return dummy;
   
   FreeTrajectoryState tt1State  = tt1 -> impactPointTSCP().theState();
   FreeTrajectoryState tt2State  = tt2  -> impactPointTSCP().theState();
   // Closest approach and crossing points of the tracks
   ClosestApproachInRPhi cApp;
   GlobalPoint cxPt;
   // track1- track2
   cApp.calculate(tt1State, tt2State);
   if ( ! cApp.status() ) return dummy;
   dca = fabs( cApp.distance() );
   if ( dca < 0. || dca > trkDCAMax_ ) return dummy;
   cxPt = cApp.crossingPoint();
   if ( cxPt.transverse() > 120. || fabs(cxPt.z()) > 300.) return dummy;
   TrajectoryStateClosestToPoint tt1TSCP = tt1 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt1TSCP.isValid() ) return dummy;
   TrajectoryStateClosestToPoint tt2TSCP = tt2 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt2TSCP.isValid() ) return dummy;
   
   return dca;
   
}

bool KpipiTrueAnalyzer::IsGoodCandidate( const pat::PATObject<reco::VertexCompositeCandidate> & kpipi , const Vertex & primaryVertex, const MagneticField * bField  )
{
   
   if ( kpipi.mass() < KpipiMassWindow_[0] || kpipi.mass() > KpipiMassWindow_[1] ) return false;
   if ( kpipi.pt() < KpipiPtMin_ ) return false;
   if ( kpipi.eta() < KpipiEtaRange_[0] || kpipi.eta() > KpipiEtaRange_[1] ) return false; 
   
   // decay vertex
   Vertex decayVertex(kpipi.vertex(),kpipi.vertexCovariance(),kpipi.vertexChi2(),kpipi.vertexNdof(),kpipi.numberOfDaughters());
   double decayVertexProb = ChiSquaredProbability(decayVertex.chi2(),decayVertex.ndof());
   if ( decayVertexProb < KpipiDecayVertexProbMin_ ) return false;
   GlobalVector kpipiVector(kpipi.px(),kpipi.py(),kpipi.pz());
   SecondaryVertex sv(primaryVertex,decayVertex,kpipiVector,true);
   if ( sv.dist3d().value() < KpipiDecayLengthMin_ ) return false;
   if ( sv.dist3d().value()/sv.dist3d().error() < KpipiDecayLengthSigMin_ ) return false;
   if ( sv.dist2d().value() < KpipiDecayLength2DMin_ ) return false;
   if ( sv.dist2d().value()/sv.dist2d().error() < KpipiDecayLengthSig2DMin_ ) return false;
   
   
   GlobalVector decayVertexVector(kpipi.vx(),kpipi.vy(),kpipi.vz());
   GlobalVector primaryVertexVector(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
   GlobalVector decayLengthVector = decayVertexVector - primaryVertexVector;
   double kpipiPvAngle = acos(fabs(decayLengthVector.dot(kpipiVector)/decayLengthVector.mag()/kpipiVector.mag()));
   if ( kpipiPvAngle > KpipiPVAngleMax_ ) return false;
   
   // Proper time
   double properDecayLength = kpipi.mass()/kpipi.pt()*sv.dist2d().value();
   if ( properDecayLength < KpipiProperDecayLengthMin_ ) return false;
   
   // daughters
   std::vector<RecoChargedCandidate> dplusDaughters;
   std::vector<TrackRef> dplusTracksRef;
   bool isGoodTrack[3];
   
   for ( unsigned int j = 0; j < 3; ++j )
   {
      dplusDaughters.push_back( *(dynamic_cast<const RecoChargedCandidate *> (kpipi.daughter(j))) );
      dplusTracksRef.push_back( dplusDaughters.at(j).track() );
      
      isGoodTrack[j] =
      fabs(dplusTracksRef.at(j) -> dxy(primaryVertex.position())) < trkDxyMax_
      && fabs(dplusTracksRef.at(j) -> dz(primaryVertex.position())) < trkDzMax_
      && dplusTracksRef.at(j) -> normalizedChi2() < trkNormChi2Max_
      && dplusTracksRef.at(j) -> eta() > trkEtaRange_[0] && dplusTracksRef.at(j) -> eta() < trkEtaRange_[1] 
      && dplusTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
      && dplusTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_;
      
      if ( j == 0 ) isGoodTrack[j] += dplusTracksRef.at(j) -> pt() > trkKaonPtMin_ ;
      else          isGoodTrack[j] += dplusTracksRef.at(j) -> pt() > trkPionPtMin_;
   }
   if ( ! isGoodTrack[0] || ! isGoodTrack[1] || ! isGoodTrack[2] ) return false;
   
   TrackRef kaonTrackRef =  dplusDaughters.at(0).track();
   TrackRef pi1TrackRef  =  dplusDaughters.at(1).track();
   TrackRef pi2TrackRef  =  dplusDaughters.at(2).track();
   // Closest approach and crossing points of the tracks
   double dcaKpi1   = DistanceOfClosestApproach(kaonTrackRef,pi1TrackRef,&(*bField)); // Kaon - pion1
   double dcaKpi2   = DistanceOfClosestApproach(kaonTrackRef,pi2TrackRef,&(*bField)); // Kaon - pion2
   double dcapi1pi2 = DistanceOfClosestApproach(pi1TrackRef,pi2TrackRef,&(*bField));  // pion1 - pion2
   
   if ( dcaKpi1 > trkDCAMax_ || dcaKpi2 > trkDCAMax_ || dcapi1pi2 > trkDCAMax_ ) return false;
   
   TLorentzVector kaonVec(dplusDaughters.at(0).px(),dplusDaughters.at(0).py(),dplusDaughters.at(0).pz(),dplusDaughters.at(0).energy());
   TLorentzVector kpipiVec(kpipi.px(),kpipi.py(),kpipi.pz(),kpipi.energy());
   kaonVec.Boost(-kpipiVec.BoostVector());
   
   double cosThStarK = cos(kaonVec.Vect().Angle(kpipiVec.Vect()));
   if ( cosThStarK < KpipiCosThStarKMin_ ) return false;
   
   
   return true;
}

bool KpipiTrueAnalyzer::CheckParentage ( const GenParticle * mcp, int pdg )
{
   bool parentFound = false;
   
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      if ( mom -> pdgId() == pdg ) return true; 
   }
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      parentFound = CheckParentage (mom, pdg);
      if ( parentFound ) return true;
   }
   
   return parentFound;
   
}

bool KpipiTrueAnalyzer::CheckParentage ( const GenParticleRef & mcpRef, int pdg )
{
   bool parentFound = false;
   const GenParticle * mcp = dynamic_cast<const GenParticle *>(mcpRef.get());
   if ( mcp != 0 ) parentFound = CheckParentage (mcp, pdg);
   
   return parentFound; 
}

bool KpipiTrueAnalyzer::CheckParentage ( const Candidate * mcp, int pdg )
{
   bool parentFound = false;
   
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      if ( mom -> pdgId() == pdg ) return true; 
   }
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      parentFound = CheckParentage (mom, pdg);
      if ( parentFound ) return true;
   }
   
   return parentFound; 
}

bool KpipiTrueAnalyzer::FromBdecay ( const Candidate * mcp )
{
   bool fromB = false;
   
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      if ( (abs(mom -> pdgId()) ) == 5 ) return true; 
   }
   
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      fromB = FromBdecay(mom); 
   }
   
   return fromB; 
}

bool KpipiTrueAnalyzer::FromCdecay ( const Candidate * mcp )
{
   bool fromC = false;
   
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      if ( (abs(mom -> pdgId()) ) == 4 && ! FromBdecay(mcp) ) return true; 
   }
   
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      fromC = FromCdecay(mom); 
   }
   
   return fromC; 
}

Vertex KpipiTrueAnalyzer::TruePrimaryVertex ( const Candidate * mcp, bool & foundVertex )
{
   Vertex pv;
   CovarianceMatrix error;
   foundVertex = false;
   
   size_t nMom = mcp -> numberOfMothers();
   for ( size_t i = 0; i < nMom; ++i )
   {
      const Candidate * mom = mcp -> mother(i);
      if ( (abs(mom -> pdgId()) ) <= 6 )
      {
         Vertex vertex ( mom -> vertex(), error , mom -> vertexChi2(), mom -> vertexNdof(), 0);
         pv = vertex;
         foundVertex = true; 
      }
   }
   
   if ( foundVertex ) return pv;
   
   while ( ! foundVertex && mcp -> mother() != 0 )
   {
      for ( size_t i = 0; i < nMom; ++i )
      {
         const Candidate * mom = mcp -> mother(i);
         pv = TruePrimaryVertex(mom, foundVertex); 
         if ( foundVertex ) return pv;
      }
   }
   
   return pv; 
}


bool KpipiTrueAnalyzer::AreSiblings (
                                     const Candidate * mom ,
                                     const GenParticleRef & kaonRef ,
                                     const GenParticleRef & pi1Ref ,
                                     const GenParticleRef & pi2Ref )
{
   
   if ( mom -> numberOfDaughters() != 3 ) return false;
   
   const Candidate * kaon = dynamic_cast<const Candidate *>(kaonRef.get());
   const Candidate * pi1  = dynamic_cast<const Candidate *>(pi1Ref.get());
   const Candidate * pi2  = dynamic_cast<const Candidate *>(pi2Ref.get());
   
   if ( abs(kaon -> pdgId()) != 321 || abs(pi1 -> pdgId()) != 211 || abs(pi2 -> pdgId()) != 211 )
      std::cout << " ***  KpipiTrueAnalyzer::AreSiblings: Not correct PDG Ids ***" << std::endl;
   
   int cntSib = 0;
   
   for ( Candidate::const_iterator dau = mom->begin(); dau != mom->end(); ++dau )
      if ( kaon == &*(dau) || pi1 == &*(dau) || pi2 == &*(dau) ) ++cntSib;
   
   if ( cntSib != 3 ) return false;

   return true;
   
}

VertexCompositeCandidate KpipiTrueAnalyzer::TrueKpipi (
                                                       const GenParticleRef & kaonRef ,
                                                       const GenParticleRef & pi1Ref ,
                                                       const GenParticleRef & pi2Ref ,
                                                       const Candidate * dplus
                                                       )
{
   
   const Candidate * kaon = dynamic_cast<const Candidate *>(kaonRef.get());
   const Candidate * pi1  = dynamic_cast<const Candidate *>(pi1Ref.get());
   const Candidate * pi2  = dynamic_cast<const Candidate *>(pi2Ref.get());

   VertexCompositeCandidate kpipi;

   AddFourMomenta addp4;
   kpipi.addDaughter(*kaon);
   kpipi.addDaughter(*pi1);
   kpipi.addDaughter(*pi2);
   kpipi.setPdgId(dplus -> pdgId());
   kpipi.setCharge(dplus -> charge());
   kpipi.setVertex(kaon -> vertex());
   addp4.set( kpipi );

   return kpipi;
   
}

void KpipiTrueAnalyzer::FillHistograms(
                                       const pat::PATObject<reco::VertexCompositeCandidate> & kpipi, 
                                       const VertexCompositeCandidate & trueKpipi, 
                                       const Vertex & primaryVertex,  
                                       const Vertex & truePrimaryVertex,  
                                       const MagneticField * bField )
{
   Histogram_["Mass"]   ->  Fill(kpipi.mass());
   Histogram_["Pt"]     ->  Fill(kpipi.pt());
   Histogram_["Eta"]    ->  Fill(kpipi.eta());
   

   // decay length 
   Vertex decayVertex(kpipi.vertex(),kpipi.vertexCovariance(),kpipi.vertexChi2(),kpipi.vertexNdof(),kpipi.numberOfDaughters());
   double decayVertexProb = ChiSquaredProbability(decayVertex.chi2(),decayVertex.ndof());
   GlobalVector decayVertexVector(kpipi.vx(),kpipi.vy(),kpipi.vz());
   GlobalVector primaryVertexVector(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
   GlobalVector kpipiVector(kpipi.px(),kpipi.py(),kpipi.pz());
   GlobalVector decayLengthVector = decayVertexVector - primaryVertexVector;
   SecondaryVertex sv(primaryVertex,decayVertex,kpipiVector,true);
   double kpipiPvAngle = acos(fabs(decayLengthVector.dot(kpipiVector)/decayLengthVector.mag()/kpipiVector.mag()));
   double properDecayLength = kpipi.mass()/kpipi.pt()*sv.dist2d().value();
   double decayLength2D = sv.dist2d().value();
   double decayLength = sv.dist3d().value();
   
   Histogram_["DecayLength"] -> Fill(sv.dist3d().value());
   Histogram_["DecayLengthSignificance"] -> Fill(sv.dist3d().value()/sv.dist3d().error());
   Histogram_["DecayLength2D"] -> Fill(sv.dist2d().value());
   Histogram_["DecayLengthSignificance2D"] -> Fill(sv.dist2d().value()/sv.dist2d().error());
   Histogram_["DecayVertexProb"] -> Fill(decayVertexProb);
   Histogram_["PVAngle"] -> Fill(kpipiPvAngle);
   Histogram_["ProperDecayLength"] -> Fill(properDecayLength);
   Histogram_["ProperDecayLengthLargeBins"] -> Fill(properDecayLength);
   
   // daughters
   std::vector<RecoChargedCandidate> dplusDaughters;
   std::vector<TrackRef> dplusTracksRef;
   
   for ( unsigned int j = 0; j < 3; ++j )
   {
      dplusDaughters.push_back( *(dynamic_cast<const RecoChargedCandidate *> (kpipi.daughter(j))) );
      dplusTracksRef.push_back( dplusDaughters.at(j).track() );
      
      Histogram_["TrkDxy"] -> Fill(dplusTracksRef.at(j) -> dxy(primaryVertex.position()) );
      Histogram_["TrkDz"] -> Fill(dplusTracksRef.at(j) -> dz(primaryVertex.position()) );
      Histogram_["TrkNormChi2"] -> Fill(dplusTracksRef.at(j) -> normalizedChi2() );
      Histogram_["TrkPt"] -> Fill(dplusTracksRef.at(j) -> pt() );
      Histogram_["TrkEta"] -> Fill(dplusTracksRef.at(j) -> eta() );
      Histogram_["TrkPixelHits"] -> Fill(dplusTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() );
      Histogram_["TrkSiliconHits"] -> Fill(dplusTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() );
   }
   TrackRef kaonTrackRef =  dplusDaughters.at(0).track();
   TrackRef pi1TrackRef  =  dplusDaughters.at(1).track();
   TrackRef pi2TrackRef  =  dplusDaughters.at(2).track();
   
   // Closest approach and crossing points of the tracks
   double dcaKpi1   = DistanceOfClosestApproach(kaonTrackRef,pi1TrackRef,&(*bField)); // Kaon - pion1
   double dcaKpi2   = DistanceOfClosestApproach(kaonTrackRef,pi2TrackRef,&(*bField)); // Kaon - pion2
   double dcapi1pi2 = DistanceOfClosestApproach(pi1TrackRef,pi2TrackRef,&(*bField));  // pion1 - pion2
   
   Histogram_["TrkDCA"] -> Fill(dcaKpi1);
   Histogram_["TrkDCA"] -> Fill(dcaKpi2);
   Histogram_["TrkDCA"] -> Fill(dcapi1pi2);
   
   TLorentzVector kaonVec(dplusDaughters.at(0).px(),dplusDaughters.at(0).py(),dplusDaughters.at(0).pz(),dplusDaughters.at(0).energy());
   TLorentzVector kpipiVec(kpipi.px(),kpipi.py(),kpipi.pz(),kpipi.energy());
   kaonVec.Boost(-kpipiVec.BoostVector());
   double cosThStarK = cos(kaonVec.Vect().Angle(kpipiVec.Vect()));
   Histogram_["CosThStarK"] -> Fill(cosThStarK);

   // True stuff
   
   // decay length 
   CovarianceMatrix error;
   Vertex trueDecayVertex(trueKpipi.vertex(),error,0,0,3);
   GlobalVector trueDecayVertexVector(trueKpipi.vx(),trueKpipi.vy(),trueKpipi.vz());
   GlobalVector truePrimaryVertexVector(truePrimaryVertex.x(),truePrimaryVertex.y(),truePrimaryVertex.z());
   GlobalVector trueKpipiVector(trueKpipi.px(),trueKpipi.py(),trueKpipi.pz());
   GlobalVector trueDecayLengthVector = trueDecayVertexVector - truePrimaryVertexVector;
   SecondaryVertex trueSV(truePrimaryVertex,trueDecayVertex,trueKpipiVector,false);
   double trueDecayLength2D = trueSV.dist2d().value();
   double trueDecayLength = trueSV.dist3d().value();
   double trueProperDecayLength = trueKpipi.mass()/trueKpipi.pt()*trueSV.dist2d().value();
   Histogram_["ProperDecayLengthTrue"] -> Fill(trueProperDecayLength);
   
   // Tree variables
   pt_true_ = trueKpipi.pt();
   pt_rec_ = kpipi.pt();
   m_true_ = trueKpipi.mass();
   m_rec_ = kpipi.mass();
   lxy_true_ = trueDecayLength2D;
   lxy_rec_ = decayLength2D;
   lxyz_true_ = trueDecayLength;
   lxyz_rec_ = decayLength;
   ct_true_ = trueProperDecayLength;
   ct_rec_ = properDecayLength;
   
   Histogram_["ResolutionPt"] -> Fill(pt_true_ - pt_rec_);
   Histogram_["ResolutionMass"] -> Fill(m_true_ - m_rec_);
   Histogram_["ResolutionDecayLength2D"]  -> Fill(lxy_true_ - lxy_rec_);
   Histogram_["ResolutionDecayLength"]  -> Fill(lxyz_true_ - lxyz_rec_);
   Histogram_["ResolutionProperDecayLength"]  -> Fill(ct_true_ - ct_rec_);
   
}

void KpipiTrueAnalyzer::CreateHistograms()
{
   int nMassBins = int((KpipiMassWindow_[1] - KpipiMassWindow_[0])*1000/2);
   
   // Histograms
   //   // D+
   Histogram_["NDplus"] = fs_ -> make<TH1D>( "NDplus", "", 20, 0., 20. );
   Histogram_["Mass"]   = fs_ -> make<TH1D>( "Mass", "", nMassBins, KpipiMassWindow_[0], KpipiMassWindow_[1] );
   Histogram_["Pt"]     = fs_ -> make<TH1D>( "Pt", "", 200, 0., 20. );
   Histogram_["Eta"]    = fs_ -> make<TH1D>( "Eta", "", 200, -3.5, 3.5 );
   //   
   //   
   //   // D+ decay vertex
   Histogram_["DecayVertexProb"]              = fs_ -> make<TH1D>( "DecayVertexProb", "",          1000,  0.,   1. );
   Histogram_["DecayLength"]                  = fs_ -> make<TH1D>( "DecayLength", "",              1000,  0.,   1. );
   Histogram_["DecayLengthSignificance"]      = fs_ -> make<TH1D>( "DecayLengthSignificance", "",   100,  0.,  10. );
   Histogram_["DecayLength2D"]                = fs_ -> make<TH1D>( "DecayLength2D", "",            1000,  0.,   1. );
   Histogram_["DecayLengthSignificance2D"]    = fs_ -> make<TH1D>( "DecayLengthSignificance2D", "", 100,  0.,  10. );
   Histogram_["ProperDecayLength"]            = fs_ -> make<TH1D>( "ProperDecayLength", "",         480, -0.04, 0.2 );
   
   
   //   // D+ tracks
   Histogram_["TrkDxy"]          = fs_ -> make<TH1D>( "TrkDxy", "",        200, -0.2, 0.2 );
   Histogram_["TrkDz"]           = fs_ -> make<TH1D>( "TrkDz", "",        1000, -1.,  1. );
   Histogram_["TrkNormChi2"]     = fs_ -> make<TH1D>( "TrkNormChi2", "",    60,  0.,  3. );
   Histogram_["TrkPt"]           = fs_ -> make<TH1D>( "TrkPt", "",         500,  0., 10. );
   Histogram_["TrkEta"]          = fs_ -> make<TH1D>( "TrkEta", "",        250, -2.5, 2.5 );
   Histogram_["TrkPixelHits"]    = fs_ -> make<TH1D>( "TrkPixelHits", "",   10,  0., 10. );
   Histogram_["TrkSiliconHits"]  = fs_ -> make<TH1D>( "TrkSiliconHits", "", 40,  0., 40. );
   Histogram_["TrkDCA"]          = fs_ -> make<TH1D>( "TrkDCA", "",        100,  0.,  1. );
   
   
   //   Extras
   Histogram_["PVAngle"] = fs_ -> make<TH1D>( "PVAngle", "", 50, 0., 0.5 );
   Histogram_["CosThStarK"] = fs_ -> make<TH1D>( "CosThStarK", "", 100, -1, 1 );
   
   Histogram_["ProperDecayLengthLargeBins"]   = fs_ -> make<TH1D>( "ProperDecayLengthLargeBins", "", 48, -0.04, 0.2 );
   Histogram_["ProperDecayLengthTrue"]        = fs_ -> make<TH1D>( "ProperDecayLengthTrue", "", 48, -0.04, 0.2 );
   
   // Resolution
   Histogram_["ResolutionDecayLength2D"]  = fs_ -> make<TH1D>( "ResolutionDecayLength2D", "", 200, -0.2, 0.2);
   Histogram_["ResolutionDecayLength"]  = fs_ -> make<TH1D>( "ResolutionDecayLength", "", 200, -0.2, 0.2);
   Histogram_["ResolutionProperDecayLength"]  = fs_ -> make<TH1D>( "ResolutionProperDecayLength", "", 200, -0.2, 0.2);
   Histogram_["ResolutionMass"] = fs_ -> make<TH1D>( "ResolutionMass", "", 1000, -0.5,0.5 );
   Histogram_["ResolutionPt"] = fs_ -> make<TH1D>( "ResolutionPt", "", 1000, -1, 1.);
   
   // Generated
   Histogram_["NGen"]      = fs_ -> make<TH1D>( "NGen", "", 10, 0., 10. );
   Histogram_["MassGen"]   = fs_ -> make<TH1D>( "MassGen", "", 40, 1.7, 2.1 );
   Histogram_["PtGen"]     = fs_ -> make<TH1D>( "PtGen", "", 300, 0., 30. );
   Histogram_["EtaGen"]    = fs_ -> make<TH1D>( "EtaGen", "", 200, -10, 10. );
   //   
   
   //-- Create new Tree to store variables necessary for jet resolution study --//
   Tree_ = new TTree("Tree","Tree");
   Tree_ -> Branch("pt_true",  &pt_true_,   "pt_true/D");
   Tree_ -> Branch("m_true",   &m_true_,    "m_true/D");
   Tree_ -> Branch("lxy_true", &lxy_true_,  "lxy_true/D");
   Tree_ -> Branch("lxyz_true",&lxyz_true_, "lxyz_true/D");
   Tree_ -> Branch("ct_true",  &ct_true_,   "ct_true/D");
   
   Tree_ -> Branch("pt_rec",  &pt_rec_,   "pt_rec/D");
   Tree_ -> Branch("m_rec",   &m_rec_,    "m_rec/D");
   Tree_ -> Branch("lxy_rec", &lxy_rec_,  "lxy_rec/D");
   Tree_ -> Branch("lxyz_rec",&lxyz_rec_, "lxyz_rec/D");
   Tree_ -> Branch("ct_rec",  &ct_rec_,   "ct_rec/D");
   
   
}

void KpipiTrueAnalyzer:: SetCuts(const edm::ParameterSet& iConfig)
{
   if ( useCuts_ )
   {
      // Track selection
      trkDxyMax_         = iConfig.getParameter<double> ("trkDxyMax");
      trkDzMax_          = iConfig.getParameter<double> ("trkDzMax");
      trkPixelHitsMin_   = iConfig.getParameter<double> ("trkPixelHitsMin");
      trkSiliconHitsMin_ = iConfig.getParameter<double> ("trkSiliconHitsMin");
      trkNormChi2Max_    = iConfig.getParameter<double> ("trkNormChi2Max");
      trkEtaRange_       = iConfig.getParameter< std::vector<double> > ("trkEtaRange");
      trkKaonPtMin_      = iConfig.getParameter<double> ("trkKaonPtMin");
      trkPionPtMin_      = iConfig.getParameter<double> ("trkKaonPtMin");
      trkDCAMax_         = iConfig.getParameter<double> ("trkDCAMax");
      
      // D+ selection
      KpipiMassWindow_         = iConfig.getParameter< std::vector<double> > ("MassWindow");
      KpipiPtMin_              = iConfig.getParameter<double> ("PtMin");
      KpipiEtaRange_           = iConfig.getParameter< std::vector<double> > ("EtaRange");
      KpipiDecayVertexProbMin_ = iConfig.getParameter<double> ("DecayVertexProbMin");
      KpipiPVAngleMax_         = iConfig.getParameter<double> ("PVAngleMax");
      
      KpipiDecayLengthSigMin_  = iConfig.getParameter<double> ("DecayLengthSigMin");
      KpipiDecayLengthMin_     = iConfig.getParameter<double> ("DecayLengthMin");
      KpipiDecayLengthSig2DMin_  = iConfig.getParameter<double> ("DecayLengthSig2DMin");
      KpipiDecayLength2DMin_     = iConfig.getParameter<double> ("DecayLength2DMin");
      
      KpipiProperDecayLengthMin_     = iConfig.getParameter<double> ("ProperDecayLengthMin");
      
      KpipiCosThStarKMin_     = iConfig.getParameter<double> ("CosThStarKMin");
      
   }
   else
   {
      SetDummyCuts();
   }
}

void KpipiTrueAnalyzer:: SetDummyCuts()
{
   
   double dummyMax = 99999.99;
   double dummyMin = -99999.99;
   int dummyIntMin = -9999;
   
   trkDxyMax_         = dummyMax;
   trkDzMax_          = dummyMax;
   trkPixelHitsMin_   = dummyIntMin;
   trkSiliconHitsMin_ = dummyIntMin;
   trkNormChi2Max_    = dummyMax;
   trkEtaRange_.push_back(dummyMin);
   trkEtaRange_.push_back(dummyMax);
   trkKaonPtMin_      = dummyMin;
   trkPionPtMin_      = dummyMin;
   trkDCAMax_         = dummyMax;
   
   // D+ selection
   KpipiMassWindow_.push_back(1.7);
   KpipiMassWindow_.push_back(2.1);
   KpipiPtMin_               = dummyMin;
   KpipiEtaRange_.push_back(dummyMin);
   KpipiEtaRange_.push_back(dummyMax);
   KpipiDecayVertexProbMin_  = dummyMin;
   KpipiPVAngleMax_          = dummyMax;
   
   KpipiDecayLengthSigMin_    = dummyMin;
   KpipiDecayLengthMin_       = dummyMin;
   KpipiDecayLengthSig2DMin_  = dummyMin;
   KpipiDecayLength2DMin_     = dummyMin;
   
   KpipiProperDecayLengthMin_ = dummyMin;
   
   KpipiCosThStarKMin_        = dummyMin;
   
}



void KpipiTrueAnalyzer::PrintPtEtaPhi ( const Candidate * part )
{
   std::cout << "        PDG " << part -> pdgId() << ": " << part -> pt() << ", " << part -> eta() << ", " << part -> phi() << std::endl;
}

// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(KpipiTrueAnalyzer);

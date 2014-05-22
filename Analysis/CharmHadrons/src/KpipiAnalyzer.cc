// -*- C++ -*-
//
// Package:    KpipiAnalyzer
// Class:      KpipiAnalyzer
// 
/**\class KpipiAnalyzer KpipiAnalyzer.cc Analysis/CharmHadrons/src/KpipiAnalyzer.cc

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"


#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
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

class KpipiAnalyzer : public edm::EDAnalyzer
{
public:
   explicit KpipiAnalyzer(const edm::ParameterSet&);
   ~KpipiAnalyzer();
   
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   std::vector<double> DistanceOfClosestApproach( const TrackRef & , const TrackRef &, const MagneticField * );
   
   void CreateHistograms();
   void SetDummyCuts();
   void SetCuts(const edm::ParameterSet& );
   
   
   // ----------member data ---------------------------
   
   // Config parameters
   InputTag srcKpipi_;
   InputTag srcPVs_;
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
   double trkXpTMax_;
   double trkXpZMax_;
   
   // D+ candidates selection
   double PtMin_;
   std::vector<double> EtaRange_;
   std::vector<double> MassWindow_;
   double DecayVertexProbMin_;
   double PVAngleMax_;
   
   double DecayLengthSigMin_;
   double DecayLengthMin_;
   double DecayLengthSig2DMin_;
   double DecayLength2DMin_;
   
   double ProperDecayLengthMin_;
   
   double CosThStarKMin_;
   
	double firstCTbin_;
	double lastCTbin_;
	double step_;
	int nCTbins_;
	
	
   bool useCuts_;
   
   // Get the file service
   Service<TFileService> fs_;
   
   // Histogram handler
   std::map<std::string, TH1D *> Histogram_;
   std::map<std::string, TH2D *> Histogram2D_;
   //
   
};

KpipiAnalyzer::KpipiAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcKpipi_            = iConfig.getParameter<InputTag> ("KpipiCandidates");
   srcPVs_              = iConfig.getParameter<InputTag> ("PrimaryVertex");
   srcPocas_            = iConfig.getParameter<InputTag> ("PointOfClosestApproach");
   srcPocaErrors_       = iConfig.getParameter<InputTag> ("PointOfClosestApproachError");
   
   //
   useCuts_ = iConfig.getParameter<bool> ("UseCuts");
   
   SetCuts( iConfig );
   
	firstCTbin_ = iConfig.getParameter<double> ("lowFirstCTbin");
	lastCTbin_ = iConfig.getParameter<double> ("upLastCTbin");
	step_ = iConfig.getParameter<double> ("stepCTbin");
	nCTbins_ = int((lastCTbin_-firstCTbin_)/step_);
	
   CreateHistograms();
   
}

KpipiAnalyzer::~KpipiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to for each event  ------------
void KpipiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   ESHandle<MagneticField> bFieldHandler;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandler);
   
   // original primary vertex
   Handle<VertexCollection> pvHandler;
   iEvent.getByLabel(srcPVs_, pvHandler);
   Vertex primaryVertex = pvHandler -> at(0);
   
   // Kpipi
   Handle<VertexCompositeCandidateCollection> kpipiHandler;
   iEvent.getByLabel (srcKpipi_, kpipiHandler);
   
   // Point of closest approach
   Handle<PointCollection> pocaHandler;
   iEvent.getByLabel (srcPocas_, pocaHandler);
   Handle<CovarianceMatrixCollection> pocaErrHandler;
   iEvent.getByLabel (srcPocaErrors_, pocaErrHandler);
   
   Histogram_["Ncands"] -> Fill( kpipiHandler -> size() );
   
   for ( unsigned int i = 0; i < kpipiHandler -> size() ; ++i )
   {
      VertexCompositeCandidate kpipi = kpipiHandler -> at(i);
      
      if ( kpipi.mass() < MassWindow_[0] || kpipi.mass() > MassWindow_[1] ) continue;
      
      // decay vertex
      Vertex decayVertex(kpipi.vertex(),kpipi.vertexCovariance(),kpipi.vertexChi2(),kpipi.vertexNdof(),kpipi.numberOfDaughters());
      double decayVertexProb = ChiSquaredProbability(decayVertex.chi2(),decayVertex.ndof());
      GlobalVector decayVertexVector(kpipi.vx(),kpipi.vy(),kpipi.vz());
      GlobalVector primaryVertexVector(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
      GlobalVector kpipiVector(kpipi.px(),kpipi.py(),kpipi.pz());
      GlobalVector decayLengthVector = decayVertexVector - primaryVertexVector;
      SecondaryVertex sv(primaryVertex,decayVertex,kpipiVector,true);
      double kpipiPvAngle = acos(fabs(decayLengthVector.dot(kpipiVector)/decayLengthVector.mag()/kpipiVector.mag()));
      // Proper time
      double properDecayLength = kpipi.mass()/kpipi.pt()*sv.dist2d().value();
      
      // point of closest approach
      Point poca = pocaHandler -> at(i);
      CovarianceMatrix pocaError = pocaErrHandler -> at(i);
      
      // daughters
      std::vector<RecoChargedCandidate> dplusDaughters;
      std::vector<TrackRef> dplusTracksRef;
      bool isGoodTrack[3];
      
      for ( unsigned int j = 0; j < 3; ++j )
      {
         dplusDaughters.push_back( *(dynamic_cast<const RecoChargedCandidate *> (kpipi.daughter(j))) );
         dplusTracksRef.push_back( dplusDaughters.at(j).track() );
         
         if ( j == 0 )
         isGoodTrack[j] =
         fabs(dplusTracksRef.at(j) -> dxy(primaryVertex.position())) < trkDxyMax_
         && fabs(dplusTracksRef.at(j) -> dz(primaryVertex.position())) < trkDzMax_
         && dplusTracksRef.at(j) -> normalizedChi2() < trkNormChi2Max_
         && dplusTracksRef.at(j) -> eta() > trkEtaRange_[0] && dplusTracksRef.at(j) -> eta() < trkEtaRange_[1] 
         && dplusTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
         && dplusTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_
         && dplusTracksRef.at(j) -> pt() > trkKaonPtMin_;
         else
         isGoodTrack[j] =
         fabs(dplusTracksRef.at(j) -> dxy(primaryVertex.position())) < trkDxyMax_
         && fabs(dplusTracksRef.at(j) -> dz(primaryVertex.position())) < trkDzMax_
         && dplusTracksRef.at(j) -> normalizedChi2() < trkNormChi2Max_
         && dplusTracksRef.at(j) -> eta() > trkEtaRange_[0] && dplusTracksRef.at(j) -> eta() < trkEtaRange_[1]
         && dplusTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
         && dplusTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_
         && dplusTracksRef.at(j) -> pt() > trkPionPtMin_;

         
//         if ( j == 0 ) isGoodTrack[j] += dplusTracksRef.at(j) -> pt() > trkKaonPtMin_ ;
//         else          isGoodTrack[j] += dplusTracksRef.at(j) -> pt() > trkPionPtMin_;
         
      }
      if ( ! isGoodTrack[0] || ! isGoodTrack[1] || ! isGoodTrack[2] ) continue;

      
      TrackRef kaonTrackRef =  dplusDaughters.at(0).track();
      TrackRef pi1TrackRef  =  dplusDaughters.at(1).track();
      TrackRef pi2TrackRef  =  dplusDaughters.at(2).track();

      // Closest approach and crossing points of the tracks
      std::vector<double> dcaKpi1   = DistanceOfClosestApproach(kaonTrackRef,pi1TrackRef,&(*bFieldHandler)); // Kaon - pion1
      if ( dcaKpi1.size() == 1 ) continue;
      std::vector<double> dcaKpi2   = DistanceOfClosestApproach(kaonTrackRef,pi2TrackRef,&(*bFieldHandler)); // Kaon - pion2
      if ( dcaKpi2.size() == 1 ) continue;
      std::vector<double> dcapi1pi2 = DistanceOfClosestApproach(pi1TrackRef,pi2TrackRef,&(*bFieldHandler));  // pion1 - pion2
      if ( dcapi1pi2.size() == 1 ) continue;
      
      TLorentzVector kaonVec(dplusDaughters.at(0).px(),dplusDaughters.at(0).py(),dplusDaughters.at(0).pz(),dplusDaughters.at(0).energy());
      TLorentzVector kpipiVec(kpipi.px(),kpipi.py(),kpipi.pz(),kpipi.energy());
      kaonVec.Boost(-kpipiVec.BoostVector());
      double cosThStarK = cos(kaonVec.Vect().Angle(kpipiVec.Vect()));
      
      // More cuts!!!
      if ( decayVertexProb < DecayVertexProbMin_ ) continue;
      if ( kpipi.pt() < PtMin_ ) continue;
      if ( kpipi.eta() < EtaRange_[0] || kpipi.eta() > EtaRange_[1] ) continue; 
      if ( cosThStarK < CosThStarKMin_ ) continue;
      if ( properDecayLength < ProperDecayLengthMin_ ) continue;
      if ( kpipiPvAngle > PVAngleMax_ )  continue;
      
      // Distributions after the cuts
      for ( unsigned int j = 0; j < 3; ++j )
      {
         Histogram_["TrkDxy"] -> Fill(dplusTracksRef.at(j) -> dxy(primaryVertex.position()) );
         Histogram_["TrkDz"] -> Fill(dplusTracksRef.at(j) -> dz(primaryVertex.position()) );
         Histogram_["TrkNormChi2"] -> Fill(dplusTracksRef.at(j) -> normalizedChi2() );
         Histogram_["TrkEta"] -> Fill(dplusTracksRef.at(j) -> eta() );
         Histogram_["TrkPixelHits"] -> Fill(dplusTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() );
         Histogram_["TrkSiliconHits"] -> Fill(dplusTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() );
      }
      Histogram_["TrkKPt"] -> Fill(dplusTracksRef.at(0) -> pt() );
      Histogram_["TrkPiPt"] -> Fill(dplusTracksRef.at(1) -> pt() );
      Histogram_["TrkPiPt"] -> Fill(dplusTracksRef.at(2) -> pt() );
      Histogram_["TrkDCA"]  -> Fill(dcaKpi1.at(0));
      Histogram_["TrkDCA"]  -> Fill(dcaKpi2.at(0));
      Histogram_["TrkDCA"]  -> Fill(dcapi1pi2.at(0));
      Histogram_["TrkXpT"]  -> Fill(dcaKpi1.at(1));
      Histogram_["TrkXpT"]  -> Fill(dcaKpi2.at(1));
      Histogram_["TrkXpT"]  -> Fill(dcapi1pi2.at(1));
      Histogram_["TrkXpZ"]  -> Fill(dcaKpi1.at(2));
      Histogram_["TrkXpZ"]  -> Fill(dcaKpi2.at(2));
      Histogram_["TrkXpZ"]  -> Fill(dcapi1pi2.at(2));
      Histogram_["DecayLength"] -> Fill(sv.dist3d().value());
      Histogram_["DecayLengthSignificance"] -> Fill(sv.dist3d().value()/sv.dist3d().error());
      Histogram_["DecayVertexProb"] -> Fill(decayVertexProb);
      Histogram_["PVAngle"] -> Fill(kpipiPvAngle);
      Histogram_["DecayLength2D"] -> Fill(sv.dist2d().value());
      Histogram_["DecayLengthError2D"] -> Fill(sv.dist2d().error());
      Histogram2D_["DecayLength2DValuexError"] -> Fill(sv.dist2d().value(),sv.dist2d().error());
      Histogram_["DecayLengthSignificance2D"] -> Fill(sv.dist2d().value()/sv.dist2d().error());
      Histogram_["ProperDecayLength"] -> Fill(properDecayLength);
      Histogram_["Pt"] -> Fill( kpipi.pt() );
      Histogram_["Eta"] -> Fill( kpipi.eta() );
      Histogram_["Mass"] -> Fill( kpipi.mass() );
      
	   for ( int i = 0 ; i < nCTbins_ ; ++i )
	   {
		   double lower = firstCTbin_ + step_*double(i);
		   double upper = lower + step_;
		   if ( properDecayLength > lower && properDecayLength <= upper )
            Histogram_[Form("Mass_CTbin_%03d",i)] -> Fill( kpipi.mass() );
	   }
      // Decay length cuts (for the mass signal only)
      if ( sv.dist3d().value() < DecayLengthMin_ ) continue;
      if ( sv.dist3d().value()/sv.dist3d().error() < DecayLengthSigMin_ ) continue;
      if ( sv.dist2d().value() < DecayLength2DMin_ ) continue;
      if ( sv.dist2d().value()/sv.dist2d().error() < DecayLengthSig2DMin_ ) continue;
      
      Histogram_["MassSignal"] -> Fill( kpipi.mass() );
      
   }
   

}


// ------------ method called once each job just before starting event loop  ------------
void KpipiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpipiAnalyzer::endJob()
{
}

void KpipiAnalyzer::CreateHistograms()
{
   int nMassBins = int((MassWindow_[1] - MassWindow_[0])*1000/2);
   
   // Histograms
   // D+
   Histogram_["Ncands"] = fs_ -> make<TH1D>( "Ncands", "", 20, 0., 20. );
   Histogram_["Mass"]   = fs_ -> make<TH1D>( "Mass", "", nMassBins, MassWindow_[0], MassWindow_[1] );
   Histogram_["MassSignal"]   = fs_ -> make<TH1D>( "MassSignal", "", nMassBins, MassWindow_[0], MassWindow_[1] );
   Histogram_["Pt"]     = fs_ -> make<TH1D>( "Pt", "", 200, 0., 20. );
   Histogram_["Eta"]    = fs_ -> make<TH1D>( "Eta", "", 200, -2.5, 2.5 );
   
   // D+ decay vertex
   Histogram_["DecayVertexProb"]              = fs_ -> make<TH1D>( "DecayVertexProb", "",          1000,  0.,   1. );
   Histogram_["DecayLength"]                  = fs_ -> make<TH1D>( "DecayLength", "",              1000,  0.,   1. );
   Histogram_["DecayLengthSignificance"]      = fs_ -> make<TH1D>( "DecayLengthSignificance", "",   100,  0.,  10. );
   Histogram_["DecayLength2D"]                = fs_ -> make<TH1D>( "DecayLength2D", "",            1000,  0.,   1. );
   Histogram_["DecayLengthError2D"]           = fs_ -> make<TH1D>( "DecayLengthError2D", "",       1000,  0.,   0.15 );
   Histogram_["DecayLengthSignificance2D"]    = fs_ -> make<TH1D>( "DecayLengthSignificance2D", "", 100,  0.,  10. );
   Histogram_["ProperDecayLength"]            = fs_ -> make<TH1D>( "ProperDecayLength", "",         480, -0.04, 0.2 );
   
   //   // D+ tracks
   Histogram_["TrkDxy"]          = fs_ -> make<TH1D>( "TrkDxy", "",        200, -0.2, 0.2 );
   Histogram_["TrkDz"]           = fs_ -> make<TH1D>( "TrkDz", "",        1000, -1.,  1. );
   Histogram_["TrkNormChi2"]     = fs_ -> make<TH1D>( "TrkNormChi2", "",    60,  0.,  3. );
   Histogram_["TrkKPt"]          = fs_ -> make<TH1D>( "TrkKPt", "",        500,  0., 10. );
   Histogram_["TrkPiPt"]         = fs_ -> make<TH1D>( "TrkPiPt", "",       500,  0., 10. );
   Histogram_["TrkEta"]          = fs_ -> make<TH1D>( "TrkEta", "",        250, -2.5, 2.5 );
   Histogram_["TrkPixelHits"]    = fs_ -> make<TH1D>( "TrkPixelHits", "",   10,  0., 10. );
   Histogram_["TrkSiliconHits"]  = fs_ -> make<TH1D>( "TrkSiliconHits", "", 40,  0., 40. );
   Histogram_["TrkDCA"]          = fs_ -> make<TH1D>( "TrkDCA", "",        100,  0.,  1. );
   Histogram_["TrkXpT"]          = fs_ -> make<TH1D>( "TrkXpT", "",        500,  0.,  500. );
   Histogram_["TrkXpZ"]          = fs_ -> make<TH1D>( "TrkXpZ", "",        500,  0.,  500. );

   
   // Bins of proper time
	for ( int i = 0 ; i < nCTbins_; ++i )
	{
		double lower = firstCTbin_ + step_*double(i);
		double upper = lower + step_;
		Histogram_[Form("Mass_CTbin_%03d",i)] = fs_ -> make<TH1D>( Form("Mass_CTbin_%03d",i),
																Form("Kpipi mass distribution in ct bin (%5.0f, %5.0f] um",lower*10000,upper*10000),
																nMassBins, MassWindow_[0], MassWindow_[1] );
	}
	
   // Extras
   Histogram_["PVAngle"] = fs_ -> make<TH1D>( "PVAngle", "", 50, 0., 0.5 );
   Histogram_["CosThStarK"] = fs_ -> make<TH1D>( "CosThStarK", "", 100, -1, 1 );
   
   Histogram2D_["DecayLength2DValuexError"] =  fs_ -> make<TH2D>( "DecayLength2DValuexError", "", 1000,  0.,   1., 1000,  0.,   1.5 );
   
}

void KpipiAnalyzer:: SetCuts(const edm::ParameterSet& iConfig)
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
      trkPionPtMin_      = iConfig.getParameter<double> ("trkPionPtMin");
      trkDCAMax_         = iConfig.getParameter<double> ("trkDCAMax");
      
      trkXpTMax_         = iConfig.getParameter<double> ("trkXpTMax");
      trkXpZMax_         = iConfig.getParameter<double> ("trkXpZMax");
      
      // D+ selection
      MassWindow_         = iConfig.getParameter< std::vector<double> > ("MassWindow");
      PtMin_              = iConfig.getParameter<double> ("PtMin");
      EtaRange_           = iConfig.getParameter< std::vector<double> > ("EtaRange");
      DecayVertexProbMin_ = iConfig.getParameter<double> ("DecayVertexProbMin");
      PVAngleMax_         = iConfig.getParameter<double> ("PVAngleMax");
      
      DecayLengthSigMin_    = iConfig.getParameter<double> ("DecayLengthSigMin");
      DecayLengthMin_       = iConfig.getParameter<double> ("DecayLengthMin");
      DecayLengthSig2DMin_  = iConfig.getParameter<double> ("DecayLengthSig2DMin");
      DecayLength2DMin_     = iConfig.getParameter<double> ("DecayLength2DMin");
      
      ProperDecayLengthMin_ = iConfig.getParameter<double> ("ProperDecayLengthMin");
      
      CosThStarKMin_        = iConfig.getParameter<double> ("CosThStarKMin");
   }
   else
   {
      SetDummyCuts();
   }
}

void KpipiAnalyzer:: SetDummyCuts()
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
   
   trkXpTMax_         = dummyMax;
   trkXpZMax_         = dummyMax;
   
   // D+ selection
   MassWindow_.push_back(1.7);
   MassWindow_.push_back(2.1);
   PtMin_               = dummyMin;
   EtaRange_.push_back(dummyMin);
   EtaRange_.push_back(dummyMax);
   DecayVertexProbMin_  = dummyMin;
   PVAngleMax_          = dummyMax;
   
   DecayLengthSigMin_    = dummyMin;
   DecayLengthMin_       = dummyMin;
   DecayLengthSig2DMin_  = dummyMin;
   DecayLength2DMin_     = dummyMin;
   
   ProperDecayLengthMin_ = dummyMin;
   
   CosThStarKMin_        = dummyMin;
   
}

std::vector<double> KpipiAnalyzer::DistanceOfClosestApproach( const TrackRef & tRef1, const TrackRef & tRef2, const MagneticField * bField )
{
   double dca;
   std::vector<double> dummy;
   dummy.push_back(-99999.99);
   
   TransientTrack t1 (tRef1, &(*bField) );
   TransientTrack* tt1  = &t1;
   if ( ! tt1 -> impactPointTSCP().isValid() ) return dummy;
   TransientTrack t2 (tRef2, &(*bField) );
   TransientTrack* tt2  = &t2;
   if ( ! tt2 -> impactPointTSCP().isValid() ) return dummy;

   FreeTrajectoryState tt1State  = tt1 -> impactPointTSCP().theState();
   FreeTrajectoryState tt2State  = tt2 -> impactPointTSCP().theState();
   // Closest approach and crossing points of the tracks
   ClosestApproachInRPhi cApp;
   GlobalPoint cxPt;
   // track1- track2
   cApp.calculate(tt1State, tt2State);
   if ( ! cApp.status() ) return dummy;
   dca = fabs( cApp.distance() );
   if ( dca < 0. || dca > trkDCAMax_ ) return dummy;
   cxPt = cApp.crossingPoint();
   if ( cxPt.transverse() > trkXpTMax_ || fabs(cxPt.z()) > trkXpZMax_ ) return dummy;
   TrajectoryStateClosestToPoint tt1TSCP = tt1 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt1TSCP.isValid() ) return dummy;
   TrajectoryStateClosestToPoint tt2TSCP = tt2 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt2TSCP.isValid() ) return dummy;
   
   std::vector<double> results;
   results.push_back(dca);
   results.push_back(cxPt.transverse());
   results.push_back(cxPt.z());
   return results;
   
}


// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(KpipiAnalyzer);

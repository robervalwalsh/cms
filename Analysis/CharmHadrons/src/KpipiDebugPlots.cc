// -*- C++ -*-
//
// Package:    KpipiDebugPlots
// Class:      KpipiDebugPlots
// 
/**\class KpipiDebugPlots KpipiDebugPlots.cc Analysis/CharmHadrons/src/KpipiDebugPlots.cc

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

class KpipiDebugPlots : public edm::EDAnalyzer
{
public:
   explicit KpipiDebugPlots(const edm::ParameterSet&);
   ~KpipiDebugPlots();
   
private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   double DistanceOfClosestApproach( const TrackRef & , const TrackRef &, const MagneticField * );
   
   void CreateHistograms();
   
   
   // ----------member data ---------------------------
   
   // Config parameters
   InputTag srcKpipi_;
   InputTag srcPVs_;
   InputTag srcPocas_;
   InputTag srcPocaErrors_;
      
   // Get the file service
   Service<TFileService> fs_;
   
   // Histogram handler
   std::map<std::string, TH1D *> Histogram_;
   std::map<std::string, TH2D *> Histogram2D_;
   //
   
};

KpipiDebugPlots::KpipiDebugPlots(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcKpipi_            = iConfig.getParameter<InputTag> ("KpipiCandidates");
   srcPVs_              = iConfig.getParameter<InputTag> ("PrimaryVertex");
   srcPocas_            = iConfig.getParameter<InputTag> ("PointOfClosestApproach");
   srcPocaErrors_       = iConfig.getParameter<InputTag> ("PointOfClosestApproachError");
   //
   CreateHistograms();
   
}

KpipiDebugPlots::~KpipiDebugPlots()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to for each event  ------------
void KpipiDebugPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
      
      TLorentzVector kaonVec(dplusDaughters.at(0).px(),dplusDaughters.at(0).py(),dplusDaughters.at(0).pz(),dplusDaughters.at(0).energy());
      TLorentzVector kpipiVec(kpipi.px(),kpipi.py(),kpipi.pz(),kpipi.energy());
      kaonVec.Boost(-kpipiVec.BoostVector());
      double cosThStarK = cos(kaonVec.Vect().Angle(kpipiVec.Vect()));
      
      // Closest approach and crossing points of the tracks
      DistanceOfClosestApproach(kaonTrackRef,pi1TrackRef,&(*bFieldHandler)); // Kaon - pion1
      DistanceOfClosestApproach(kaonTrackRef,pi2TrackRef,&(*bFieldHandler)); // Kaon - pion2
      DistanceOfClosestApproach(pi1TrackRef,pi2TrackRef,&(*bFieldHandler));  // pion1 - pion2
      
      // Distributions after the cuts
      Histogram_["CosThStarK"]  -> Fill(cosThStarK);
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
   }
   

}


// ------------ method called once each job just before starting event loop  ------------
void KpipiDebugPlots::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpipiDebugPlots::endJob()
{
}

void KpipiDebugPlots::CreateHistograms()
{
   // Histograms
   // D+
   Histogram_["Ncands"] = fs_ -> make<TH1D>( "Ncands", "", 20, 0., 20. );
   Histogram_["Mass"]   = fs_ -> make<TH1D>( "Mass", "", 400, 1.5, 2.3 );
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
   Histogram_["TrkPt"]           = fs_ -> make<TH1D>( "TrkPt", "",         500,  0., 10. );
   Histogram_["TrkEta"]          = fs_ -> make<TH1D>( "TrkEta", "",        250, -2.5, 2.5 );
   Histogram_["TrkPixelHits"]    = fs_ -> make<TH1D>( "TrkPixelHits", "",   10,  0., 10. );
   Histogram_["TrkSiliconHits"]  = fs_ -> make<TH1D>( "TrkSiliconHits", "", 40,  0., 40. );
   Histogram_["TrkDCA"]          = fs_ -> make<TH1D>( "TrkDCA", "",        100,  0.,  1. );
   Histogram_["TrkXpT"]          = fs_ -> make<TH1D>( "TrkXpT", "",        500,  0.,  500. );
   Histogram_["TrkXpZ"]          = fs_ -> make<TH1D>( "TrkXpZ", "",        500,  0.,  500. );
   
   // Extras
   Histogram_["PVAngle"] = fs_ -> make<TH1D>( "PVAngle", "", 50, 0., 0.5 );
   Histogram_["CosThStarK"] = fs_ -> make<TH1D>( "CosThStarK", "", 100, -1, 1 );
   
   Histogram2D_["DecayLength2DValuexError"] =  fs_ -> make<TH2D>( "DecayLength2DValuexError", "", 1000,  0.,   1., 1000,  0.,   1.5 );
   
}


double KpipiDebugPlots::DistanceOfClosestApproach( const TrackRef & tRef1, const TrackRef & tRef2, const MagneticField * bField )
{
   double dca;
   double dummy = -99999.;
   
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
   cxPt = cApp.crossingPoint();
   TrajectoryStateClosestToPoint tt1TSCP = tt1 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt1TSCP.isValid() ) return dummy;
   TrajectoryStateClosestToPoint tt2TSCP = tt2 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt2TSCP.isValid() ) return dummy;
   
   Histogram_["TrkDCA"] -> Fill(dca);
   Histogram_["TrkXpT"] -> Fill(cxPt.transverse());
   Histogram_["TrkXpZ"] -> Fill(cxPt.z());
   
   return dca;
   
}


// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(KpipiDebugPlots);

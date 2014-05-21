// -*- C++ -*-
//
// Package:    KpiAnalyzer
// Class:      KpiAnalyzer
// 
/**\class KpiAnalyzer KpiAnalyzer.cc Analysis/CharmHadrons/src/KpiAnalyzer.cc

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

using namespace edm;
using namespace reco;

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > CovarianceMatrix;
typedef std::vector<CovarianceMatrix> CovarianceMatrixCollection;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;
typedef std::vector<Point> PointCollection;
typedef ROOT::Math::SVector<double, 3> Vector3D;


//
// class declaration

class KpiAnalyzer : public edm::EDAnalyzer
{
   public:
      explicit KpiAnalyzer(const edm::ParameterSet&);
      ~KpiAnalyzer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
   
      double DistanceOfClosestApproach( const TrackRef & , const TrackRef &, const MagneticField * );


      // ----------member data ---------------------------
      
      // Config parameters
      InputTag srcKpi_;
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
      
      // D+ candidates selection
      double KpiPtMin_;
      std::vector<double> KpiEtaRange_;
      std::vector<double> KpiMassWindow_;
      double KpiDecayVertexProbMin_;
      double KpiPVAngleMax_;
      
      double KpiDecayLengthSigMin_;
      double KpiDecayLengthMin_;
   
 
      // Get the file service
      Service<TFileService> fs_;
      
      // Histogram handler
      std::map<std::string, TH1D *> Histogram_;
      //
      
};

KpiAnalyzer::KpiAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   srcKpi_            = iConfig.getParameter<InputTag> ("KpiCandidates");
   srcPVs_              = iConfig.getParameter<InputTag> ("PrimaryVertex");
   srcPocas_            = iConfig.getParameter<InputTag> ("PointOfClosestApproach");
   srcPocaErrors_       = iConfig.getParameter<InputTag> ("PointOfClosestApproachError");
   
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
   KpiMassWindow_         = iConfig.getParameter< std::vector<double> > ("MassWindow");
   KpiPtMin_              = iConfig.getParameter<double> ("PtMin");
   KpiEtaRange_           = iConfig.getParameter< std::vector<double> > ("EtaRange");
   KpiDecayVertexProbMin_ = iConfig.getParameter<double> ("DecayVertexProbMin");
   KpiPVAngleMax_         = iConfig.getParameter<double> ("PVAngleMax");
   
   KpiDecayLengthSigMin_  = iConfig.getParameter<double> ("DecayLengthSigMin");
   KpiDecayLengthMin_     = iConfig.getParameter<double> ("DecayLengthMin");
   
// Histograms
   // D+
   Histogram_["Ncands"] = fs_ -> make<TH1D>( "Ncands", "", 20, 0., 20. );
   Histogram_["Mass"]   = fs_ -> make<TH1D>( "Mass", "", 210, 1.5, 2.2 );
   Histogram_["Pt"]     = fs_ -> make<TH1D>( "Pt", "", 200, 0., 20. );
   Histogram_["Eta"]    = fs_ -> make<TH1D>( "Eta", "", 200, -2.5, 2.5 );
   
   Histogram_["PVAngle"] = fs_ -> make<TH1D>( "PVAngle", "", 50, 0., 0.5 );
   
   // D+ decay vertex
   Histogram_["DecayVertexProb"]    = fs_ -> make<TH1D>( "DecayVertexProb", "", 1000, 0., 1. );
   Histogram_["DecayLength"]        = fs_ -> make<TH1D>( "DecayLength", "", 1000, 0., 1. );
   Histogram_["DecayLengthSignificance"]        = fs_ -> make<TH1D>( "DecayLengthSignificance", "", 100, 0., 10. );

   // D+ tracks
   Histogram_["TrkDxy"]      = fs_ -> make<TH1D>( "TrkDxy", "", 200, -0.2, 0.2 );
   Histogram_["TrkDz"]       = fs_ -> make<TH1D>( "TrkDz", "", 1000, -1., 1. );
   Histogram_["TrkNormChi2"] = fs_ -> make<TH1D>( "TrkNormChi2", "", 60, 0., 3. );
   Histogram_["TrkPt"]       = fs_ -> make<TH1D>( "TrkPt", "", 500, 0., 10. );
   Histogram_["TrkEta"]      = fs_ -> make<TH1D>( "TrkEta", "", 250, -2.5, 2.5 );
   Histogram_["TrkPixelHits"]    = fs_ -> make<TH1D>( "TrkPixelHits", "", 10, 0., 10. );
   Histogram_["TrkSiliconHits"]  = fs_ -> make<TH1D>( "TrkSiliconHits", "", 40, 0., 40. );
   
   Histogram_["TrkDCAkpi"]   = fs_ -> make<TH1D>( "TrkDCAkpi", "", 100, 0., 1. );
   
   
}

KpiAnalyzer::~KpiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to for each event  ------------
void KpiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   GlobalVector nullVector(0,0,0);
   
   ESHandle<MagneticField> bFieldHandler;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandler);
   
   try
   {
      
      // original primary vertex
      Handle<VertexCollection> pvHandler;
      iEvent.getByLabel(srcPVs_, pvHandler);
      Vertex primaryVertex = pvHandler -> at(0);
      
      
      // Kpi
      Handle<VertexCompositeCandidateCollection> kpiHandler;
      iEvent.getByLabel (srcKpi_, kpiHandler);
      
      // Point of closest approach
      Handle<PointCollection> pocaHandler;
      iEvent.getByLabel (srcPocas_, pocaHandler);
      Handle<CovarianceMatrixCollection> pocaErrHandler;
      iEvent.getByLabel (srcPocaErrors_, pocaErrHandler);
      
      Histogram_["Ncands"] -> Fill( kpiHandler -> size() );
      
      for ( unsigned int i = 0; i < kpiHandler -> size() ; ++i )
      {
         VertexCompositeCandidate kpi = kpiHandler -> at(i);
         
         // decay vertex
         Vertex decayVertex(kpi.vertex(),kpi.vertexCovariance(),kpi.vertexChi2(),kpi.vertexNdof(),kpi.numberOfDaughters());
         SecondaryVertex sv(primaryVertex,decayVertex,nullVector,true);
         double decayVertexProb = ChiSquaredProbability(decayVertex.chi2(),decayVertex.ndof());
         GlobalVector decayVertexVector(kpi.vx(),kpi.vy(),kpi.vz());
         GlobalVector primaryVertexVector(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
         GlobalVector kpiVector(kpi.px(),kpi.py(),kpi.pz());
         GlobalVector decayLengthVector = decayVertexVector - primaryVertexVector;
         double kpiPvAngle = acos(decayLengthVector.dot(kpiVector)/decayLengthVector.mag()/kpiVector.mag());
         Histogram_["DecayLength"] -> Fill(sv.dist3d().value());
         Histogram_["DecayLengthSignificance"] -> Fill(sv.dist3d().value()/sv.dist3d().error());
         Histogram_["DecayVertexProb"] -> Fill(decayVertexProb);
         Histogram_["PVAngle"] -> Fill(kpiPvAngle);
         
         // point of closest approach
         Point poca = pocaHandler -> at(i);
         CovarianceMatrix pocaError = pocaErrHandler -> at(i);
         
         // daughters
         std::vector<RecoChargedCandidate> kpiDaughters;
         std::vector<TrackRef> kpiTracksRef;
         bool isGoodTrack[2];
         for ( unsigned int j = 0; j < 2; ++j )
         {
            kpiDaughters.push_back( *(dynamic_cast<const RecoChargedCandidate *> (kpi.daughter(j))) );
            kpiTracksRef.push_back( kpiDaughters.at(j).track() );
            
            Histogram_["TrkDxy"] -> Fill(kpiTracksRef.at(j) -> dxy(primaryVertex.position()) );
            Histogram_["TrkDz"] -> Fill(kpiTracksRef.at(j) -> dz(primaryVertex.position()) );
            Histogram_["TrkNormChi2"] -> Fill(kpiTracksRef.at(j) -> normalizedChi2() );
            Histogram_["TrkPt"] -> Fill(kpiTracksRef.at(j) -> pt() );
            Histogram_["TrkEta"] -> Fill(kpiTracksRef.at(j) -> eta() );
            Histogram_["TrkPixelHits"] -> Fill(kpiTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() );
            Histogram_["TrkSiliconHits"] -> Fill(kpiTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() );
            
            isGoodTrack[j] =
            fabs(kpiTracksRef.at(j) -> dxy(primaryVertex.position())) < trkDxyMax_
            && fabs(kpiTracksRef.at(j) -> dz(primaryVertex.position())) < trkDzMax_
            && kpiTracksRef.at(j) -> normalizedChi2() < trkNormChi2Max_
            && kpiTracksRef.at(j) -> eta() > trkEtaRange_[0] && kpiTracksRef.at(j) -> eta() < trkEtaRange_[1] 
            && kpiTracksRef.at(j) -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
            && kpiTracksRef.at(j) -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_;
            
            
         }
         
         isGoodTrack[0] = isGoodTrack[0] && kpiTracksRef.at(0) -> pt() > trkKaonPtMin_ ;
         isGoodTrack[1] = isGoodTrack[1] && kpiTracksRef.at(1) -> pt() > trkPionPtMin_;
         
         TrackRef kaonTrackRef =  kpiDaughters.at(0).track();
         TrackRef pionTrackRef =  kpiDaughters.at(1).track();
         
         // Closest approach and crossing points of the tracks
         double dcaKpi   = DistanceOfClosestApproach(kaonTrackRef,pionTrackRef,&(*bFieldHandler)); // Kaon - pion1
         
         Histogram_["TrkDCAkpi"] -> Fill(dcaKpi);
         
         
         // CUTS!!!
         if ( ! isGoodTrack[0] || ! isGoodTrack[1] ) continue;
         if ( decayVertexProb < KpiDecayVertexProbMin_ ) continue;
         if ( kpi.pt() < KpiPtMin_ ) continue;
         if ( kpi.mass() < KpiMassWindow_[0] || kpi.mass() > KpiMassWindow_[1] ) continue;
         if ( kpi.eta() < KpiEtaRange_[0] || kpi.eta() > KpiEtaRange_[1] ) continue; 
         if ( sv.dist3d().value() < KpiDecayLengthMin_ ) continue;
         if ( sv.dist3d().value()/sv.dist3d().error() < KpiDecayLengthSigMin_ ) continue;
         if ( kpiPvAngle > KpiPVAngleMax_ ) continue;
         if ( dcaKpi > trkDCAMax_ ) continue;
         
         
         // Distributions after the cuts
         
         Histogram_["Pt"] -> Fill( kpi.pt() );
         Histogram_["Eta"] -> Fill( kpi.eta() );
         Histogram_["Mass"] -> Fill( kpi.mass() );
         
      }
   } catch (...) {
      
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void KpiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpiAnalyzer::endJob()
{
}


double KpiAnalyzer::DistanceOfClosestApproach( const TrackRef & tRef1, const TrackRef & tRef2, const MagneticField * bField )
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


// _____________________________________________________________________________________________________________


//define this as a plug-in
DEFINE_FWK_MODULE(KpiAnalyzer);

// -*- C++ -*-
//
// Package:    KpiProducer
// Class:      KpiProducer
// 
/**\class KpiProducer KpiProducer.cc CharmHadrons/D0/src/KpiProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Roberval Walsh,01b/352,4849,
// Modifications  :  Alexei Raspereza,01b/352
//         Created:  Wed Jan 27 10:40:28 CET 2010
//        Modified:  Wed Apr 28 11:30:16 CET 2010 
// $Id$
//
//


// system include files
#include <memory>

#include <cmath>

#include <Math/GenVector/PxPyPzM4D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h" 
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"


//#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "DataFormats/GeometrySurface/interface/Line.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "TMath.h"

using namespace edm;
using namespace reco;

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > CovarianceMatrix;
typedef std::vector<CovarianceMatrix> CovarianceMatrixCollection;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;
typedef std::vector<Point> PointCollection;
typedef ROOT::Math::SVector<double, 3> Vector3D;

//
// class decleration
//

class KpiProducer : public edm::EDProducer {
   public:
      explicit KpiProducer(const edm::ParameterSet&);
      ~KpiProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
   
   bool IsGoodTrack ( TrackRef, Vertex );
   bool IsGoodBtwTwoTracks( const TransientTrack*, const TransientTrack*);
   VertexCompositeCandidate RefittedKpi( const TrackRef&, const TrackRef&, const MagneticField *, RefCountedKinematicParticle &, bool & );
   Vertex RefittedPrimaryVertex ( const Vertex&, const TrackRef&, const TrackRef&, const MagneticField *, const BeamSpot&, bool & );
   
      
  // ----------member data ---------------------------
      edm::InputTag srcTracks_;
      edm::InputTag srcPVs_;
      edm::InputTag srcBS_;
   
  // Track selection
      double trkDxyMax_;
      double trkDzMax_;
      int    trkPixelHitsMin_;
      int    trkSiliconHitsMin_;
      double trkNormChi2Max_;
      std::vector<double> trkEtaRange_;
      double trkPtMin_;
      double trkDCAMax_;
   
  // D0 candidates selection
      double KpiPtMin_;
      std::vector<double> KpiEtaRange_;
      std::vector<double> KpiMassWindow_;
      double KpiDecayVertexProbMin_;
   
      bool refitPrimaryVertex_;
      double KpiDecayLengthSigMin_;
      double KpiDecayLengthMin_;
      int KpiPDG_;

   
      ParticleMass pion_mass_;
      ParticleMass proton_mass_;
      ParticleMass kaon_mass_;
      ParticleMass lambda_mass_;
      float pion_sigma_ ;
      float proton_sigma_;
      float kaon_sigma_;
      float lambda_sigma_;
   
   
};

KpiProducer::KpiProducer(const edm::ParameterSet& iConfig)
{
   srcTracks_  = iConfig.getParameter<InputTag> ("Tracks");
   srcPVs_     = iConfig.getParameter<InputTag> ("PrimaryVertices");
   srcBS_      = iConfig.getParameter<InputTag> ("BeamSpot");
   
   // Track selection
   trkDxyMax_         = iConfig.getParameter<double> ("trkDxyMax");
   trkDzMax_          = iConfig.getParameter<double> ("trkDzMax");
   trkPixelHitsMin_   = iConfig.getParameter<int> ("trkPixelHitsMin");
   trkSiliconHitsMin_ = iConfig.getParameter<int> ("trkSiliconHitsMin");
   trkNormChi2Max_    = iConfig.getParameter<double> ("trkNormChi2Max");
   trkEtaRange_       = iConfig.getParameter< std::vector<double> > ("trkEtaRange");
   trkPtMin_          = iConfig.getParameter<double> ("trkPtMin");
   trkDCAMax_         = iConfig.getParameter<double> ("trkDCAMax");
   
   // Kpi selection
   KpiMassWindow_         = iConfig.getParameter< std::vector<double> > ("MassWindow");
   KpiPtMin_              = iConfig.getParameter<double> ("PtMin");
   KpiEtaRange_           = iConfig.getParameter< std::vector<double> > ("EtaRange");
   KpiDecayVertexProbMin_ = iConfig.getParameter<double> ("DecayVertexProbMin");
   
   KpiDecayLengthSigMin_  = iConfig.getParameter<double> ("DecayLengthSigMin");
   KpiDecayLengthMin_     = iConfig.getParameter<double> ("DecayLengthMin");
   
   KpiPDG_                = iConfig.getParameter<int> ("PDG");
   
   refitPrimaryVertex_   = iConfig.getParameter<bool> ("RefitPrimaryVertex");
   
   pion_mass_    = 0.13957018;
   proton_mass_  = 0.938272013;
   kaon_mass_    = 0.493677;
   lambda_mass_  = 1.115683;
   pion_sigma_   = pion_mass_*1.e-6;
   proton_sigma_ = proton_mass_*1.e-6;
   kaon_sigma_   = kaon_mass_*1.e-6;
   lambda_sigma_ = 0.000006;
   
   //register your products
   produces<VertexCompositeCandidateCollection>();
   produces<VertexCollection>("RefittedPrimaryVertex");
   produces<VertexCollection>("OriginalPrimaryVertex");
   produces<CovarianceMatrixCollection>("PointOfClosestApproachError");
   produces<PointCollection>("PointOfClosestApproach");
   
}


KpiProducer::~KpiProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
KpiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   bool valid;
   
   ESHandle<MagneticField> bFieldHandler;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandler);
   
   try {

      // Primary vertices
      Handle<VertexCollection> pvHandler;
      iEvent.getByLabel(srcPVs_, pvHandler);
      unsigned int nPVs = pvHandler -> size();
      Vertex bestPV;
      Vertex refittedPrimaryVertex;
      unsigned int nPVTrks = 0;
      for (unsigned int i = 0; i < nPVs; ++i )
      {
         Vertex pv = pvHandler -> at(i);
         if ( nPVTrks < pv.tracksSize() )
         {
            bestPV = pv;
            nPVTrks = pv.tracksSize();
         }
      }
      
      // Beam Spot
      BeamSpot beamSpot;
      Handle<reco::BeamSpot> beamSpotHandler;
      iEvent.getByLabel(srcBS_, beamSpotHandler);
      if ( beamSpotHandler.isValid() ) beamSpot = *beamSpotHandler; 
      
      // Tracks
      Handle<TrackCollection> trackHandler;
      iEvent.getByLabel(srcTracks_,trackHandler);
      unsigned int nTracks = trackHandler -> size();
      TrackCollection tracks;
      tracks.clear();
      tracks.insert( tracks.end(), trackHandler->begin(), trackHandler->end() );
      
      std::vector<VertexCompositeCandidate> Kpis;
      std::vector<Vertex> OriginalPrimaryVertices;
      std::vector<Vertex> RefittedPrimaryVertices;
      std::vector<Point> PoCAs;
      std::vector<CovarianceMatrix> PoCAErrors;
      
      
      for ( size_t i = 0; i < nTracks; ++i ) // kaon track
      {
         // Kaon candidate
         TrackRef kaonTrackRef( trackHandler, i );
         if ( ! IsGoodTrack(kaonTrackRef,bestPV) ) continue;
         // Transient tracks
         TransientTrack kaonTT(kaonTrackRef, &(*bFieldHandler) );
         TransientTrack* kaonTTPtr = &kaonTT;
         if ( ! kaonTTPtr -> impactPointTSCP().isValid() ) continue;
         // pion candidate
         for ( size_t j = 0; j < nTracks; ++j ) // pion track
         {
            if ( j == i ) continue;  // same track
            TrackRef pionTrackRef( trackHandler, j );
            if ( ! IsGoodTrack(pionTrackRef,bestPV) ) continue;
            if ( kaonTrackRef->charge()*pionTrackRef->charge() > 0 ) continue;  // remove wrong charge combintations
            // Transient tracks
            TransientTrack pionTT (pionTrackRef, &(*bFieldHandler) );
            TransientTrack* pionTTPtr  = &pionTT;
            if( ! pionTTPtr  -> impactPointTSCP().isValid() ) continue;
            // check status between 2 tracks, dca etc.
            if ( ! IsGoodBtwTwoTracks(kaonTTPtr,pionTTPtr) ) continue;
            
            // Reconstruct K-pi candidate
            valid = true;
            RefCountedKinematicParticle KpiCand;
            VertexCompositeCandidate Kpi = RefittedKpi(kaonTrackRef, pionTrackRef, &(*bFieldHandler), KpiCand, valid);
            if ( ! valid ) continue;
            CovarianceMatrix KpiDecayVertexError = Kpi.vertexCovariance();
            
            // Refit primary vertex
            valid = true;
            refittedPrimaryVertex = bestPV;
            if ( refitPrimaryVertex_ )
               refittedPrimaryVertex = RefittedPrimaryVertex( bestPV, kaonTrackRef, pionTrackRef, &(*bFieldHandler), beamSpot,  valid);
            if ( ! valid ) continue;
            
            // Decay length
            GlobalVector decayVertexVector(Kpi.vx(),Kpi.vy(),Kpi.vz());
            GlobalVector primaryVertexVector(refittedPrimaryVertex.x(),refittedPrimaryVertex.y(),refittedPrimaryVertex.z());
            GlobalVector KpiVector(Kpi.px(),Kpi.py(),Kpi.pz());
            if ( KpiVector.dot(decayVertexVector-primaryVertexVector) < 0. ) continue; //
            Vertex decayVertex(Kpi.vertex(),Kpi.vertexCovariance(),Kpi.vertexChi2(),Kpi.vertexNdof(),Kpi.numberOfDaughters());
            GlobalVector nullVector(0.,0.,0.);
            SecondaryVertex sv(refittedPrimaryVertex,decayVertex,nullVector,true);
            if ( sv.dist3d().value() < KpiDecayLengthMin_ ) continue;
            if ( sv.dist3d().value()/sv.dist3d().error() < KpiDecayLengthSigMin_ ) continue;
            
            // Point/Distance of closest approach
            GlobalPoint primaryVertexPoint(refittedPrimaryVertex.x(),refittedPrimaryVertex.y(),refittedPrimaryVertex.z());
            GlobalPoint decayVertexPoint(Kpi.vx(),Kpi.vy(),Kpi.vz());
            Line lineOfFlight(decayVertexPoint,KpiVector);  // line of flight
            GlobalVector pocaVector = primaryVertexVector + lineOfFlight.distance(primaryVertexPoint);
            GlobalPoint pocaPoint(pocaVector.x(),pocaVector.y(),pocaVector.z());
            Point PoCA(KpiCand -> stateAtPoint(pocaPoint).globalPosition().x(),
                       KpiCand -> stateAtPoint(pocaPoint).globalPosition().y(),
                       KpiCand -> stateAtPoint(pocaPoint).globalPosition().z() );
            CovarianceMatrix PoCAError;
            PoCAError.At(0,0) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(0,0);
            PoCAError.At(1,0) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(1,0);
            PoCAError.At(1,1) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(1,1);
            PoCAError.At(2,0) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,0);
            PoCAError.At(2,1) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,1);
            PoCAError.At(2,2) = KpiCand -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,2);
            
            PoCAs.push_back(PoCA);
            PoCAErrors.push_back(PoCAError);
            RefittedPrimaryVertices.push_back(refittedPrimaryVertex);
            Kpis.push_back(Kpi);
         }
      }
      // _______________________________________________   
      
      // save collections in Event
      {
         // save original selected primary vertex in Event
         OriginalPrimaryVertices.push_back(bestPV);
         std::auto_ptr<VertexCollection> OriginalPVCollection(new VertexCollection);
         OriginalPVCollection -> reserve(OriginalPrimaryVertices.size());
         std::copy( OriginalPrimaryVertices.begin(), OriginalPrimaryVertices.end(), std::back_inserter(*OriginalPVCollection) );
         iEvent.put(OriginalPVCollection, "OriginalPrimaryVertex");
         
         // save refitted primary vertex in Event
         std::auto_ptr<VertexCollection> RefittedPVCollection(new VertexCollection);
         RefittedPVCollection -> reserve(RefittedPrimaryVertices.size());
         std::copy( RefittedPrimaryVertices.begin(), RefittedPrimaryVertices.end(), std::back_inserter(*RefittedPVCollection) );
         iEvent.put(RefittedPVCollection, "RefittedPrimaryVertex");
         
         // save PoCA in Event
         std::auto_ptr<PointCollection> PoCACollection(new PointCollection);
         PoCACollection -> reserve(PoCAs.size());
         std::copy( PoCAs.begin(), PoCAs.end(), std::back_inserter(*PoCACollection) );
         iEvent.put(PoCACollection, "PointOfClosestApproach");
         // save PoCA error in Event
         std::auto_ptr<CovarianceMatrixCollection> PoCAErrorCollection(new CovarianceMatrixCollection);
         PoCAErrorCollection -> reserve(PoCAErrors.size());
         std::copy( PoCAErrors.begin(), PoCAErrors.end(), std::back_inserter(*PoCAErrorCollection) );
         iEvent.put(PoCAErrorCollection, "PointOfClosestApproachError");
      
         // save K-pi candidates in Event
         std::auto_ptr<VertexCompositeCandidateCollection> KpiCollection(new VertexCompositeCandidateCollection);
         KpiCollection -> reserve(Kpis.size());
         std::copy( Kpis.begin(), Kpis.end(), std::back_inserter(*KpiCollection) );
         iEvent.put(KpiCollection);
      }
   }
   catch(...){
      std::cout << "Some of the input collections is absent in an event container " << std::endl;
   }
   
}

// ------------ method called once each job just before starting event loop  ------------
void KpiProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpiProducer::endJob()
{
}

bool KpiProducer::IsGoodTrack ( TrackRef track, Vertex vertex )
{   
   bool isGoodTrack = 
         fabs(track -> dxy(vertex.position())) < trkDxyMax_
      && fabs(track -> dz(vertex.position())) < trkDzMax_
      && track -> normalizedChi2() < trkNormChi2Max_
      && track -> eta() > trkEtaRange_[0] && track -> eta() < trkEtaRange_[1] 
      && track -> pt() > trkPtMin_
      && track -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
      && track -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_;
//   && track -> numberOfValidHits() >= trkPixelHitsMin_
//   && track -> hitPattern().trackerLayersWithMeasurement() >= trkSiliconHitsMin_;
   
   return isGoodTrack;
}

// ======= REFIT K-PI =======  
VertexCompositeCandidate KpiProducer::RefittedKpi
( const TrackRef & kaonTrackRef, const TrackRef & pionTrackRef, const MagneticField * bField, RefCountedKinematicParticle & KpiRef, bool & status )
{
   status = false;
   VertexCompositeCandidate dummy;
   
   TransientTrack kaonTT(kaonTrackRef, bField );
   TransientTrack pionTT(pionTrackRef, bField );
   
   //initial chi2 and ndof before kinematic fits.
   float chi2 = 0.;
   float ndof = 0.;
   
   //Creating a KinematicParticleFactory
   KinematicParticleFactoryFromTransientTrack pFactory;
   std::vector<RefCountedKinematicParticle> KpiParticles;
   KpiParticles.push_back(pFactory.particle(kaonTT,kaon_mass_,chi2,ndof,kaon_sigma_));
   KpiParticles.push_back(pFactory.particle(pionTT,pion_mass_,chi2,ndof,pion_sigma_));
   
   KinematicParticleVertexFitter fitter;   
   RefCountedKinematicTree vertexFitTree = fitter.fit(KpiParticles); 
   if (!vertexFitTree->isValid()) return dummy;
   
   // refitted Kaon
   vertexFitTree->movePointerToTheFirstChild();
   RefCountedKinematicParticle kaonCand = vertexFitTree -> currentParticle();
   KinematicParameters kaonCandKP = kaonCand -> currentState().kinematicParameters();
   // refitted pion(1)
   vertexFitTree->movePointerToTheNextChild();
   RefCountedKinematicParticle pionCand = vertexFitTree -> currentParticle();
   KinematicParameters pionCandKP = pionCand -> currentState().kinematicParameters();
   // refitted D0
   vertexFitTree -> movePointerToTheTop();
   RefCountedKinematicParticle kpiCand = vertexFitTree -> currentParticle();
   
   // Reconstructed D0 candidate
   Candidate::LorentzVector kpi4Momentum;
   double kpiEnergy;
   double kpiMass = kpiCand -> currentState().mass();
   int kpiCharge = kpiCand -> currentState().particleCharge();
   if ( kpiMass < KpiMassWindow_[0] || kpiMass > KpiMassWindow_[1] ) return dummy;
   kpi4Momentum.SetPx(kpiCand -> currentState().globalMomentum().x());
   kpi4Momentum.SetPy(kpiCand -> currentState().globalMomentum().y());
   kpi4Momentum.SetPz(kpiCand -> currentState().globalMomentum().z());
   kpiEnergy = sqrt(kpi4Momentum.P2() + kpiMass*kpiMass );
   kpi4Momentum.SetE (kpiEnergy);
   if ( kpi4Momentum.Pt() < KpiPtMin_ ) return dummy;
   if ( kpi4Momentum.Eta() < KpiEtaRange_[0] || kpi4Momentum.Eta() > KpiEtaRange_[1] ) return dummy; 
   
   // Refitted decay vertex
   RefCountedKinematicVertex kpiDecayVertex = vertexFitTree -> currentDecayVertex();
   float kpiDecayVertexChi2 = kpiDecayVertex -> chiSquared();
   float kpiDecayVertexNdof = kpiDecayVertex -> degreesOfFreedom();
   float kpiDecayVertexProb = ChiSquaredProbability(kpiDecayVertexChi2,kpiDecayVertexNdof);
   if ( kpiDecayVertexProb < KpiDecayVertexProbMin_ ) return dummy;
   // Vertex position & error
   Point kpiDecayVertexPosition;
   CovarianceMatrix kpiDecayVertexError;
   kpiDecayVertexPosition.SetX(kpiDecayVertex -> position().x());
   kpiDecayVertexPosition.SetY(kpiDecayVertex -> position().y());
   kpiDecayVertexPosition.SetZ(kpiDecayVertex -> position().z());
   //
   kpiDecayVertexError.At(0,0) = kpiDecayVertex -> error().cxx();
   kpiDecayVertexError.At(1,0) = kpiDecayVertex -> error().cyx();
   kpiDecayVertexError.At(1,1) = kpiDecayVertex -> error().cyy();
   kpiDecayVertexError.At(2,0) = kpiDecayVertex -> error().czx();
   kpiDecayVertexError.At(2,1) = kpiDecayVertex -> error().czy();
   kpiDecayVertexError.At(2,2) = kpiDecayVertex -> error().czz();
   
   // Reconstruct a kpi candidate            
   VertexCompositeCandidate kpi( kpiCharge, kpi4Momentum,
                                  kpiDecayVertexPosition, kpiDecayVertexError,
                                  kpiDecayVertexChi2, kpiDecayVertexNdof,
                                  KpiPDG_ );
   
   // Create daughter candidates for the VertexCompositeCandidates
   // Kaon candidates
   Candidate::LorentzVector kaon4Momentum;
   kaon4Momentum.SetPx(kaonCand -> currentState().globalMomentum().x());
   kaon4Momentum.SetPy(kaonCand -> currentState().globalMomentum().y());
   kaon4Momentum.SetPz(kaonCand -> currentState().globalMomentum().z());
   double kaonEnergy = sqrt((kaonCand -> currentState().globalMomentum().mag2()) +
                            (kaonCand -> currentState().mass() * kaonCand -> currentState().mass()));
   kaon4Momentum.SetE (kaonEnergy);
   int kaonCharge = kaonCand -> currentState().particleCharge();
   RecoChargedCandidate theKaonCand(kaonCharge, kaon4Momentum, kpiDecayVertexPosition);
   theKaonCand.setTrack(kaonTrackRef);
   // pion1 candidates
   Candidate::LorentzVector pion4Momentum;
   pion4Momentum.SetPx(pionCand -> currentState().globalMomentum().x());
   pion4Momentum.SetPy(pionCand -> currentState().globalMomentum().y());
   pion4Momentum.SetPz(pionCand -> currentState().globalMomentum().z());
   double pionEnergy = sqrt((pionCand -> currentState().globalMomentum().mag2()) +
                           (pionCand -> currentState().mass() * pionCand -> currentState().mass()));
   pion4Momentum.SetE (pionEnergy);
   int pionCharge = pionCand -> currentState().particleCharge();
   RecoChargedCandidate thePionCand(pionCharge, pion4Momentum, kpiDecayVertexPosition);
   thePionCand.setTrack(pionTrackRef);
   
   kpi.addDaughter(theKaonCand,"kaon");
   kpi.addDaughter(thePionCand,"pion");
   
   status = true;
   KpiRef = kpiCand;
   
   return kpi;
}


// ======= REFIT PRIMARY VERTEX =======            
Vertex KpiProducer::RefittedPrimaryVertex
( const Vertex & bestPV, const TrackRef & kaonTrackRef, const TrackRef & pionTrackRef,
 const MagneticField  * bField, const BeamSpot & beamSpot, bool & status )
{
   status = false;
   Vertex dummy;
   Vertex refittedPrimaryVertex;
   //   Check for K-pi tracks in primary vertex, remove them end refit primary vertex
   vector<reco::TransientTrack> newVertexTracks;
   float weightKaon = 0.;
   float weightPion = 0.;
   
   for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestPV.tracks_begin(); iTrack != bestPV.tracks_end(); ++iTrack)
   {
      // compare primary tracks to check for matches with Xi cand
      TrackRef trackRef = iTrack->castTo<TrackRef>();
      
      // the 3 tracks in the D+ candidate are theDaughterTracks[0] theDaughterTracks[1] glbTrack
      if ( kaonTrackRef == trackRef ) weightKaon = bestPV.trackWeight ( *iTrack );
      if ( pionTrackRef == trackRef ) weightPion = bestPV.trackWeight ( *iTrack );
      
      // try refitting the primary without the tracks in the D+ candidate	   
      if ( ! ( (kaonTrackRef == trackRef) || (pionTrackRef == trackRef) ) )
      {
         TransientTrack newVertexTT (trackRef, bField );
         newVertexTracks.push_back(newVertexTT);
      }
   }
   refittedPrimaryVertex = bestPV;
   
   // if no tracks in primary or no reco track included in primary then don't do anything
   if ( refitPrimaryVertex_ && newVertexTracks.size() > 0 && bestPV.tracksSize() != newVertexTracks.size() )
   {
      AdaptiveVertexFitter pvFitter;
      TransientVertex pvTrans = pvFitter.vertex( newVertexTracks, beamSpot );
      if ( pvTrans.isValid() ) refittedPrimaryVertex = pvTrans;
   }
   status = true;
   return refittedPrimaryVertex;
}

// 
bool KpiProducer:: IsGoodBtwTwoTracks( const TransientTrack* tt1, const TransientTrack* tt2)
{
   
   FreeTrajectoryState tt1State  = tt1 -> impactPointTSCP().theState();
   FreeTrajectoryState tt2State  = tt2  -> impactPointTSCP().theState();
   // Closest approach and crossing points of the tracks
   ClosestApproachInRPhi cApp;
   GlobalPoint cxPt;
   double dca;
   // track1- track2
   cApp.calculate(tt1State, tt2State);
   if ( ! cApp.status() ) return false;
   dca = fabs( cApp.distance() );
   if ( dca < 0. || dca > trkDCAMax_ ) return false;
   cxPt = cApp.crossingPoint();
   if ( cxPt.transverse() > 120. || fabs(cxPt.z()) > 300.) return false;
   TrajectoryStateClosestToPoint tt1TSCP = tt1 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt1TSCP.isValid() ) return false;
   TrajectoryStateClosestToPoint tt2TSCP = tt2 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt2TSCP.isValid() ) return false;
   
   return true;
   
}

//define this as a plug-in
DEFINE_FWK_MODULE(KpiProducer);

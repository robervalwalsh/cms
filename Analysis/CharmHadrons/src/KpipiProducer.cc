// -*- C++ -*-
//
// Package:    KpipiProducer
// Class:      KpipiProducer
// 
/**\class KpipiProducer KpipiProducer.cc CharmHadrons/Dplus/src/KpipiProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Roberval Walsh,01b/352,4849,
//         Created:  Wed Jan 27 10:40:28 CET 2010
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

class KpipiProducer : public edm::EDProducer
{
public:
   explicit KpipiProducer(const edm::ParameterSet&);
   ~KpipiProducer();
   
private:
   virtual void beginJob() ;
   virtual void produce(edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   
   bool IsGoodTrack ( TrackRef, Vertex );
   bool AreGoodTwoTracks( const TransientTrack*, const TransientTrack*);
   VertexCompositeCandidate RefittedKpipi( const TrackRef&, const TrackRef&, const TrackRef&, const MagneticField *, RefCountedKinematicParticle &, bool & );
   
   
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
   double DecayLengthSigMin_;
   double DecayLengthMin_;
   double DecayLengthSig2DMin_;
   double DecayLength2DMin_;
   double ProperDecayLengthMin_;
   double DecayLengthAngleMax_;
   
   int PDG_;
   
   ParticleMass pion_mass_;
   ParticleMass proton_mass_;
   ParticleMass kaon_mass_;
   ParticleMass lambda_mass_;
   float pion_sigma_ ;
   float proton_sigma_;
   float kaon_sigma_;
   float lambda_sigma_;
   
   
};

KpipiProducer::KpipiProducer(const edm::ParameterSet& iConfig)
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
   trkDCAMax_         = iConfig.getParameter<double> ("trkDCAMax");
   trkKaonPtMin_      = iConfig.getParameter<double> ("trkKaonPtMin");
   trkPionPtMin_      = iConfig.getParameter<double> ("trkPionPtMin");

   trkXpTMax_         = iConfig.getParameter<double> ("trkXpTMax");
   trkXpZMax_         = iConfig.getParameter<double> ("trkXpZMax");
   
   // D+ selection
   MassWindow_           = iConfig.getParameter< std::vector<double> > ("MassWindow");
   PtMin_                = iConfig.getParameter<double> ("PtMin");
   EtaRange_             = iConfig.getParameter< std::vector<double> > ("EtaRange");
   DecayVertexProbMin_   = iConfig.getParameter<double> ("DecayVertexProbMin");
   DecayLengthSigMin_    = iConfig.getParameter<double> ("DecayLengthSigMin");
   DecayLengthMin_       = iConfig.getParameter<double> ("DecayLengthMin");
   DecayLengthSig2DMin_  = iConfig.getParameter<double> ("DecayLengthSig2DMin");
   DecayLength2DMin_     = iConfig.getParameter<double> ("DecayLength2DMin");
   ProperDecayLengthMin_ = iConfig.getParameter<double> ("ProperDecayLengthMin");
   DecayLengthAngleMax_  = iConfig.getParameter<double> ("DecayLengthAngleMax");

   PDG_                  = iConfig.getParameter<int> ("PDG");

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
   produces<CovarianceMatrixCollection>("PointOfClosestApproachError");
   produces<PointCollection>("PointOfClosestApproach");
   
}

KpipiProducer::~KpipiProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
KpipiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   bool valid;
   
   ESHandle<MagneticField> bFieldHandler;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandler);
   
   // Primary vertices (primary vertex chosen in the config by PrimaryVertexSelector - vertex with largest #tracks)
   Handle<VertexCollection> pvHandler;
   iEvent.getByLabel(srcPVs_, pvHandler);
   Vertex bestPV = pvHandler -> at(0);

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
   
   std::vector<VertexCompositeCandidate> Kpipis;
   std::vector<Point> PoCAs;
   std::vector<CovarianceMatrix> PoCAErrors;
   
   for ( size_t i = 0; i < nTracks; ++i ) // kaon track
   {
      // Kaon candidate
      TrackRef kaonTrackRef( trackHandler, i );
      if ( ! IsGoodTrack(kaonTrackRef,bestPV) && kaonTrackRef -> pt() < trkKaonPtMin_  ) continue;
      // Transient tracks
      TransientTrack kaonTT(kaonTrackRef, &(*bFieldHandler) );
      TransientTrack* kaonTTPtr = &kaonTT;
      if ( ! kaonTTPtr -> impactPointTSCP().isValid() ) continue;
      // pion1 candidate
      for ( size_t j = 0; j < nTracks-1; ++j ) // pion(1) track
      {
         if ( j == i ) continue;  // same track
         TrackRef pi1TrackRef( trackHandler, j );
         if ( ! IsGoodTrack(pi1TrackRef,bestPV) && pi1TrackRef -> pt() < trkPionPtMin_ ) continue;
         if ( kaonTrackRef->charge()*pi1TrackRef->charge() > 0 ) continue;  // remove wrong charge combintations
         // Transient tracks
         TransientTrack pi1TT (pi1TrackRef, &(*bFieldHandler) );
         TransientTrack* pi1TTPtr  = &pi1TT;
         if ( ! pi1TTPtr -> impactPointTSCP().isValid() ) continue;
         // check status between 2 tracks, dca etc.
         if ( ! AreGoodTwoTracks(kaonTTPtr,pi1TTPtr) ) continue;
         
         // pion2 candidate
         for ( size_t k = j+1; k < nTracks; ++k ) // pion(2) track
         {
            if ( k == i ) continue; // same track
            TrackRef pi2TrackRef( trackHandler, k );
            if ( ! IsGoodTrack(pi2TrackRef,bestPV) && pi2TrackRef -> pt() < trkPionPtMin_) continue;
            if ( kaonTrackRef->charge()*pi2TrackRef->charge() > 0 ) continue;  // remove wrong charge combintations
            // Transient tracks
            TransientTrack pi2TT (pi2TrackRef, &(*bFieldHandler) );
            TransientTrack* pi2TTPtr  = &pi2TT;
            if ( ! pi2TTPtr -> impactPointTSCP().isValid() ) continue;
            // check status between 2 tracks, dca etc.
            if ( ! AreGoodTwoTracks(kaonTTPtr,pi2TTPtr) ) continue;
            if ( ! AreGoodTwoTracks(pi1TTPtr,pi2TTPtr) ) continue;
            
            // Reconstruct K-pi-pi candidate
            valid = true;
            RefCountedKinematicParticle KpipiKP;
            
            VertexCompositeCandidate Kpipi = RefittedKpipi(kaonTrackRef, pi1TrackRef, pi2TrackRef, &(*bFieldHandler), KpipiKP, valid);
            if ( ! valid ) continue;
            CovarianceMatrix DecayVertexError = Kpipi.vertexCovariance();
            
            // Decay length
            GlobalVector decayVertexVector(Kpipi.vx(),Kpipi.vy(),Kpipi.vz());
            GlobalVector primaryVertexVector(bestPV.x(),bestPV.y(),bestPV.z());
            GlobalVector dplusVector(Kpipi.px(),Kpipi.py(),Kpipi.pz());
            GlobalVector decayLengthVector = decayVertexVector - primaryVertexVector;

            Vertex decayVertex(Kpipi.vertex(),Kpipi.vertexCovariance(),Kpipi.vertexChi2(),Kpipi.vertexNdof(),Kpipi.numberOfDaughters());
            GlobalVector nullVector(0.,0.,0.);
            SecondaryVertex sv(bestPV,decayVertex,dplusVector,true);
            if ( sv.dist3d().value() < DecayLengthMin_ ) continue;
            if ( sv.dist3d().value()/sv.dist3d().error() < DecayLengthSigMin_ ) continue;
            if ( sv.dist2d().value() < DecayLength2DMin_ ) continue;
            if ( sv.dist2d().value()/sv.dist2d().error() < DecayLengthSig2DMin_ ) continue;
            // Proper time
            double properDecayLength = Kpipi.mass()/Kpipi.pt()*sv.dist2d().value();
            if ( properDecayLength < ProperDecayLengthMin_ ) continue;
            
            double DecayLengthAngle = acos(fabs(decayLengthVector.dot(dplusVector)/decayLengthVector.mag()/dplusVector.mag()));
            if ( DecayLengthAngle > DecayLengthAngleMax_ ) continue;

            // Point/Distance of closest approach
            GlobalPoint primaryVertexPoint(bestPV.x(),bestPV.y(),bestPV.z());
            GlobalPoint decayVertexPoint(Kpipi.vx(),Kpipi.vy(),Kpipi.vz());
            Line lineOfFlight(decayVertexPoint,dplusVector);  // line of flight
            GlobalVector pocaVector = primaryVertexVector + lineOfFlight.distance(primaryVertexPoint);
            GlobalPoint pocaPoint(pocaVector.x(),pocaVector.y(),pocaVector.z());
            Point PoCA(KpipiKP -> stateAtPoint(pocaPoint).globalPosition().x(),
                       KpipiKP -> stateAtPoint(pocaPoint).globalPosition().y(),
                       KpipiKP -> stateAtPoint(pocaPoint).globalPosition().z() );
            CovarianceMatrix PoCAError;
            PoCAError.At(0,0) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(0,0);
            PoCAError.At(1,0) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(1,0);
            PoCAError.At(1,1) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(1,1);
            PoCAError.At(2,0) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,0);
            PoCAError.At(2,1) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,1);
            PoCAError.At(2,2) = KpipiKP -> stateAtPoint(pocaPoint).kinematicParametersError().matrix()(2,2);
            
            PoCAs.push_back(PoCA);
            PoCAErrors.push_back(PoCAError);
            Kpipis.push_back(Kpipi);
         }
      }
   }

   // save collections in Event
   {
      
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
      
      // save Dplus candidates in Event
      std::auto_ptr<VertexCompositeCandidateCollection> KpipiCollection(new VertexCompositeCandidateCollection);
      KpipiCollection -> reserve(Kpipis.size());
      std::copy( Kpipis.begin(), Kpipis.end(), std::back_inserter(*KpipiCollection) );
      iEvent.put(KpipiCollection);
   }
   
}

// ------------ method called once each job just before starting event loop  ------------
void KpipiProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void KpipiProducer::endJob()
{
}


// ======= REFIT K-PI-PI =======  
VertexCompositeCandidate KpipiProducer::RefittedKpipi
( const TrackRef & kaonTrackRef, const TrackRef & pi1TrackRef, const TrackRef & pi2TrackRef,
  const MagneticField * bField, RefCountedKinematicParticle & KpipiRef, bool & status )
{
   status = false;
   VertexCompositeCandidate dummy;

   TransientTrack kaonTT(kaonTrackRef, bField );
   TransientTrack pi1TT (pi1TrackRef,  bField );
   TransientTrack pi2TT (pi2TrackRef,  bField );
   
   //initial chi2 and ndof before kinematic fits.
   float chi2 = 0.;
   float ndof = 0.;
   
   //Creating a KinematicParticleFactory
   KinematicParticleFactoryFromTransientTrack pFactory;
   std::vector<RefCountedKinematicParticle> KpipiParticles;
   KpipiParticles.push_back(pFactory.particle(kaonTT,kaon_mass_,chi2,ndof,kaon_sigma_));
   KpipiParticles.push_back(pFactory.particle(pi1TT, pion_mass_,chi2,ndof,pion_sigma_));
   KpipiParticles.push_back(pFactory.particle(pi2TT ,pion_mass_,chi2,ndof,pion_sigma_));
   
   KinematicParticleVertexFitter fitter;   
   RefCountedKinematicTree vertexFitTree = fitter.fit(KpipiParticles); 
   if (!vertexFitTree->isValid()) return dummy;
   
   // refitted Kaon
   vertexFitTree->movePointerToTheFirstChild();
   RefCountedKinematicParticle kaonCand = vertexFitTree -> currentParticle();
   KinematicParameters kaonCandKP = kaonCand -> currentState().kinematicParameters();
   // refitted pion(1)
   vertexFitTree->movePointerToTheNextChild();
   RefCountedKinematicParticle pi1Cand = vertexFitTree -> currentParticle();
   KinematicParameters pi1CandKP = pi1Cand -> currentState().kinematicParameters();
   // refitted pion(2)
   vertexFitTree->movePointerToTheNextChild();
   RefCountedKinematicParticle pi2Cand = vertexFitTree -> currentParticle();
   KinematicParameters pi2CandKP = pi2Cand -> currentState().kinematicParameters();
   // refitted D+
   vertexFitTree -> movePointerToTheTop();
   RefCountedKinematicParticle KpipiKP = vertexFitTree -> currentParticle();
   
   // Reconstructed D+ candidate
   Candidate::LorentzVector FourMomentum;
   double Energy;
   double Mass = KpipiKP -> currentState().mass();
   if ( Mass < MassWindow_[0] || Mass > MassWindow_[1] ) return dummy;
   FourMomentum.SetPx(KpipiKP -> currentState().globalMomentum().x());
   FourMomentum.SetPy(KpipiKP -> currentState().globalMomentum().y());
   FourMomentum.SetPz(KpipiKP -> currentState().globalMomentum().z());
   Energy = sqrt(FourMomentum.P2() + Mass*Mass );
   FourMomentum.SetE (Energy);
   if ( FourMomentum.Pt() < PtMin_ ) return dummy;
   if ( FourMomentum.Eta() < EtaRange_[0] || FourMomentum.Eta() > EtaRange_[1] ) return dummy; 
   int Charge  = KpipiKP -> currentState().particleCharge();
   
   // Refitted decay vertex
   RefCountedKinematicVertex DecayVertex = vertexFitTree -> currentDecayVertex();
   float DecayVertexChi2 = DecayVertex -> chiSquared();
   float DecayVertexNdof = DecayVertex -> degreesOfFreedom();
   float DecayVertexProb = ChiSquaredProbability(DecayVertexChi2,DecayVertexNdof);
   if ( DecayVertexProb < DecayVertexProbMin_ ) return dummy;
   // Vertex position & error
   Point DecayVertexPosition;
   CovarianceMatrix DecayVertexError;
   DecayVertexPosition.SetX(DecayVertex -> position().x());
   DecayVertexPosition.SetY(DecayVertex -> position().y());
   DecayVertexPosition.SetZ(DecayVertex -> position().z());
   //
   DecayVertexError.At(0,0) = DecayVertex -> error().cxx();
   DecayVertexError.At(1,0) = DecayVertex -> error().cyx();
   DecayVertexError.At(1,1) = DecayVertex -> error().cyy();
   DecayVertexError.At(2,0) = DecayVertex -> error().czx();
   DecayVertexError.At(2,1) = DecayVertex -> error().czy();
   DecayVertexError.At(2,2) = DecayVertex -> error().czz();
   
   // Reconstruct a Dplus candidate            
   VertexCompositeCandidate Kpipi( Charge, FourMomentum,
                                   DecayVertexPosition, DecayVertexError,
                                   DecayVertexChi2, DecayVertexNdof,
                                   PDG_*Charge );
   
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
   RecoChargedCandidate theKaonCand(kaonCharge, kaon4Momentum, DecayVertexPosition);
   theKaonCand.setTrack(kaonTrackRef);
   // pion1 candidates
   Candidate::LorentzVector pi14Momentum;
   pi14Momentum.SetPx(pi1Cand -> currentState().globalMomentum().x());
   pi14Momentum.SetPy(pi1Cand -> currentState().globalMomentum().y());
   pi14Momentum.SetPz(pi1Cand -> currentState().globalMomentum().z());
   double pi1Energy = sqrt((pi1Cand -> currentState().globalMomentum().mag2()) +
                           (pi1Cand -> currentState().mass() * pi1Cand -> currentState().mass()));
   pi14Momentum.SetE (pi1Energy);
   int pi1Charge = pi1Cand -> currentState().particleCharge();
   RecoChargedCandidate thePi1Cand(pi1Charge, pi14Momentum, DecayVertexPosition);
   thePi1Cand.setTrack(pi1TrackRef);
   // pion2 candidates
   Candidate::LorentzVector pi24Momentum;
   pi24Momentum.SetPx(pi2Cand -> currentState().globalMomentum().x());
   pi24Momentum.SetPy(pi2Cand -> currentState().globalMomentum().y());
   pi24Momentum.SetPz(pi2Cand -> currentState().globalMomentum().z());
   double pi2Energy = sqrt((pi2Cand -> currentState().globalMomentum().mag2()) +
                           (pi2Cand -> currentState().mass() * pi2Cand -> currentState().mass()));
   pi24Momentum.SetE (pi2Energy);
   int pi2Charge = pi2Cand -> currentState().particleCharge();
   RecoChargedCandidate thePi2Cand(pi2Charge, pi24Momentum, DecayVertexPosition);
   thePi2Cand.setTrack(pi2TrackRef);
   
   Kpipi.addDaughter(theKaonCand,"kaon");
   Kpipi.addDaughter(thePi1Cand,"pion1");
   Kpipi.addDaughter(thePi2Cand,"pion2");
   
   status = true;
   KpipiRef = KpipiKP;
   
   return Kpipi;
}

// ======= Track quality cuts =======  
bool KpipiProducer::IsGoodTrack ( TrackRef track, Vertex vertex )
{   
   bool isGoodTrack = 
   fabs(track -> dxy(vertex.position())) < trkDxyMax_
   && fabs(track -> dz(vertex.position())) < trkDzMax_
   && track -> normalizedChi2() < trkNormChi2Max_
   && track -> eta() > trkEtaRange_[0] && track -> eta() < trkEtaRange_[1] 
   && track -> hitPattern().numberOfValidPixelHits() >= trkPixelHitsMin_
   && track -> hitPattern().numberOfValidTrackerHits() >= trkSiliconHitsMin_;
   
   return isGoodTrack;
}


bool KpipiProducer:: AreGoodTwoTracks( const TransientTrack* tt1, const TransientTrack* tt2)
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
//   if ( cxPt.transverse() > 120. || fabs(cxPt.z()) > 300.) return false;
   if ( cxPt.transverse() > trkXpTMax_ || fabs(cxPt.z()) > trkXpZMax_ ) return false;
   TrajectoryStateClosestToPoint tt1TSCP = tt1 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt1TSCP.isValid() ) return false;
   TrajectoryStateClosestToPoint tt2TSCP = tt2 -> trajectoryStateClosestToPoint( cxPt );
   if ( ! tt2TSCP.isValid() ) return false;
   
   return true;
   
}

//define this as a plug-in
DEFINE_FWK_MODULE(KpipiProducer);

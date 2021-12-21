/*
 * DataDumper.cc
 *
 *  Created on: Nov 5, 2021
 *      Author: kbunkow
 */

//The includes in the PreTrackMatchedMuona and other places have no the proper includes of the below files,
//so they must be before  DataDumper.h
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "L1Trigger/Phase2L1GMT/interface/Constants.h"

#include "L1Trigger/Phase2L1GMT/interface/DataDumper.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

namespace Phase2L1GMT {

DataDumper::DataDumper(const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& ttTrackMCTruthToken,
    const edm::EDGetTokenT< std::vector< TrackingParticle > >& trackingParticleToken,
    std::string rootFileName) :
    ttTrackMCTruthToken(ttTrackMCTruthToken), trackingParticleToken(trackingParticleToken)
{
  initializeTTree(rootFileName);
}

DataDumper::~DataDumper() {
  rootFile->Write();

  delete rootTree;
  delete rootFile;
}


void DataDumper::initializeTTree(std::string rootFileName) {
  rootFile = new TFile(rootFileName.c_str(), "RECREATE");
  rootTree = new TTree("GmtMuonTree", "");

  rootTree->Branch("tpPt", &record.tpPt);
  rootTree->Branch("tpEta", &record.tpEta);
  rootTree->Branch("tpPhi", &record.tpPhi);
  rootTree->Branch("tpBeta", &record.tpBeta);

  //rootTree->Branch("tpCharge", &record.tpCharge); charge is in type
  rootTree->Branch("tpType", &record.tpType);

  rootTree->Branch("matching", &record.matching);

  rootTree->Branch("tttCurvature", &record.tttCurvature);
  rootTree->Branch("tttCharge", &record.tttCharge);
  rootTree->Branch("tttPt", &record.tttPt);
  rootTree->Branch("tttEta", &record.tttEta);

  rootTree->Branch("tttPhi", &record.tttPhi);
  rootTree->Branch("tttZ0", &record.tttZ0);
  rootTree->Branch("tttD0", &record.tttD0);

  rootTree->Branch("gmtBeta", &record.gmtBeta);

  rootTree->Branch("isGlobal", &record.isGlobal);

  rootTree->Branch("quality", &record.quality);

  rootTree->Branch("deltaCoords1", &record.deltaCoords1);
  rootTree->Branch("deltaCoords2", &record.deltaCoords2);

  rootTree->Branch("deltaEtas1", &record.deltaEtas1);
  rootTree->Branch("deltaEtas2", &record.deltaEtas2);


  rootTree->Branch("stubTiming", &record.stubTiming);

  rootTree->Branch("stubType", &record.stubType);

}

void DataDumper::getHandles(const edm::Event& event) {
  // MC truth association maps


  event.getByToken(ttTrackMCTruthToken, mcTruthTTTrackHandle);

  // tracking particles
  event.getByToken(trackingParticleToken, trackingParticleHandle);

  //LogTrace("gmtDataDumper")<<"DataDumper::getHandles ";
  muonTrackingParticlesFilled = false;
  muonTrackingParticles.clear();
}

void DataDumper::process(PreTrackMatchedMuon& preTrackMatchedMuon) {
  auto& ttTrackPtr = preTrackMatchedMuon.trkPtr();

  if(ttTrackPtr.isNull())
    return;

  //from ttTrack
  record.tttCurvature = preTrackMatchedMuon.curvature();
  record.tttCharge = preTrackMatchedMuon.charge();
  record.tttPt = preTrackMatchedMuon.pt();
  record.tttEta = preTrackMatchedMuon.eta();
  record.tttPhi = preTrackMatchedMuon.phi();
  record.tttZ0 = preTrackMatchedMuon.z0();
  record.tttD0 = preTrackMatchedMuon.d0();

  LogTrace("gmtDataDumper")<<"DataDumper::process(): preTrackMatchedMuon pt: "<<record.tttPt<<" eta "<<record.tttEta<<" phi "<<record.tttPhi;
  //from GMT
  record.gmtBeta = preTrackMatchedMuon.beta();
  record.isGlobal =  preTrackMatchedMuon.isGlobalMuon();
  record.quality = preTrackMatchedMuon.quality();

  edm::Ptr< TrackingParticle > tpMatchedToL1MuCand = mcTruthTTTrackHandle->findTrackingParticlePtr(ttTrackPtr);

  if(tpMatchedToL1MuCand.isNonnull() ) {
    LogTrace("gmtDataDumper")<<" findTrackingParticlePtr() - found matching TrackingParticle";

    //something not good here, crashing
/*    if(mcTruthTTTrackHandle->isGenuine(ttTrackPtr))
      record.matching = 3;
    else if(mcTruthTTTrackHandle->isLooselyGenuine(ttTrackPtr))
      record.matching = 2;*/

    record.matching = 2;
  }
  else {
    LogTrace("gmtDataDumper")<<" findTrackingParticlePtr() - nothing found";
    if(!muonTrackingParticlesFilled) {
      for (unsigned int iTP = 0; iTP < trackingParticleHandle->size(); ++iTP) {
        edm::Ptr< TrackingParticle > tpPtr(trackingParticleHandle, iTP);
        if(abs(tpPtr->pdgId()) == 13 || abs(tpPtr->pdgId()) == 1000015) {
          muonTrackingParticles.push_back(tpPtr);
        }
      }
      LogTrace("gmtDataDumper")<<"filling muonTrackingParticles: muonTrackingParticles.size() = "<<muonTrackingParticles.size();
      muonTrackingParticlesFilled = true;
    }

    bool isVeryLoose = false;
    for(auto& muonTrackingPart : muonTrackingParticles) {
      //here we have ttTracks tagged as muon by correlator that have no matching genuine/loose genuine tracking particle
      //so we go over all muonTrackingParticles and check if muonTrackingParticle has given ttTrack matched,
      //here, "match" means ttTracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters - so it is very loose match
      std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = mcTruthTTTrackHandle->findTTTrackPtrs(muonTrackingPart);
      for(auto& matchedTTTrack : matchedTracks) {
        if(matchedTTTrack == ttTrackPtr) {
          isVeryLoose = true;

          tpMatchedToL1MuCand = muonTrackingPart;
          LogTrace("l1tMuBayesEventPrint") <<" veryLoose matching muonTrackingPart found";
          break;
        }
      }
      if(isVeryLoose) {
        record.matching = 1;
        break;
      }
    }
  }

  if(tpMatchedToL1MuCand.isNonnull() ) {
    //if(abs(tpMatchedToL1MuCand->pdgId()) == 13 || abs(tpMatchedToL1MuCand->pdgId()) == 1000015) {
    record.tpType = tpMatchedToL1MuCand->pdgId();
    record.tpPt = tpMatchedToL1MuCand->pt();
    record.tpEta = tpMatchedToL1MuCand->eta();
    record.tpPhi = tpMatchedToL1MuCand->phi();
    record.tpBeta = tpMatchedToL1MuCand->p4().Beta();

    LogTrace("gmtDataDumper")<<"ttTrack matched to the TrackingParticle";
    LogTrace("gmtDataDumper")<<" TrackingParticle type"<<(int)record.tpType<<" tpPt "<<record.tpPt<<" tpEta "<<record.tpEta<<" tpPhi "<<record.tpPhi ;

  }
  else {
    record.tpType = 0;
    record.tpPt = 0;
    record.tpEta = 0;
    record.tpPhi = 0;

    record.matching = 0;
  }


  LogTrace("gmtDataDumper")<<"gmtDataDumper preTrackMatchedMuon.getDeltaCoords1().size "<<preTrackMatchedMuon.getDeltaCoords1().size()
      <<" preTrackMatchedMuon.getDeltaCoords1().size "<<preTrackMatchedMuon.getDeltaCoords1().size();
  for (const auto& stub : preTrackMatchedMuon.stubs()) {
    LogTrace("gmtDataDumper")<<"gmtDataDumper stub: tfLayer "<<stub->tfLayer()<<" time() "<<stub->time()<<" type " <<stub->type();
    record.stubTiming.at(stub->tfLayer()) = stub->time();
    record.stubType.at(stub->tfLayer()) = stub->type();

  }

  for(unsigned int iLayer = 0; iLayer < preTrackMatchedMuon.getDeltaCoords1().size(); iLayer++) {
    record.deltaCoords1.at(iLayer) = preTrackMatchedMuon.getDeltaCoords1().at(iLayer).toRaw();
    record.deltaCoords2.at(iLayer) = preTrackMatchedMuon.getDeltaCoords2().at(iLayer).toRaw();

    record.deltaEtas1.at(iLayer) = preTrackMatchedMuon.getDeltaEtas1().at(iLayer).toRaw();
    record.deltaEtas2.at(iLayer) = preTrackMatchedMuon.getDeltaEtas2().at(iLayer).toRaw();

    LogTrace("gmtDataDumper")<<"gmtDataDumper record: tfLayer "<<iLayer
        <<" deltaCoords1 "<<(int)record.deltaCoords1.at(iLayer)
        <<" deltaCoords2 "<<(int)record.deltaCoords2.at(iLayer)
        <<" deltaEtas1 "<<(int)record.deltaEtas1.at(iLayer)
        <<" deltaEtas2 "<<(int)record.deltaEtas2.at(iLayer)
        <<" stubTiming "<<(int)record.stubTiming.at(iLayer)
        <<" stubType "<<(int)record.stubType.at(iLayer) ;
  }

  rootTree->Fill();
  evntCnt++;
}

} /* namespace Phase2L1GMT */

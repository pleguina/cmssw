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

#include "DataDumper.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include <iostream>
#include <fstream>


namespace Phase2L1GMT {

DataDumper::DataDumper(const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& ttTrackMCTruthToken,
    const edm::EDGetTokenT< std::vector< TrackingParticle > >& trackingParticleToken,
	bool dumpToRoot, bool dumpToXml) :
    ttTrackMCTruthToken(ttTrackMCTruthToken), trackingParticleToken(trackingParticleToken), dumpToRoot(dumpToRoot), dumpToXml(dumpToXml)
{
  if(dumpToRoot)
	  initializeTTree();
}

DataDumper::~DataDumper() {
  //rootFile->Write();

  //delete rootTree;
  //delete rootFile;
}


void DataDumper::initializeTTree() {

  edm::Service<TFileService> fs;

  TFileDirectory subDir = fs->mkdir("GmtDataDumper");

  rootTree = subDir.make<TTree>("GmtMuonTree", "");

  //rootTree = new TTree("GmtMuonTree", "");

  rootTree->Branch("eventNum", &record.eventNum);
  rootTree->Branch("tpEvent", &record.tpEvent);
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

  rootTree->Branch("tttChi2rphi", &record.tttChi2rphi);
  rootTree->Branch("tttChi2rz", &record.tttChi2rz);
  rootTree->Branch("tttBendChi2", &record.tttBendChi2);
  //rootTree->Branch("tttQualityMVA", &record.tttQualityMVA); lokks that is not set in the tracking trigger
  //rootTree->Branch("tttOtherMVA", &record.tttOtherMVA); is not set in the tracking trigger


  rootTree->Branch("gmtBeta", &record.gmtBeta);

  rootTree->Branch("isGlobal", &record.isGlobal);

  rootTree->Branch("quality", &record.quality);


  rootTree->Branch("hitsValid", &record.hitsValid);

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

/*
 * TODO
 * this DataDumper does not collect the muonTrackingPart the were not matched to the ttTrack
 * - so it is not possible to calculate the tracking trigger efficiency from the collected data
 * also it does not collect the muon ttTracks that were not tagged by the GMT as muons, a
 * - so it is not possible to calculate the GMT efficiency.
 * In principle this can be achieved but would require rewriting this DataDumper
 */
void DataDumper::process(PreTrackMatchedMuon& preTrackMatchedMuon) {
  if(!dumpToRoot)
    return;

  auto& ttTrackPtr = preTrackMatchedMuon.trkPtr();

  if(ttTrackPtr.isNull())
    return;

  //from ttTrack
  record.tttCurvature = preTrackMatchedMuon.curvature();
  record.tttCharge = preTrackMatchedMuon.charge();
  record.tttPt = preTrackMatchedMuon.pt(); //loks that this pt is pt in GeV i.e. ttTrackPtrpt momentum().perp() * 32
  record.tttEta = preTrackMatchedMuon.eta();
  record.tttPhi = preTrackMatchedMuon.phi();
  record.tttZ0 = preTrackMatchedMuon.z0();
  record.tttD0 = preTrackMatchedMuon.d0();

  LogTrace("gmtDataDumper")<<"DataDumper::process(): preTrackMatchedMuon pt: "<<record.tttPt<<" eta "<<record.tttEta<<" phi "<<record.tttPhi;
  LogTrace("gmtDataDumper")<<"ttTrackPtrpt momentum().perp() "<<ttTrackPtr->momentum().perp();

  //from GMT
  record.gmtBeta = preTrackMatchedMuon.beta();
  record.isGlobal =  preTrackMatchedMuon.isGlobalMuon();
  record.quality = preTrackMatchedMuon.quality();

  record.tttChi2rphi = ttTrackPtr->getChi2RPhiBits();
  record.tttChi2rz = ttTrackPtr->getChi2RZBits();
  record.tttBendChi2 = ttTrackPtr->getBendChi2Bits();
  record.tttQualityMVA = ttTrackPtr->getMVAQualityBits();
  record.tttOtherMVA = ttTrackPtr->getMVAOtherBits();

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
      //here we have ttTracks tagged as muon by GMT that have no matching genuine/loose genuine tracking particle
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

  record.eventNum = evntCnt;

  if(tpMatchedToL1MuCand.isNonnull() ) {
    record.tpEvent = tpMatchedToL1MuCand->eventId().event();
    //if(abs(tpMatchedToL1MuCand->pdgId()) == 13 || abs(tpMatchedToL1MuCand->pdgId()) == 1000015) {
    record.tpType = tpMatchedToL1MuCand->pdgId();
    record.tpPt = tpMatchedToL1MuCand->pt();
    record.tpEta = tpMatchedToL1MuCand->eta();
    record.tpPhi = tpMatchedToL1MuCand->phi();
    record.tpBeta = tpMatchedToL1MuCand->p4().Beta();

    LogTrace("gmtDataDumper")<<"ttTrack matched to the TrackingParticle";
    LogTrace("gmtDataDumper")<<" TrackingParticle event "<<record.tpEvent<<" type "<<record.tpType<<" tpPt "<<record.tpPt<<" tpEta "<<record.tpEta<<" tpPhi "<<record.tpPhi ;

  }
  else {
    record.tpType = 0;
    record.tpPt = 0;
    record.tpEta = 0;
    record.tpPhi = 0;

    record.matching = 0;
  }

  record.hitsValid =  preTrackMatchedMuon.getHitsValid();

  LogTrace("gmtDataDumper")<<"gmtDataDumper preTrackMatchedMuon.getDeltaCoords1().size "<<preTrackMatchedMuon.getDeltaCoords1().size()
      <<" preTrackMatchedMuon.getDeltaCoords1().size "<<preTrackMatchedMuon.getDeltaCoords1().size();

  for(unsigned int iLayer = 0; iLayer < preTrackMatchedMuon.getDeltaCoords1().size(); iLayer++) {
    record.deltaCoords1.at(iLayer) = preTrackMatchedMuon.getDeltaCoords1().at(iLayer).toRaw();
    record.deltaCoords2.at(iLayer) = preTrackMatchedMuon.getDeltaCoords2().at(iLayer).toRaw();

    record.deltaEtas1.at(iLayer) = preTrackMatchedMuon.getDeltaEtas1().at(iLayer).toRaw();
    record.deltaEtas2.at(iLayer) = preTrackMatchedMuon.getDeltaEtas2().at(iLayer).toRaw();

    //reseting here as in the next loop only these for which stub exists are set
    record.stubTiming.at(iLayer) = std::numeric_limits<decltype(record.stubTiming)::value_type>::min();
    record.stubType.at(iLayer) = 0xff; //-1
  }

  for (const auto& stub : preTrackMatchedMuon.stubs()) {
    LogTrace("gmtDataDumper")<<"gmtDataDumper stub: tfLayer "<<stub->tfLayer()<<" time() "<<stub->time()<<" type " <<stub->type();
    record.stubTiming.at(stub->tfLayer()) = stub->time();
    record.stubType.at(stub->tfLayer()) = stub->type();
  }

  //debug printout
  for(unsigned int iLayer = 0; iLayer < preTrackMatchedMuon.getDeltaCoords1().size(); iLayer++) {
    LogTrace("gmtDataDumper")<<"gmtDataDumper record: tfLayer "<<iLayer
        <<" deltaCoords1 "<<(int)record.deltaCoords1.at(iLayer)
        <<" deltaCoords2 "<<(int)record.deltaCoords2.at(iLayer)
        <<" deltaEtas1 "<<(int)record.deltaEtas1.at(iLayer)
        <<" deltaEtas2 "<<(int)record.deltaEtas2.at(iLayer)
        <<" stubTiming "<<(int)record.stubTiming.at(iLayer)
        <<" stubType "<<(int)record.stubType.at(iLayer) ;
  }

  rootTree->Fill();
  //evntCnt++;
}

void DataDumper::collectNonant(int nonant, std::vector<MuonROI>& rois,
    std::vector<ConvertedTTTrack>& convertedTracks,
    std::vector<PreTrackMatchedMuon>& muons,
    std::vector<PreTrackMatchedMuon>& muCleaned) {
  if(!dumpToXml)
    return;

  eventXml.nonants.at(nonant).rois.clear();

  for(auto& roi : rois) {
    auto stubInRoi =  std::vector<l1t::MuonStub>();
    for(auto& stub : roi.stubs()) {
      stubInRoi.push_back(*stub);
      eventXml.empty = false;
    }
    eventXml.nonants.at(nonant).rois.push_back(stubInRoi);
  }

  eventXml.nonants[nonant].convertedTracks = convertedTracks;
  eventXml.nonants[nonant].muons = muons;
  eventXml.nonants[nonant].muCleaned = muCleaned;
}

void DataDumper::writeToXml() {
  if(!dumpToXml)
    return;

  if(eventXml.empty == false) {
    // create and open a character archive for output
    std::ofstream ofs("event" + std::to_string(evntCnt) + ".xml");
    // save data to archive
    {
      boost::archive::xml_oarchive oa(ofs);
      // write class instance to archive
      oa << BOOST_SERIALIZATION_NVP(eventXml);
      // archive and stream closed when destructors are called
    }
  }
}

} /* namespace Phase2L1GMT */

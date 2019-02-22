/*
 * MuCorrelatorProcessor.cpp
 *
 *  Created on: Jan 18, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/MuCorrelatorProcessor.h"
#include "L1Trigger/L1TMuonBayes/interface/MuCorrelator/PdfModuleWithStats.h"

MuCorrelatorProcessor::MuCorrelatorProcessor(MuCorrelatorConfigPtr& config, std::string pdfModuleType): config(config) {
  if(pdfModuleType== "PdfModuleWithStats")
    pdfModule = std::make_unique<PdfModuleWithStats>(config);
  else
    pdfModule = std::make_unique<PdfModule>(config);

  ghostBustFunc = ghostBust3;
}

MuCorrelatorProcessor::MuCorrelatorProcessor(MuCorrelatorConfigPtr& config, unique_ptr<IPdfModule> pdfModule): config(config), pdfModule(std::move(pdfModule) ) {
  ghostBustFunc = ghostBust3;
}

MuCorrelatorProcessor::~MuCorrelatorProcessor() {
  // TODO Auto-generated destructor stub
}


AlgoTTMuons MuCorrelatorProcessor::processTracks(const MuonStubsInput& muonStubs, const TrackingTriggerTracks& ttTracks) {
  AlgoTTMuons algoTTMuons;

  for(auto& ttTrack : ttTracks) {
    //TODO add switch to use either processTrack or processTrackUsingRefStubs

    //cout<<"\n"<<__FUNCTION__<<":"<<__LINE__<<" "<<*ttTrack<<std::endl;
    auto algoTTMuon = processTrack(muonStubs, ttTrack);

    if(algoTTMuon->getFiredLayerCnt() >= config->nMinFiredLayers()) {
      algoTTMuon->setValid(true);
      algoTTMuons.emplace_back(algoTTMuon);
      //cout<<"\n"<<__FUNCTION__<<":"<<__LINE__<<" "<<*algoTTMuon<<endl;
    }
    //cout<<"\n";
  }

  auto ghostBustedTTmuons = ghostBust(algoTTMuons);

  /*
  for(auto& algoTTMuon : ghostBustedTTmuons) {
    cout<<__FUNCTION__<<":"<<__LINE__<<" ghostBusted "<<*algoTTMuon<<endl;
  }*/

  return ghostBustedTTmuons;
}


AlgoTTMuonPtr MuCorrelatorProcessor::processTrack(const MuonStubsInput& muonStubs, const TrackingTriggerTrackPtr& ttTrack) {
  //Selecting stubs that fit coarsely to the ttTrack, e.g. the full chambers
  MuonStubsInput selectedMuonStubs = selectStubs(muonStubs, ttTrack);

  AlgoTTMuonPtr algoTTMuon = std::make_shared<AlgoTTMuon>(ttTrack, config);
  for(unsigned int iLayer = 0; iLayer < config->nLayers(); ++iLayer) {
      processStubs(selectedMuonStubs, iLayer, ttTrack, MuonStubPtr(), algoTTMuon);
      //the muonStubs has no stubs in the banding layer, the phiB is processed wheh the corresponding phi layer is processed
  }

  return algoTTMuon;
}


AlgoTTMuonPtr MuCorrelatorProcessor::processTrackUsingRefStubs(const MuonStubsInput& muonStubs, const TrackingTriggerTrackPtr& ttTrack) {
  AlgoTTMuonPtr bestAlgoTTMuon;

  MuonStubsInput selectedMuonStubs = selectStubs(muonStubs, ttTrack);

  MuonStubPtrs1D refStubs = selectRefStubs(selectedMuonStubs, ttTrack);

  for(unsigned int iRefStub = 0; iRefStub < refStubs.size(); ++iRefStub) {

    AlgoTTMuonPtr algoTTMuon = std::make_shared<AlgoTTMuon>(ttTrack, config, refStubs[iRefStub]);
    for(unsigned int iLayer = 0; iLayer < config->nLayers(); ++iLayer) {
      processStubs(selectedMuonStubs, iLayer, ttTrack, refStubs[iRefStub], algoTTMuon);
    }

    //TODO do something better?
    if(algoTTMuon->getFiredLayerCnt() >= config->nMinFiredLayers()) {
      algoTTMuon->setValid(true);
    }

    if(algoTTMuon->isValid()) {
      if( !bestAlgoTTMuon || //TODO maybe better just use the pdfSum - check and optimize
          (algoTTMuon->getFiredLayerCnt() >  bestAlgoTTMuon->getFiredLayerCnt() ) ||
          (algoTTMuon->getFiredLayerCnt() == bestAlgoTTMuon->getFiredLayerCnt() && algoTTMuon->getPdfSum() > bestAlgoTTMuon->getPdfSum() )
      )
      {
        bestAlgoTTMuon = algoTTMuon;
      }
    }

  }

  return bestAlgoTTMuon;
}

MuonStubsInput MuCorrelatorProcessor::selectStubs(const MuonStubsInput& muonStubs, const TrackingTriggerTrackPtr& ttTrack) {
  MuonStubsInput selectedMuonStubs(config);

  //TODO this implementation is rather not possible in the hardware, a different approach would be needed
  for(unsigned int iLayer = 0; iLayer < muonStubs.getMuonStubs().size(); ++iLayer) {
    for(auto& stub : muonStubs.getMuonStubs()[iLayer] ) {
      if(config->isPhiLayer(iLayer)) {//phi stubs, TODO add a better way to distinguish eta and phi stubs
        //for phi stubs, we simply check that the eta extent of the chamber is more or less compatible with the eta of the ttTrack
        //TODO for the moment any phi extrapolation is ignored, implement it (it should be something simple at this stage, i.e. e.g. not including the pt
        int etaMargin = stub->etaSigmaHw + 10; //adding some additional margin, TODO optimize margin and move to config
        if( abs(stub->etaHw - ttTrack->getEtaHw()) < etaMargin ) {
          selectedMuonStubs.addStub(iLayer, stub);
        }
      }
      else {
        //TODO implement something for the eta stubs
      }
    }

  }

  return selectedMuonStubs;
}

MuonStubPtrs1D MuCorrelatorProcessor::selectRefStubs(const MuonStubsInput& muonStubs, const TrackingTriggerTrackPtr& ttTrack) {
  MuonStubPtrs1D refStubs;
  //TODO implement
  return refStubs;
}

void MuCorrelatorProcessor::processStubs(const MuonStubsInput& muonStubs, unsigned int layer, const TrackingTriggerTrackPtr& ttTrack, const MuonStubPtr refStub, AlgoTTMuonPtr algoTTMuon) {
  pdfModule->processStubs(muonStubs, layer, ttTrack, refStub, algoTTMuon);
}


AlgoTTMuons MuCorrelatorProcessor::ghostBust(AlgoTTMuons& algoTTMuons) {
  AlgoTTMuons selectedTTMuons;
  for(auto& ttMuon : algoTTMuons) {
    for(auto& selected : selectedTTMuons) {
      int ghostBustResult = ghostBustFunc(selected, ttMuon);
      if(ghostBustResult == 0) { //selected kills ttMuon
        ttMuon->kill();
      }
      else if(ghostBustResult == 1) {//ttMuon kills selected
        selected->kill();
      }
      else {
        //ttMuon neither kills nor is killed
      }
    }

    selectedTTMuons.erase(std::remove_if(selectedTTMuons.begin(), selectedTTMuons.end(),
                                  [](AlgoTTMuonPtr& x){return x->isKilled();} ), selectedTTMuons.end());

    if( !ttMuon->isKilled()) {
      selectedTTMuons.push_back(ttMuon);
    }
  }

  return selectedTTMuons;
}


/**should return:
 * 0 if first kills second
 * 1 if second kills first
 * 2 otherwise (none is killed)
 */
int MuCorrelatorProcessor::ghostBust3(std::shared_ptr<AlgoTTMuon> first, std::shared_ptr<AlgoTTMuon> second) {
  //cout<<__FUNCTION__<<":"<<__LINE__<<endl;
  //good ghost bust function looks on the hits indexes in each candidate and check how many hits are common, kill one of them if more then e.g. 1
  int commonHits = 0;
  for(unsigned int iLayer=0; iLayer < first->getStubResults().size(); ++iLayer) {
    if( first->isLayerFired(iLayer) &&
       second->isLayerFired(iLayer) &&
       first->getStubResult(iLayer).getMuonStub() == second->getStubResult(iLayer).getMuonStub() ) { //TODO comparing here just the pointer to the muon stub, in hardware probably it should be an index of the stub
      commonHits++;
    }
  }

  if(commonHits >= 1) { //probably to sharp...
    if(      first->getPdfSum() > second->getPdfSum() ) {
      return 0;
    }
    else if( first->getPdfSum() < second->getPdfSum() ) {
      return 1;
    }
    else {// first->getGpResult().getPdfSum()== second->getGpResult().getPdfSum()
      if( first->getTTTrack()->getPtHw() > second->getTTTrack()->getPtHw() )
        return 0;
      else
        return 1;
    }
  }
  return 2;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

AlgoTTMuons MuCorrelatorProcessor::processTracks(const StandaloneCandWithStubsVec& candsWithStubs, const TrackingTriggerTracks& ttTracks) {
  AlgoTTMuons algoTTMuons;

  for(auto& ttTrack : ttTracks) {
    //TODO add switch to use either processTrack or processTrackUsingRefStubs
    auto algoTTMuon = processTrack(candsWithStubs, ttTrack);

    if(algoTTMuon->getFiredLayerCnt() >= config->nMinFiredLayers())
      algoTTMuons.emplace_back(algoTTMuon);
  }

  auto ghostBustedTTmuons = ghostBust(algoTTMuons);

  return ghostBustedTTmuons;
}

AlgoTTMuonPtr MuCorrelatorProcessor::processTrack(const StandaloneCandWithStubsVec& candsWithStubs, const TrackingTriggerTrackPtr& ttTrack) {
  AlgoTTMuonPtr algoTTMuon = std::make_shared<AlgoTTMuon>(ttTrack, config);

  StandaloneCandWithStubsVec selectedStandaloneCands = selectCandsWithStubs(candsWithStubs, ttTrack);

  for(auto& candWithStubs : selectedStandaloneCands) {
    processTrack(candWithStubs.stubs, ttTrack);
  }

  return algoTTMuon;
}

StandaloneCandWithStubsVec MuCorrelatorProcessor::selectCandsWithStubs(const StandaloneCandWithStubsVec& candsWithStubs, const TrackingTriggerTrackPtr& ttTrack) {
  StandaloneCandWithStubsVec selectedStandaloneCands;

  for(auto& candWithStubs : candsWithStubs) {
//TODO implement
  }

  return selectedStandaloneCands;
}



std::vector<l1t::RegionalMuonCand> MuCorrelatorProcessor::getFinalCandidates(unsigned int iProcessor, l1t::tftype mtfType, AlgoTTMuons& algoTTMuons) {
  std::vector<l1t::RegionalMuonCand> candidates;

  for(auto& algoTTMuon: algoTTMuons) {
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(algoTTMuon->getTTTrack()->getPtHw());
    candidate.setHwEta(algoTTMuon->getTTTrack()->getEtaHw());
    candidate.setHwPhi(config->phiToGlobalHwPhi(algoTTMuon->getTTTrack()->getPhi())); //TODO use hw phi


    candidate.setHwSign(algoTTMuon->getTTTrack()->getCharge() < 0 ? 1 : 0  );
    candidate.setHwSignValid(1);

    unsigned int quality = checkHitPatternValidity(algoTTMuon->getFiredLayerBits() );                      //=4

    candidate.setHwQual (quality);

    std::map<int, int> trackAddr;
    trackAddr[0] = algoTTMuon->getFiredLayerBits().to_ulong();
    trackAddr[1] = 1234;//algoTTMuon->getRefLayer();
    trackAddr[2] = algoTTMuon->getPdfSum();
    trackAddr[3] = algoTTMuon->getTTTrack()->getIndex();
    candidate.setTrackAddress(trackAddr);
    candidate.setTFIdentifiers(iProcessor,mtfType);

    if (candidate.hwPt() > 0 && quality > 0)  //rejecting here the candidates with eta 121, i.e. > 1.31
      candidates.push_back(candidate);
  }
  return candidates;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

bool MuCorrelatorProcessor::checkHitPatternValidity(const boost::dynamic_bitset<>& firedLayerBits) {
  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<boost::dynamic_bitset<> > badPatterns = { //move to the class definition
      //TODO add something if needed
  };

  for(auto badPattern : badPatterns){
    if(firedLayerBits == badPattern) return false;
  }

  return true;
}


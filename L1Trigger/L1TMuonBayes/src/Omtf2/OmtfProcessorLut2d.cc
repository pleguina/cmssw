/*
 * OmtfProcessorLut2d.cc
 *
 *  Created on: May 14, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/Omtf2/OmtfProcessorLut2d.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

std::vector<MuonStubsInput> SimpleStubSelector::selelctStubs(const MuonStubsInput& muonStubs) {
  std::vector<MuonStubsInput> selectedMuonStubs;
  selectedMuonStubs.emplace_back(config);

  //TODO this implementation is rather not possible in the hardware, a different approach would be needed
  for(unsigned int iLayer = 0; iLayer < muonStubs.getMuonStubs().size(); ++iLayer) {
    if(muonStubs.getMuonStubs()[iLayer].size() == 1) {
      selectedMuonStubs.back().addStub(iLayer, muonStubs.getMuonStubs()[iLayer][0]);
    }
  }

  return selectedMuonStubs;
}

OmtfProcessorLut2d::OmtfProcessorLut2d(OMTFConfiguration* omtfConfig):
        config(omtfConfig),
        stubSelector(std::make_unique<SimpleStubSelector>(omtfConfig)),
        ptModule(std::make_shared<PtModuleLut2D>(omtfConfig))
{


}

OmtfProcessorLut2d::OmtfProcessorLut2d(OMTFConfiguration* omtfConfig, const std::shared_ptr<PtModuleLut2D>& ptModule):
        config(omtfConfig),
        stubSelector(std::make_unique<SimpleStubSelector>(omtfConfig)),
        ptModule(ptModule) {

}

OmtfProcessorLut2d::~OmtfProcessorLut2d() {
  // TODO Auto-generated destructor stub
}

AlgoMuon2s OmtfProcessorLut2d::processStubs(const MuonStubsInput& muonStubs) {
  std::vector<MuonStubsInput> selectedMuonStubs = stubSelector->selelctStubs(muonStubs);

  bool areStubs = false;
  for(auto& layer : selectedMuonStubs[0].getMuonStubs() ) {
    if(layer.size() )
      areStubs = true;
  } //TODO crate function in MuonStubsInput

  if(areStubs) {
    LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" selectedMuonStubs\n"<<selectedMuonStubs[0]<<"\n"<<endl;
  }

  AlgoMuon2s algoCands;
  //TODO add handling more then one muon in one proccessor
  AlgoMuon2Ptr algoMuon = ptModule->processStubs(selectedMuonStubs[0]);

  algoCands.emplace_back(algoMuon);

  return algoCands;
}

std::vector<l1t::RegionalMuonCand> OmtfProcessorLut2d::getFinalcandidates(unsigned int iProcessor, l1t::tftype mtfType, const AlgoMuon2s& algoCands) {

  std::vector<l1t::RegionalMuonCand> result;

  for(auto& algoMuon: algoCands) {
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(algoMuon->getPtHw());
    candidate.setHwEta(algoMuon->getEtaHw());

    int phiValue = algoMuon->getPhiHw();
/*    if(phiValue>= int(this->myOmtfConfig->nPhiBins()) ) TODO,
      phiValue -= this->myOmtfConfig->nPhiBins();*/
    ///conversion factor from OMTF to uGMT scale is  5400/576 i.e. phiValue/=9.375;
    phiValue = floor(phiValue*437./pow(2,12));    // ie. use as in hw: 9.3729977
    candidate.setHwPhi(phiValue);

    candidate.setHwSign(algoMuon->getCharge()<0 ? 1:0  );
    candidate.setHwSignValid(1);

    candidate.setHwQual (12);

    std::map<int, int> trackAddr;
    //trackAddr[0] = algoMuon->getFiredLayerBits(); //TODO
    trackAddr[1] = algoMuon->getMuonLikelihood();
    trackAddr[2] = algoMuon->getDisplacedLikelihood();
    trackAddr[3] = algoMuon->getDisplacementHw();

    candidate.setTrackAddress(trackAddr);
    candidate.setTFIdentifiers(iProcessor,mtfType);
    if (candidate.hwPt() > 0)  result.push_back(candidate);
  }
  return result;
}

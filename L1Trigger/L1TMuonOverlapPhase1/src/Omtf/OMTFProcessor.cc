/*
 * OMTFProcessor.cpp
 *
 *  Created on: Oct 7, 2017
 *      Author: kbunkow
 */
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubsInput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GhostBuster.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GhostBusterPreferRefDt.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/StubResult.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

///////////////////////////////////////////////
///////////////////////////////////////////////
template <class GoldenPatternType>
OMTFProcessor<GoldenPatternType>::OMTFProcessor(OMTFConfiguration* omtfConfig,
                                                const edm::ParameterSet& edmCfg,
                                                edm::EventSetup const& evSetup,
                                                const L1TMuonOverlapParams* omtfPatterns)
    : ProcessorBase<GoldenPatternType>(omtfConfig, omtfPatterns) {
  init(edmCfg, evSetup);
};

template <class GoldenPatternType>
OMTFProcessor<GoldenPatternType>::OMTFProcessor(OMTFConfiguration* omtfConfig,
                                                const edm::ParameterSet& edmCfg,
                                                edm::EventSetup const& evSetup,
                                                GoldenPatternVec<GoldenPatternType>&& gps)
    : ProcessorBase<GoldenPatternType>(omtfConfig, std::forward<GoldenPatternVec<GoldenPatternType> >(gps)) {
  init(edmCfg, evSetup);
};

template <class GoldenPatternType>
OMTFProcessor<GoldenPatternType>::~OMTFProcessor() {}

template <class GoldenPatternType>
void OMTFProcessor<GoldenPatternType>::init(const edm::ParameterSet& edmCfg, edm::EventSetup const& evSetup) {
  setSorter(new OMTFSorter<GoldenPatternType>(this->myOmtfConfig->getSorterType()));
  //initialize with the default sorter

  if (this->myOmtfConfig->getGhostBusterType() == "GhostBusterPreferRefDt" ||
      this->myOmtfConfig->getGhostBusterType() == "byLLH" ||
      this->myOmtfConfig->getGhostBusterType() == "byFPLLH"  ||
      this->myOmtfConfig->getGhostBusterType() == "byRefLayer" ) {
    setGhostBuster(new GhostBusterPreferRefDt(this->myOmtfConfig));
    edm::LogVerbatim("OMTFReconstruction") << "setting " << this->myOmtfConfig->getGhostBusterType() << std::endl;
  } else {
    setGhostBuster(new GhostBuster(this->myOmtfConfig));  //initialize with the default sorter
    edm::LogVerbatim("OMTFReconstruction") << "setting GhostBuster" << std::endl;
  }

  edm::LogVerbatim("OMTFReconstruction") << "fwVersion 0x" << hex << this->myOmtfConfig->fwVersion() << std::endl;
}

template <class GoldenPatternType>
std::vector<l1t::RegionalMuonCand> OMTFProcessor<GoldenPatternType>::getFinalcandidates(unsigned int iProcessor,
                                                                                        l1t::tftype mtfType,
                                                                                        const AlgoMuons& algoCands) {
  std::vector<l1t::RegionalMuonCand> result;

  for (auto& myCand : algoCands) {
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(myCand->getPt());
    candidate.setHwEta(myCand->getEtaHw());

    int phiValue = myCand->getPhi();
    if (phiValue >= int(this->myOmtfConfig->nPhiBins()))
      phiValue -= this->myOmtfConfig->nPhiBins();
    phiValue = this->myOmtfConfig->procPhiToGmtPhi(phiValue);
    candidate.setHwPhi(phiValue);

    candidate.setHwSign(myCand->getCharge() < 0 ? 1 : 0);
    candidate.setHwSignValid(1);

    if(myCand->getPtUnconstrained() >= 0) //empty PtUnconstrained is -1, maybe should be corrected on the source
      candidate.setHwPtUnconstrained(myCand->getPtUnconstrained());
    else
      candidate.setHwPtUnconstrained(0);

    unsigned int quality = 12;
    if (this->myOmtfConfig->fwVersion() <= 6)
      quality = checkHitPatternValidity(myCand->getFiredLayerBits()) ? 0 | (1 << 2) | (1 << 3) : 0 | (1 << 2); //12 : 4

    if (abs(myCand->getEtaHw()) == 115 && //115 is eta 1.25
        (static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000001110000000").to_ulong() ||
         static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000001110000000").to_ulong() ||
         static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000000110000000").to_ulong() ||
         static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000001100000000").to_ulong() ||
         static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000001010000000").to_ulong())) {
      if (this->myOmtfConfig->fwVersion() <= 6)
        quality = 4;
      else
        quality = 1;
    }

    if (this->myOmtfConfig->fwVersion() >= 5 && this->myOmtfConfig->fwVersion() <= 6) {
      if (static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000100000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000000000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000100000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000100000000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000000000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000000000110000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000000000110000").to_ulong())
        quality = 1;
    } else if (this->myOmtfConfig->fwVersion() >= 8) {  //TODO fix the fwVersion
      if (static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110000000001").to_ulong() ||

          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011000000000100").to_ulong() ||

          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000011000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000010000000001").to_ulong())
        quality = 1;
      else if (
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010001000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000011000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000011000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000011100000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100001000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100100000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110100000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000111000000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000111000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000111000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000001000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001010000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001010000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001010000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000000111").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100001000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001110000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001110000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010010000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010010000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010010000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010100000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010100000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011110000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011110000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000101000000010101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000010000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000011000000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000011000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000100000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000110000000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001001000000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001001100000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001010000000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000000010000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000000011000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000010000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000100000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000011000000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000011000000000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000010010000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001001000001000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000100000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001110000000111").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000110001000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001110000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000000001000100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000110001000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000000101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001010000001000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001100000001000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("100000010000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000010010000000").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000010100000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000110000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001000000001100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000000000000111101").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000001100000110001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000100000000010100").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000100000000011").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("001000110000000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("010000100010000001").to_ulong() ||
          static_cast<unsigned int>(myCand->getFiredLayerBits()) == std::bitset<18>("000100000000110000").to_ulong())
        quality = 8;
    }  //  if (abs(myCand->getEta()) == 121) quality = 4;
    if (abs(myCand->getEtaHw()) == 121)
      quality = 0;  // changed from 4 on request from HI

    candidate.setHwQual(quality);

    std::map<int, int> trackAddr;
    trackAddr[0] = myCand->getFiredLayerBits();
    trackAddr[1] = myCand->getRefLayer();
    trackAddr[2] = myCand->getDisc();
    trackAddr[3] = myCand->getGpResultUpt().getPdfSumUpt();
    if (candidate.hwPt() > 0) {
      if (ptAssignment) {
        auto pts = ptAssignment->getPts(myCand);
        for (unsigned int i = 0; i < pts.size(); i++) {
          trackAddr[10 + i] = this->myOmtfConfig->ptGevToHw(pts[i]);
        }
      }

      candidate.setTrackAddress(trackAddr);
      candidate.setTFIdentifiers(iProcessor, mtfType);
      result.push_back(candidate);
    }
  }
  return result;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
bool OMTFProcessor<GoldenPatternType>::checkHitPatternValidity(unsigned int hits) {
  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<unsigned int> badPatterns = {
      99840, 34304, 3075, 36928, 12300, 98816, 98944, 33408, 66688, 66176, 7171, 20528, 33856, 35840, 4156, 34880};

  /*
99840 01100001 1000 000000      011000011000000000
34304 00100001 1000 000000      001000011000000000
 3075 00000011 0000 000011      000000110000000011
36928 00100100 0001 000000      001001000001000000
12300 00001100 0000 001100      000011000000001100
98816 01100000 1000 000000      011000001000000000
98944 01100000 1010 000000      011000001010000000
33408 00100000 1010 000000      001000001010000000
66688 01000001 0010 000000      010000010010000000
66176 01000000 1010 000000      010000001010000000
 7171 00000111 0000 000011      000001110000000011
20528 00010100 0000 110000      000101000000110000
33856 00100001 0001 000000      001000010001000000
35840 00100011 0000 000000      001000110000000000
 4156 00000100 0000 111100      000001000000111100
34880 00100010 0001 000000      001000100001000000
   */
  for (auto aHitPattern : badPatterns) {
    if (hits == aHitPattern)
      return false;
  }

  return true;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
AlgoMuons OMTFProcessor<GoldenPatternType>::sortResults(unsigned int iProcessor, l1t::tftype mtfType, int charge) {
  unsigned int procIndx = this->myOmtfConfig->getProcIndx(iProcessor, mtfType);
  return sorter->sortResults(procIndx, this->getPatterns(), charge);
}

template <class GoldenPatternType>
int OMTFProcessor<GoldenPatternType>::extrapolateDtPhiB(const int& refLogicLayer, const int& refPhi, const int& refPhiB, unsigned int targetLayer, const int& targetStubPhi, const int& targetStubQuality,  const int& targetStubR, const OMTFConfiguration* omtfConfig) {

  LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" refLogicLayer "<<refLogicLayer <<" targetLayer "<<targetLayer<<std::endl;

  double hsPhiPitch = 2 * M_PI / omtfConfig->nPhiBins(); //rad/halfStrip

  int phiExtr = 0; //delta phi extrapolated

  float rRefLayer = 431.133; //MB1 i.e. refLogicLayer = 0
  if(refLogicLayer == 2)
    rRefLayer = 512.401; //MB2
  else if(refLogicLayer != 0) {
    return 0;
    //throw cms::Exception("OMTFProcessor<GoldenPatternType>::extrapolateDtPhiB: wrong refStubLogicLayer " + std::to_string(refLogicLayer) );
  }

  if(targetLayer ==  0 || targetLayer ==  2 || targetLayer ==  4 || (targetLayer >= 10 && targetLayer <= 14)) {
    float rTargetLayer = 512.401; //MB2

    if(targetLayer == 0)       rTargetLayer = 431.133; //MB1
    else if(targetLayer == 4)  rTargetLayer = 617.946; //MB3

    else if(targetLayer == 10) rTargetLayer = 413.675; //RB1in
    else if(targetLayer == 11) rTargetLayer = 448.675; //RB1out
    else if(targetLayer == 12) rTargetLayer = 494.975; //RB2in
    else if(targetLayer == 13) rTargetLayer = 529.975; //RB2out
    else if(targetLayer == 14) rTargetLayer = 602.150; //RB3

    if(targetLayer ==  0 || targetLayer ==  2 || targetLayer ==  4) {
      if(targetStubQuality == 2)
        rTargetLayer = rTargetLayer - 23.5/2; //inner superlayer
      else if(targetStubQuality == 3)
        rTargetLayer = rTargetLayer + 23.5/2; //outer superlayer
    }

    float d = rTargetLayer - rRefLayer;
    float deltaPhiExtr = d/rTargetLayer * refPhiB / 512.; //[rad]
    phiExtr = round(deltaPhiExtr / hsPhiPitch); //[halfStrip] //TODO do math as in firmware
    LogTrace("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" deltaPhiExtr "<<deltaPhiExtr<<" phiExtr "<<phiExtr<<std::endl;
  }
  else if(targetLayer ==  1 || targetLayer ==  3 || targetLayer ==  5) {
    int deltaPhi = targetStubPhi - refPhi; //[halfStrip]

    deltaPhi = round(deltaPhi * hsPhiPitch * 512.);
    phiExtr = refPhiB - deltaPhi;
    LogTrace("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" deltaPhi "<<deltaPhi<<" phiExtr "<<phiExtr<<std::endl;
  }
  else if( (targetLayer >= 6 && targetLayer <= 9) || (targetLayer >= 15 && targetLayer <= 17) ) {
    float rME = targetStubR;

    float d = rME - rRefLayer;
    float deltaPhiExtr = d/rME * refPhiB / 512.; //[rad]
    phiExtr = round(deltaPhiExtr / hsPhiPitch); //[halfStrip]
  }
//TODO restrict the range of the phiExtr and refPhiB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return phiExtr;
}

template <class GoldenPatternType>
int OMTFProcessor<GoldenPatternType>::extrapolateDtPhiB(const MuonStubPtr& refStub, const MuonStubPtr& targetStub, unsigned int targetLayer, const OMTFConfiguration* omtfConfig) {
  return OMTFProcessor<GoldenPatternType>::extrapolateDtPhiB(refStub->logicLayer, refStub->phiHw, refStub->phiBHw, targetLayer, targetStub->phiHw, targetStub->qualityHw, targetStub->etaSigmaHw, omtfConfig); //TODO do not use etaSigmaHw!!!!!!
}
///////////////////////////////////////////////
///////////////////////////////////////////////
//const std::vector<OMTFProcessor::resultsMap> &
template <class GoldenPatternType>
void OMTFProcessor<GoldenPatternType>::processInput(unsigned int iProcessor,
                                                    l1t::tftype mtfType,
                                                    const OMTFinput& aInput) {
  unsigned int procIndx = this->myOmtfConfig->getProcIndx(iProcessor, mtfType);
  for (auto& itGP : this->theGPs) {
    for (auto& result : itGP->getResults()[procIndx]) {
      result.reset();
    }
  }

  //////////////////////////////////////
  //////////////////////////////////////
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if (refHitsBits.none())
    return;  // myResults;

  for (unsigned int iLayer = 0; iLayer < this->myOmtfConfig->nLayers(); ++iLayer) {
    //debug
    /*for(auto& h : layerHits) {
      if(h != 5400)
        LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<" "<<__LINE__<<" iLayer "<<iLayer<<" layerHit "<<h<<std::endl;
    }*/
    ///Number of reference hits to be checked.
    unsigned int nTestedRefHits = this->myOmtfConfig->nTestRefHits();

    //loop over all possible refHits, i.e. 128
    for (unsigned int iRefHit = 0; iRefHit < this->myOmtfConfig->nRefHits(); ++iRefHit) {
      if (!refHitsBits[iRefHit])
        continue;
      if (nTestedRefHits-- == 0)
        break;

      const RefHitDef& aRefHitDef = this->myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      unsigned int refLayerLogicNum = this->myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      const MuonStubPtr refStub = aInput.getMuonStub(refLayerLogicNum, aRefHitDef.iInput);
      int phiRef = refStub->phiHw;
      int etaRef = refStub->etaHw;

      unsigned int iRegion = aRefHitDef.iRegion;

      if (this->myOmtfConfig->getBendingLayers().count(iLayer))  //this iLayer is a bending layer
        phiRef = 0;  //then in the delta_phi in process1Layer1RefLayer one obtains simply the iLayer phi

      MuonStubPtrs1D restrictedLayerStubs = this->restrictInput(iProcessor, iRegion, iLayer, aInput);

      //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<" "<<__LINE__<<" iLayer "<<iLayer<<" iRefLayer "<<aRefHitDef.iRefLayer<<" hits.size "<<restrictedLayerHits.size()<<std::endl;
      //LogTrace("l1tOmtfEventPrint")<<"iLayer "<<iLayer<<" refHitNum "<<myOmtfConfig->nTestRefHits()-nTestedRefHits-1<<" iRefHit "<<iRefHit;
      //LogTrace("l1tOmtfEventPrint")<<" nTestedRefHits "<<nTestedRefHits<<" aRefHitDef "<<aRefHitDef<<std::endl;

      std::vector<int> extrapolatedPhi(restrictedLayerStubs.size(), 0);

      //TODO make sure the that the iRefLayer numbers used here corresponds to this in the hwToLogicLayer_0x000X.xml
      if( (this->myOmtfConfig->getUsePhiBExtrapolationMB1() && aRefHitDef.iRefLayer == 0) ||
          (this->myOmtfConfig->getUsePhiBExtrapolationMB2() && aRefHitDef.iRefLayer == 2)    ){
        if((iLayer != refLayerLogicNum) && (iLayer != refLayerLogicNum +1)) {
          LogTrace("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" extrapolating from layer "<<refLayerLogicNum<<" - iRefLayer "<<aRefHitDef.iRefLayer<<std::endl;
          unsigned int iStub = 0;
          for(auto& targetStub : restrictedLayerStubs) {
            if(targetStub)
              extrapolatedPhi[iStub] = extrapolateDtPhiB(refStub, targetStub, iLayer, this->myOmtfConfig);
            iStub++;
          }
        }
      }

      std::vector<int> extrapolatedPhi(restrictedLayerStubs.size(), 0);

      //TODO make sure the that the iRefLayer numbers used here corresponds to this in the hwToLogicLayer_0x000X.xml
      if( (this->myOmtfConfig->getUsePhiBExtrapolationMB1() && aRefHitDef.iRefLayer == 0) ||
          (this->myOmtfConfig->getUsePhiBExtrapolationMB2() && aRefHitDef.iRefLayer == 2)    ){
        if((iLayer != refLayerLogicNum) && (iLayer != refLayerLogicNum +1)) {
          LogTrace("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" extrapolating from layer "<<refLayerLogicNum<<" - iRefLayer "<<aRefHitDef.iRefLayer<<std::endl;
          unsigned int iStub = 0;
          for(auto& targetStub : restrictedLayerStubs) {
            if(targetStub)
              extrapolatedPhi[iStub] = extrapolateDtPhiB(refStub, targetStub, iLayer, this->myOmtfConfig);
            iStub++;
          }
        }
      }

      unsigned int refHitNumber = this->myOmtfConfig->nTestRefHits() - nTestedRefHits - 1;

      int phiExtrp = 0;
      if( (this->myOmtfConfig->getUsePhiBExtrapolationMB1() && aRefHitDef.iRefLayer == 0) ||
          (this->myOmtfConfig->getUsePhiBExtrapolationMB2() && aRefHitDef.iRefLayer == 2)    ){
        phiExtrp = extrapolateDtPhiB(aRefHitDef.iRefLayer, phiRef, refStub->phiBHw, 2, 0, 6, 0, this->myOmtfConfig);
      }

      for (auto& itGP : this->theGPs) {
        if (itGP->key().thePt == 0)  //empty pattern
          continue;

        StubResult stubResult =
            itGP->process1Layer1RefLayer(aRefHitDef.iRefLayer, iLayer, restrictedLayerStubs, extrapolatedPhi, refStub);

        //fixme this unnecessary repeated  for every layer - but in this layout of loops must be like that
        int phiRefSt2 = itGP->propagateRefPhi(phiRef, etaRef, aRefHitDef.iRefLayer);
        //fixme this unnecessary repeated  for every layer

        int phiRefSt2 = itGP->propagateRefPhi(phiRef + phiExtrp, etaRef, aRefHitDef.iRefLayer);
        //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" layerResult: valid"<<layerResult.valid<<" pdfVal "<<layerResult.pdfVal<<std::endl;
        itGP->getResults()[procIndx][refHitNumber].setStubResult(iLayer, stubResult);
        //fixme this unnecessary repeated  for every layer - but in this layout of loops must be like that
        itGP->getResults()[procIndx][refHitNumber].set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef);
        itGP->getResults()[procIndx][refHitNumber].set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef);
        //fixme this unnecessary repeated  for every layer
      }
    }
  }
  //////////////////////////////////////
  //////////////////////////////////////
  {
    for (auto& itGP : this->theGPs) {
      itGP->finalise(procIndx);
      //debug
      /*for(unsigned int iRefHit = 0; iRefHit < itGP->getResults()[procIndx].size(); ++iRefHit) {
        if(itGP->getResults()[procIndx][iRefHit].isValid()) {
          LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<"__LINE__"<<itGP->getResults()[procIndx][iRefHit]<<std::endl;
        }
      }*/
    }
  }

  return;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

template <class GoldenPatternType>
std::vector<l1t::RegionalMuonCand> OMTFProcessor<GoldenPatternType>::run(
    unsigned int iProcessor,
    l1t::tftype mtfType,
    int bx,
    OMTFinputMaker* inputMaker,
    std::vector<std::unique_ptr<IOMTFEmulationObserver> >& observers) {
  //uncomment if you want to check execution time of each method
  //boost::timer::auto_cpu_timer t("%ws wall, %us user in getProcessorCandidates\n");

  std::shared_ptr<OMTFinput> input = std::make_shared<OMTFinput>(this->myOmtfConfig);
  inputMaker->buildInputForProcessor(input->getMuonStubs(), iProcessor, mtfType, bx, bx);

  //LogTrace("l1tOmtfEventPrint")<<"buildInputForProce "; t.report();
  processInput(iProcessor, mtfType, *(input.get()));

  //LogTrace("l1tOmtfEventPrint")<<"processInput       "; t.report();
  AlgoMuons algoCandidates = sortResults(iProcessor, mtfType);

  //LogTrace("l1tOmtfEventPrint")<<"sortResults        "; t.report();
  // perform GB
  AlgoMuons gbCandidates = ghostBust(algoCandidates);

  //LogTrace("l1tOmtfEventPrint")<<"ghostBust"; t.report();
  // fill RegionalMuonCand colleciton
  std::vector<l1t::RegionalMuonCand> candMuons = getFinalcandidates(iProcessor, mtfType, gbCandidates);

  //LogTrace("l1tOmtfEventPrint")<<"getFinalcandidates "; t.report();
  //fill outgoing collection
  for (auto& candMuon : candMuons) {
    candMuon.setHwQual(candMuon.hwQual());
  }

  for (auto& obs : observers) {
    obs->observeProcesorEmulation(iProcessor, mtfType, input, algoCandidates, gbCandidates, candMuons);
  }

  return candMuons;
}

template <class GoldenPatternType>
void OMTFProcessor<GoldenPatternType>::printInfo() const {
  edm::LogVerbatim("OMTFReconstruction") << __PRETTY_FUNCTION__ << std::endl;

  ProcessorBase<GoldenPatternType>::printInfo();
}

/////////////////////////////////////////////////////////

template class OMTFProcessor<GoldenPattern>;
template class OMTFProcessor<GoldenPatternWithStat>;
template class OMTFProcessor<GoldenPatternWithThresh>;

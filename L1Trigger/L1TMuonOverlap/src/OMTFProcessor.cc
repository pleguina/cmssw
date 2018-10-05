/*
 * OMTFProcessor.cpp
 *
 *  Created on: Oct 7, 2017
 *      Author: kbunkow
 */

#include <iostream>
#include <algorithm>
#include <strstream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternParametrised.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"



///////////////////////////////////////////////
///////////////////////////////////////////////
template <class GoldenPatternType>
OMTFProcessor<GoldenPatternType>::OMTFProcessor(const OMTFConfiguration* myOmtfConfig): ProcessorBase<GoldenPatternType>(myOmtfConfig)  {
  setSorter(new OMTFSorter<GoldenPatternType>()); //initialize with the default sorter
  setGhostBuster(new GhostBuster()); //initialize with the default sorter
};

template <class GoldenPatternType>
OMTFProcessor<GoldenPatternType>::~OMTFProcessor() {

}

template <class GoldenPatternType>
std::vector<l1t::RegionalMuonCand> OMTFProcessor<GoldenPatternType>::getFinalcandidates(unsigned int iProcessor, l1t::tftype mtfType, const std::vector<AlgoMuon> & algoCands) {

  std::vector<l1t::RegionalMuonCand> result;

  for(auto myCand: algoCands){
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(myCand.getPt());
    candidate.setHwEta(myCand.getEta());

    int phiValue = myCand.getPhi();
    if(phiValue>= int(this->myOmtfConfig->nPhiBins()) )
      phiValue -= this->myOmtfConfig->nPhiBins();
    ///conversion factor from OMTF to uGMT scale is  5400/576 i.e. phiValue/=9.375;
    phiValue = floor(phiValue*437./pow(2,12));    // ie. use as in hw: 9.3729977
    candidate.setHwPhi(phiValue);

    candidate.setHwSign(myCand.getCharge()<0 ? 1:0  );
    candidate.setHwSignValid(1);

    unsigned int quality = checkHitPatternValidity(myCand.getFiredLayerBits()) ? 0 | (1 << 2) | (1 << 3)
                                                                     : 0 | (1 << 2);
    if (    abs(myCand.getEta()) == 115
        && (    static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001100000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001010000000").to_ulong()
           )
       ) quality =4;
    if( this->myOmtfConfig->fwVersion() >= 5 ) {
      if (    static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000010000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000100000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000001000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000010000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000100000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("001000000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("010000000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000000000011").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000010000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000100000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000001000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000010000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000100000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("001000000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("010000000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000000001100").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000010000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000100000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000001000000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000010000000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000100000000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("001000000000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("010000000000110000").to_ulong()
           || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000000110000").to_ulong()
         ) quality = 1;
}
//  if (abs(myCand.getEta()) == 121) quality = 4;
    if (abs(myCand.getEta()) == 121) quality = 0; // changed on request from HI

    candidate.setHwQual (quality);

    std::map<int, int> trackAddr;
    trackAddr[0] = myCand.getFiredLayerBits();
    trackAddr[1] = myCand.getRefLayer();
    trackAddr[2] = myCand.getDisc();
    candidate.setTrackAddress(trackAddr);
    candidate.setTFIdentifiers(iProcessor,mtfType);
    if (candidate.hwPt() >= 0)  result.push_back(candidate);
  }
  return result;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template<class GoldenPatternType>
bool OMTFProcessor<GoldenPatternType>::checkHitPatternValidity(unsigned int hits) {
  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<unsigned int> badPatterns = {99840, 34304, 3075, 36928, 12300, 98816, 98944, 33408, 66688, 66176, 7171, 20528, 33856, 35840, 4156, 34880};

  /*
99840 01100001 1000 000000
34304 00100001 1000 000000
 3075 00000011 0000 000011
36928 00100100 0001 000000
12300 00001100 0000 001100
98816 01100000 1000 000000
98944 01100000 1010 000000
33408 00100000 1010 000000
66688 01000001 0010 000000
66176 01000000 1010 000000
 7171 00000111 0000 000011
20528 00010100 0000 110000
33856 00100001 0001 000000
35840 00100011 0000 000000
 4156 00000100 0000 111100
34880 00100010 0001 000000
   */
  for(auto aHitPattern: badPatterns){
    if(hits==aHitPattern) return false;
  }

  return true;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template<class GoldenPatternType>
std::vector<AlgoMuon> OMTFProcessor<GoldenPatternType>::sortResults(unsigned int iProcessor, l1t::tftype mtfType, int charge) {
  unsigned int procIndx = this->myOmtfConfig->getProcIndx(iProcessor, mtfType);
  std::vector<AlgoMuon> algoCandidates = sorter->sortResults(procIndx, this->getPatterns(), charge);
  return algoCandidates;
}
///////////////////////////////////////////////
///////////////////////////////////////////////
//const std::vector<OMTFProcessor::resultsMap> &
template<class GoldenPatternType>
const void OMTFProcessor<GoldenPatternType>::processInput(unsigned int iProcessor, l1t::tftype mtfType,
    const OMTFinput & aInput){
  unsigned int procIndx = this->myOmtfConfig->getProcIndx(iProcessor, mtfType);
  for(auto& itGP: this->theGPs) {
    for(auto& result : itGP->getResults()[procIndx]) {
      result.reset();
    }
  }

  //////////////////////////////////////
  //////////////////////////////////////
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if(refHitsBits.none())
    return; // myResults;

  for(unsigned int iLayer=0; iLayer < this->myOmtfConfig->nLayers(); ++iLayer) {
    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);
    /*for(auto& h : layerHits) {
      if(h != 5400)
        std::cout<<__FUNCTION__<<" "<<__LINE__<<" iLayer "<<iLayer<<" layerHit "<<h<<std::endl;
    }*/
    if(layerHits.empty()) continue; //in principle not needed, the size is always 14
    ///Number of reference hits to be checked.
    unsigned int nTestedRefHits = this->myOmtfConfig->nTestRefHits();
    for(unsigned int iRefHit = 0; iRefHit < this->myOmtfConfig->nRefHits(); ++iRefHit) { //loop over all possible refHits, i.e. 128
      if(!refHitsBits[iRefHit]) continue;
      if(nTestedRefHits-- == 0) break;

      const RefHitDef& aRefHitDef = this->myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int phiRef = aInput.getLayerData(this->myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer])[aRefHitDef.iInput];
      int etaRef = aInput.getLayerData(this->myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer],true)[aRefHitDef.iInput];
      unsigned int iRegion = aRefHitDef.iRegion;

      if(this->myOmtfConfig->getBendingLayers().count(iLayer)) //this iLayer is a banding layer
        phiRef = 0;  //then in the delta_phi in process1Layer1RefLayer one obtains simply the iLayer phi

      const OMTFinput::vector1D restrictedLayerHits = this->restrictInput(iProcessor, iRegion, iLayer, layerHits);
      //std::cout<<__FUNCTION__<<" "<<__LINE__<<" iLayer "<<iLayer<<" iRefLayer "<<aRefHitDef.iRefLayer<<" hits.size "<<restrictedLayerHits.size()<<std::endl;
      //std::cout<<"iLayer "<<iLayer<<" refHitNum "<<myOmtfConfig->nTestRefHits()-nTestedRefHits-1<<" iRefHit "<<iRefHit;
      //std::cout<<" nTestedRefHits "<<nTestedRefHits<<" aRefHitDef "<<aRefHitDef<<std::endl;

      int refLayerLogicNumber = this->myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      int refLayerPhiB = 0;
      if(refLayerLogicNumber < 6) {//is DT layer TODO - check
        refLayerPhiB = aInput.getLayerData(refLayerLogicNumber+1)[aRefHitDef.iInput]; //corresponding bending layer has number +1 versus the phi DT layer - this is configured in the XML
        if(this->myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) == 0) {
          throw cms::Exception("not good: the layer is not bending layer");
        }
      }

      for(auto& itGP: this->theGPs) {
        if(itGP->key().thePt == 0) //empty pattern
          continue;
        GoldenPatternResult::LayerResult layerResult = itGP->process1Layer1RefLayer(aRefHitDef.iRefLayer, iLayer,
            phiRef,
            restrictedLayerHits,
            refLayerPhiB);
        int phiRefSt2 = itGP->propagateRefPhi(phiRef, etaRef, aRefHitDef.iRefLayer);
/*        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].setRefPhiRHits(aRefHitDef.iRefLayer, phiRef);
        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].addResult(aRefHitDef.iRefLayer, iLayer,
            aLayerResult.first,
            phiRefSt2, etaRef);*/

        itGP->getResults()[procIndx][this->myOmtfConfig->nTestRefHits()-nTestedRefHits-1].set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef, iLayer, layerResult);
      }
    }
  }
  //////////////////////////////////////
  //////////////////////////////////////
  {
    /*unsigned int iRefHitNum  = 0;
    for(auto & itRefHit: myResults) {
      for(auto & itKey: itRefHit) {
        itKey.second.finalise();
        if(itKey.second.isValid())
          std::cout<<"iRefHitNum "<<iRefHitNum<<" "<<itKey.first<<"\t"<<itKey.second<<std::endl;
      }
      iRefHitNum++;
    }*/

    for(auto& itGP: this->theGPs) {
      itGP->finalise(procIndx);
    }
  }

/*  std::ostringstream myStr;
  myStr<<"iProcessor: "<<iProcessor<<std::endl;
  myStr<<"Input: ------------"<<std::endl;
  myStr<<aInput<<std::endl;
  edm::LogInfo("OMTF processor")<<myStr.str();*/

  return;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

template class OMTFProcessor<GoldenPattern>;
template class OMTFProcessor<GoldenPatternParametrised>;
template class OMTFProcessor<GoldenPatternWithStat>;
template class OMTFProcessor<GoldenPatternWithThresh>;

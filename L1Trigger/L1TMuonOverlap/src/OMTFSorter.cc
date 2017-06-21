#include <cassert>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <bitset>


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorter.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
/*AlgoMuon OMTFSorter::findBestSingleGpResult(const GoldenPatternResult & aResult){
  const GoldenPatternResult::RefLayerResultVec1D& refLayerResults = aResult.getRefLayerResults();

  int bestResultIndex = 0;
  for(unsigned int iRefLayer=0; iRefLayer < refLayerResults.size(); ++iRefLayer) {
    if(refLayerResults[iRefLayer].getFiredLayerCnt() > refLayerResults[bestResultIndex].getFiredLayerCnt()) {
      bestResultIndex = iRefLayer;
    }
    else if (refLayerResults[iRefLayer].getFiredLayerCnt() == refLayerResults[bestResultIndex].getFiredLayerCnt()) {
      if(refLayerResults[iRefLayer].getPdfWeigtSum() > refLayerResults[bestResultIndex].getPdfWeigtSum()) {
        bestResultIndex = iRefLayer;
      }
    }
  }

  if(refLayerResults[bestResultIndex].getFiredLayerCnt() == 0) {
    return AlgoMuon(0, 1024, -10, -1, 0, 0);
  }

  return AlgoMuon(refLayerResults[bestResultIndex], bestResultIndex);
}*/
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
AlgoMuon OMTFSorter::sortRefHitResults(unsigned int iRefHit, const std::vector<IGoldenPattern*>& gPatterns,
					  int charge){

  IGoldenPattern* bestGP = gPatterns[0]; //the GoldenPattern with the best result for this iRefHit
//  std::cout <<" ====== sortRefHitResults: " << std::endl;

  bool foundBest = false;
  for(auto& itGP: gPatterns) {
    if(!itGP->getResults().at(iRefHit).isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(itGP->getResults().at(iRefHit).getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    if(itGP->getResults().at(iRefHit).getFiredLayerCnt() > bestGP->getResults().at(iRefHit).getFiredLayerCnt() ){
      bestGP = itGP;
      foundBest = true;
      //std::cout <<" sorter, byQual, now best is: "<<bestKey << " RefLayer "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<<std::endl;
    }
    else if(itGP->getResults().at(iRefHit).getFiredLayerCnt() == bestGP->getResults().at(iRefHit).getFiredLayerCnt() ) {
      if(itGP->getResults().at(iRefHit).getPdfWeigtSum() > bestGP->getResults().at(iRefHit).getPdfWeigtSum()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGP = itGP;
        foundBest = true;
        //std::cout <<" sorter, byDisc, now best is: "<<bestKey << " "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<< std::endl;
      }
    }
  }

  if(foundBest) {
     AlgoMuon candidate(bestGP->getResults().at(iRefHit), bestGP->key(), iRefHit);
     //std::cout<<__FUNCTION__<<" line "<<__LINE__ <<" return: " << candidate << std::endl;
     return candidate;
  }
  else {
    AlgoMuon candidate;
    candidate.setRefHitNumber(iRefHit);
    return candidate;
  }
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void OMTFSorter::sortRefHitResults(const std::vector<IGoldenPattern*>& gPatterns,
              std::vector<AlgoMuon> & refHitCands,
              int charge){
  
//  for(auto itRefHit: procResults) refHitCands.push_back(sortRefHitResults(itRefHit,charge));
  for (unsigned int iRefHit = 0 ; iRefHit < gPatterns.at(0)->getResults().size(); iRefHit++) {
    AlgoMuon mu = sortRefHitResults(iRefHit, gPatterns, charge);
    mu.setRefHitNumber(iRefHit);
    refHitCands.push_back(mu);
  }
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
bool OMTFSorter::checkHitPatternValidity(unsigned int hits){

  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<unsigned int> badPatterns = {99840, 34304, 3075, 36928, 12300, 98816, 98944, 33408, 66688, 66176, 7171, 20528, 33856, 35840, 4156, 34880};

  for(auto aHitPattern: badPatterns){
    if(hits==aHitPattern) return false;
  }

  return true; 
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
std::vector<l1t::RegionalMuonCand> OMTFSorter::candidates(unsigned int iProcessor, l1t::tftype mtfType, const std::vector<AlgoMuon> & algoCands)
{

  std::vector<l1t::RegionalMuonCand> result;

  for(auto myCand: algoCands){
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(myCand.getPt());
    candidate.setHwEta(myCand.getEta());

    int phiValue = myCand.getPhi();
    if(phiValue>= int(nPhiBins) ) phiValue-=nPhiBins;
    ///conversion factor from OMTF to uGMT scale: 5400/576
//    phiValue/=9.375;
    phiValue *= (437./pow(2,12));    // ie. use as in hw: 9.3729977 
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

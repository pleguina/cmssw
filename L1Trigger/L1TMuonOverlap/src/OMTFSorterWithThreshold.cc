#include <cassert>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <bitset>


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorterWithThreshold.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
AlgoMuon OMTFSorterWithThreshold::sortRefHitResults(unsigned int iRefHit, const std::vector< std::shared_ptr<GoldenPatternWithThresh> >& gPatterns,
					  int charge){

  GoldenPatternWithThresh* bestGP = 0; //the GoldenPattern with the best result for this iRefHit
//  std::cout <<" ====== sortRefHitResults: " << std::endl;

  for(auto& itGP: gPatterns) {
    if(!itGP->getResults().at(iRefHit).isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(itGP->getResults().at(iRefHit).getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    if(itGP->getResults()[iRefHit].getPdfWeigtSum() > itGP->getTreshold(itGP->getResults()[iRefHit].getRefLayer() )) {
      if(bestGP == 0) {
        bestGP = itGP.get();
      }
      else if(itGP->getResults().at(iRefHit).getFiredLayerCnt() >= bestGP->getResults().at(iRefHit).getFiredLayerCnt() ){
        bestGP = itGP.get();
        //std::cout <<" sorter, byQual, now best is: "<<bestKey << " RefLayer "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<<std::endl;
      }//we take the on with the highest pattern number (i.e. pt) among these with the same FiredLayerCnt
    }
  }

  if(bestGP) {
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


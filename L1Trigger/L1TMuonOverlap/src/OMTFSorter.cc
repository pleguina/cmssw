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
AlgoMuon OMTFSorter::sortRefHitResults(unsigned int iRefHit, const std::vector<GoldenPattern*>& gPatterns,
					  int charge){

  GoldenPattern* bestGP = 0; //the GoldenPattern with the best result for this iRefHit
//  std::cout <<" ====== sortRefHitResults: " << std::endl;

  for(auto& itGP: gPatterns) {
    if(!itGP->getResults().at(iRefHit).isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(itGP->getResults().at(iRefHit).getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    if(bestGP == 0) {
      bestGP = itGP;
    }
    else if(itGP->getResults().at(iRefHit).getFiredLayerCnt() > bestGP->getResults().at(iRefHit).getFiredLayerCnt() ){
      bestGP = itGP;
      //std::cout <<" sorter, byQual, now best is: "<<bestKey << " RefLayer "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<<std::endl;
    }
    else if(itGP->getResults().at(iRefHit).getFiredLayerCnt() == bestGP->getResults().at(iRefHit).getFiredLayerCnt() ) {
      if(itGP->getResults().at(iRefHit).getPdfWeigtSum() > bestGP->getResults().at(iRefHit).getPdfWeigtSum()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGP = itGP;
        //std::cout <<" sorter, byDisc, now best is: "<<bestKey << " "<<itKey.second.getRefLayer()<<" FiredLayerCn "<<itKey.second.getFiredLayerCnt()<< std::endl;
      }
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


#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPattern.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"

#include <cassert>
#include <iostream>
#include <algorithm>
#include <bitset>
#include <sstream>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
AlgoMuons::value_type OMTFSorter<GoldenPatternType>::sortRefHitResults(
    unsigned int procIndx, unsigned int iRefHit, const GoldenPatternVec<GoldenPatternType>& gPatterns, int charge) {
  GoldenPatternType* bestGP = nullptr;  //the GoldenPattern with the best result for this iRefHit

  GoldenPatternType* bestGpUnconstr = nullptr;

  bool dbg1078_marker = false;
  std::vector<std::string> dbg1078_buffer;

  for (auto& itGP : gPatterns) {
    if (!itGP->getResults()[procIndx][iRefHit].isValid())
      continue;

    if (charge != 0 && itGP->key().theCharge != charge)
      continue;  //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() < 3)  //TODO - move 3 to the configuration??
      continue;

    if (procIndx == 5) {
      std::ostringstream oss;
      oss << "DBG1078 CAND procIndx=" << procIndx << " iRefHit=" << iRefHit
          << " theNumber=" << itGP->key().number() << " charge=" << itGP->key().theCharge
          << " pt=" << itGP->key().thePt
          << " pdfSum=" << itGP->getResults()[procIndx][iRefHit].getPdfSum()
          << " firedCnt=" << (int)itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt();
      dbg1078_buffer.push_back(oss.str());
    }

    if (itGP->key().thePt == 401 && itGP->key().theCharge == 1 &&
        itGP->getResults()[procIndx][iRefHit].getPdfSum() == 972) {
      dbg1078_marker = true;
    }

    if (bestGP == nullptr) {
      bestGP = itGP.get();
    } else if (myType == 0 && itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() >
                                  bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt()) {
      bestGP = itGP.get();
    } else if (myType == 1 || (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ==
                               bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt())) {
      if (itGP->getResults()[procIndx][iRefHit].getPdfSum() > bestGP->getResults()[procIndx][iRefHit].getPdfSum()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGP = itGP.get();
      }
    }

    if (bestGpUnconstr == nullptr) {
      if (itGP->getResults()[procIndx][iRefHit].getPdfSumUnconstr() > 0)
        bestGpUnconstr = itGP.get();
    } else if (myType == 0 && itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() >
                                  bestGpUnconstr->getResults()[procIndx][iRefHit].getFiredLayerCnt()) {
      bestGpUnconstr = itGP.get();
    } else if (myType == 1 || (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ==
                               bestGpUnconstr->getResults()[procIndx][iRefHit].getFiredLayerCnt())) {
      if (itGP->getResults()[procIndx][iRefHit].getPdfSumUnconstr() >
          bestGpUnconstr->getResults()[procIndx][iRefHit].getPdfSumUnconstr()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGpUnconstr = itGP.get();
      }
    }
  }

  if (dbg1078_marker && procIndx == 5) {
    for (auto& line : dbg1078_buffer) {
      std::cout << line << std::endl;
    }
  }

  if (dbg1078_marker) {
    std::cout << "DBG1078 FINAL procIndx=" << procIndx << " iRefHit=" << iRefHit
              << " winner theNumber=" << (bestGP ? std::to_string(bestGP->key().number()) : std::string("none"))
              << " charge=" << (bestGP ? std::to_string(bestGP->key().theCharge) : std::string("none"))
              << " pt=" << (bestGP ? std::to_string(bestGP->key().thePt) : std::string("none"))
              << " pdfSum="
              << (bestGP ? std::to_string(bestGP->getResults()[procIndx][iRefHit].getPdfSum()) : std::string("none"))
              << std::endl;
  }

  if (bestGP) {
    //this is needed to obtain the same results as in the firmware. for the actual performance it should not matter
    if (bestGP->getResults()[procIndx][iRefHit].getPdfSum() == 0)
      bestGP = gPatterns.at(0).get();

    AlgoMuons::value_type candidate(new AlgoMuon(bestGP->getResults()[procIndx][iRefHit], bestGP, iRefHit));

    if (bestGpUnconstr) {
      candidate->setGpResultUnconstr(bestGpUnconstr->getResults()[procIndx][iRefHit]);
      candidate->setGoldenPaternUnconstr(bestGpUnconstr);
    }

    return candidate;
  } else {
    AlgoMuons::value_type candidate(new AlgoMuon());
    candidate->setRefHitNumber(iRefHit);
    return candidate;
  }
}

template class OMTFSorter<GoldenPattern>;
template class OMTFSorter<GoldenPatternWithStat>;
template class OMTFSorter<GoldenPatternWithThresh>;

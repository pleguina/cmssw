/*
 * ISorter.h
 *
 *  Created on: Jun 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_SORTERBASE_H_
#define OMTF_SORTERBASE_H_

#include <vector>

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"

template <class GoldenPatternType>
class SorterBase {
public:
  virtual ~SorterBase() {}

  virtual std::vector<AlgoMuon> sortResults(const std::vector<GoldenPatternType*>& gPatterns, int charge=0) {
    std::vector<AlgoMuon> refHitCands(gPatterns.at(0)->getResults().size());
  //  for(auto itRefHit: procResults) refHitCands.push_back(sortRefHitResults(itRefHit,charge));
    for (unsigned int iRefHit = 0 ; iRefHit < gPatterns.at(0)->getResults().size(); iRefHit++) {
      AlgoMuon mu = sortRefHitResults(iRefHit, gPatterns, charge);
      mu.setRefHitNumber(iRefHit);
      refHitCands[iRefHit] = mu;
    }
    return refHitCands;
  }


  ///Sort results from a single reference hit.
  ///Select candidate with highest number of hit layers
  ///Then select a candidate with largest likelihood value and given charge
  ///as we allow two candidates with opposite charge from single 10deg region
  virtual AlgoMuon sortRefHitResults(unsigned int iRefHit, const std::vector<GoldenPatternType*>& gPatterns,
        int charge=0)  = 0;
};

#endif /* OMTF_SORTERBASE_H_ */

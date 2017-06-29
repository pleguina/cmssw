/*
 * ISorter.h
 *
 *  Created on: Jun 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_ISORTER_H_
#define OMTF_ISORTER_H_

#include <vector>

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"

class ISorter {
public:
  virtual ~ISorter() {}

  virtual void sortResults(const std::vector<IGoldenPattern*>& gPatterns,
        std::vector<AlgoMuon> & refHitCleanCands,
        int charge=0);

  ///Sort results from a single reference hit.
  ///Select candidate with highest number of hit layers
  ///Then select a candidate with largest likelihood value and given charge
  ///as we allow two candidates with opposite charge from single 10deg region
  virtual AlgoMuon sortRefHitResults(unsigned int iRefHit, const std::vector<IGoldenPattern*>& gPatterns,
        int charge=0)  = 0;
};

#endif /* OMTF_ISORTER_H_ */

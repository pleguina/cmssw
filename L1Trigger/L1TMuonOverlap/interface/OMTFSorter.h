#ifndef OMTF_OMTFSorter_H
#define OMTF_OMTFSorter_H

#include <vector>

#include "L1Trigger/L1TMuonOverlap/interface/ISorter.h"

class OMTFSorter: public ISorter {
public:
  virtual ~OMTFSorter() {}

  ///Sort results from a single reference hit.
  ///Select candidate with highest number of hit layers
  ///Then select a candidate with largest likelihood value and given charge
  ///as we allow two candidates with opposite charge from single 10deg region
  virtual AlgoMuon sortRefHitResults(unsigned int iRefHit, const std::vector<IGoldenPattern*>& gPatterns,
				int charge=0);
};

#endif

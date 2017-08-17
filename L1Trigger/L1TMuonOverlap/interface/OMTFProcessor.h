#ifndef OMTF_OMTFProcessor_H
#define OMTF_OMTFProcessor_H

#include <memory>

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"
#include "L1Trigger/L1TMuonOverlap/interface/SorterBase.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include <L1Trigger/L1TMuonOverlap/interface/SorterBase.h>
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"
#include "L1Trigger/L1TMuonOverlap/interface/IGhostBuster.h"


class OMTFinput;

class SimTrack;

namespace edm{
class ParameterSet;
}

class OMTFProcessor: public ProcessorBase<GoldenPattern> {
 public:

  OMTFProcessor();

  virtual ~OMTFProcessor();
   ///Process input data from a single event
  ///Input data is represented by hits in logic layers expressed in local coordinates
  ///Vector index: number of the ref hit (from 0 to nTestRefHits i.e. 4)
  ///Map key: GoldenPattern key
  //const std::vector<OMTFProcessor::resultsMap> &
  const void processInput(unsigned int iProcessor,
							      const OMTFinput & aInput);
  
  std::vector<AlgoMuon> sortResults(int charge=0);

  std::vector<AlgoMuon> ghostBust(std::vector<AlgoMuon> refHitCands, int charge=0) {
    return ghostBuster->select(refHitCands, charge);
  }

  //convert algo muon to outgoing Candidates
  std::vector<l1t::RegionalMuonCand> getFinalcandidates(
                 unsigned int iProcessor, l1t::tftype mtfType,
                 const std::vector<AlgoMuon> & algoCands);

  ///allows to use other sorter implementation than the default one
  void setSorter(SorterBase<GoldenPattern>* sorter) {
    this->sorter.reset(sorter);
  }

  ///allows to use other IGhostBuster implementation than the default one
  void setGhostBuster(IGhostBuster* ghostBuster) {
    this->ghostBuster.reset(ghostBuster);
  }

 protected:

 private:
  ///Check if the hit pattern of given OMTF candite is not on the list
  ///of invalid hit patterns. Invalid hit patterns provode very little
  ///to efficiency, but gives high contribution to rate.
  ///Candidate with invalid hit patterns is assigned quality=0.
  ///Currently the list of invalid patterns is hardcoded.
  ///This has to be read from configuration.
  bool checkHitPatternValidity(unsigned int hits);

  std::unique_ptr<SorterBase<GoldenPattern> > sorter;

  std::unique_ptr<IGhostBuster> ghostBuster;

};


#endif

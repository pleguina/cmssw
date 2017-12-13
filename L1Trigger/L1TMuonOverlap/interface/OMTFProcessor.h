#ifndef OMTF_OMTFProcessor_H
#define OMTF_OMTFProcessor_H

#include <memory>

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlap/interface/IProcessorEmulator.h"
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
};

template <class GoldenPatternType>
class OMTFProcessor: public ProcessorBase<GoldenPatternType>, public IProcessorEmulator {
 public:

  OMTFProcessor(const OMTFConfiguration* myOmtfConfig);

  virtual ~OMTFProcessor();

  ///Fill GP vec with patterns from CondFormats object
  virtual bool configure(const OMTFConfiguration* omtfParams, const L1TMuonOverlapParams* omtfPatterns) {
    return ProcessorBase<GoldenPatternType>::configure(omtfParams, omtfPatterns);
  }

   ///Process input data from a single event
  ///Input data is represented by hits in logic layers expressed in local coordinates
  ///Vector index: number of the ref hit (from 0 to nTestRefHits i.e. 4)
  ///Map key: GoldenPattern key
  //const std::vector<OMTFProcessor::resultsMap> &
  virtual const void processInput(unsigned int iProcessor, l1t::tftype mtfType,
							      const OMTFinput & aInput);
  
  virtual std::vector<AlgoMuon> sortResults(unsigned int iProcessor, l1t::tftype mtfType, int charge=0);

  virtual std::vector<AlgoMuon> ghostBust(std::vector<AlgoMuon> refHitCands, int charge=0) {
    return ghostBuster->select(refHitCands, charge);
  }

  //convert algo muon to outgoing Candidates
  virtual std::vector<l1t::RegionalMuonCand> getFinalcandidates(
                 unsigned int iProcessor, l1t::tftype mtfType,
                 const std::vector<AlgoMuon> & algoCands);

  ///allows to use other sorter implementation than the default one
  virtual void setSorter(SorterBase<GoldenPatternType>* sorter) {
    this->sorter.reset(sorter);
  }

  ///allows to use other IGhostBuster implementation than the default one
  virtual void setGhostBuster(IGhostBuster* ghostBuster) {
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
  virtual bool checkHitPatternValidity(unsigned int hits);

  std::unique_ptr<SorterBase<GoldenPatternType> > sorter;

  std::unique_ptr<IGhostBuster> ghostBuster;

};

#endif

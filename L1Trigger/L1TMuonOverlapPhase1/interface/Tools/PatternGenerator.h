/*
 * PatternGenerator.h
 *
 *  Created on: Nov 8, 2019
 *      Author: kbunkow
 */

#ifndef OMTF_PATTERNGENERATOR_H_
#define OMTF_PATTERNGENERATOR_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/PatternOptimizerBase.h"

class PatternGenerator: public PatternOptimizerBase {
public:
  PatternGenerator(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps);

  virtual ~PatternGenerator();

  virtual void observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates);

  void endJob();
protected:
  void updateStat();

  void upadatePdfs();

  virtual void saveHists(TFile& outfile);

  //[charge][iLayer]
  std::vector<std::vector<TH2I*> > ptDeltaPhiHists;

  std::vector<unsigned int> eventCntPerGp;
};

#endif /* OMTF_PATTERNGENERATOR_H_ */

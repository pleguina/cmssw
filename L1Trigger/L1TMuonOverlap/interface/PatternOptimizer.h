/*
 * PatternOptimizer.h
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PATTERNOPTIMIZER_H_
#define OMTF_PATTERNOPTIMIZER_H_

#include <functional>

#include "L1Trigger/L1TMuonOverlap/interface/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"


class PatternOptimizer: public IOMTFEmulationObserver {
public:
  PatternOptimizer(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps);
  virtual ~PatternOptimizer();

  virtual void observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
      const std::vector<AlgoMuon>& algoCandidates,
      std::vector<AlgoMuon>& gbCandidates,
      const std::vector<l1t::RegionalMuonCand> & candMuons);

  virtual void observeEventBegin(const edm::Event& iEvent);

  virtual void observeEventEnd(const edm::Event& iEvent);

  virtual void endJob();

  const SimTrack* findSimMuon(const edm::Event &event, const SimTrack* previous = 0);

  static const unsigned int whatExptVal = 0;
  static const unsigned int whatExptNorm = 1;
  static const unsigned int whatOmtfVal = 2;
  static const unsigned int whatOmtfNorm = 3;
private:
  //candidate found by omtf in a given event
  AlgoMuon omtfCand;
  GoldenPatternResult omtfResult;
  GoldenPatternResult exptResult;

  edm::ParameterSet edmCfg;
  //edm::Handle<edm::SimTrackContainer> simTks;

  const OMTFConfiguration* omtfConfig;
  std::vector<std::shared_ptr<GoldenPatternWithStat> > goldenPatterns;

  const SimTrack* simMuon;

  std::function<void (GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp)> updateStatFunc;

  std::function<void (GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate) > updatePdfsFunc;

  void updateStat(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp, double delta, double norm);

  //void updateStatPtDiff_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp);
  void updateStatVoter_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp);

  void updateStatPtDiff2_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp);

  void updateStatPtLogDiff_1(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp);
  void updateStatPtLogDiff_2(GoldenPatternWithStat* omtfCandGp, GoldenPatternWithStat* exptCandGp);

  void updatePdfsMean_1(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate);
  void updatePdfsMean_2(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate);
  void updatePdfsVoter_1(GoldenPatternWithStat* gp, unsigned int& iLayer, unsigned int& iRefLayer, double& learingRate);
};

#endif /* OMTF_PATTERNOPTIMIZER_H_ */

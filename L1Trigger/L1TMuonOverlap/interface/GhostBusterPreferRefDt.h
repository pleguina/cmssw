#ifndef OMTF_GhostBusterPreferRefDt_H
#define OMTF_GhostBusterPreferRefDt_H

#include <vector>
#include <ostream>

#include <map>
#include <set>

#include <memory>

#include "L1Trigger/L1TMuonOverlap/interface/IGhostBuster.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"

class GhostBusterPreferRefDt: public IGhostBuster {
private:
  const OMTFConfiguration* omtfConfig;
public:
  GhostBusterPreferRefDt(const OMTFConfiguration* omtfConfig):omtfConfig(omtfConfig) {};

  ~GhostBusterPreferRefDt() override {};

  AlgoMuons select(AlgoMuons refHitCands, int charge=0) override;

};
#endif

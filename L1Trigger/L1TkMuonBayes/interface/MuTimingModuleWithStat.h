/*
 * MuTimingModuleWithStat.h
 *
 *  Created on: Mar 8, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef L1TkMuonBayes_MUTIMINGMODULEWITHSTAT_H_
#define L1TkMuonBayes_MUTIMINGMODULEWITHSTAT_H_

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TkMuonBayes/interface/MuTimingModule.h"

#include "TH2I.h"

class MuTimingModuleWithStat : public MuTimingModule {
public:
  MuTimingModuleWithStat(const ProcConfigurationBase* config);
  ~MuTimingModuleWithStat() override;

  void process(AlgoMuonBase* algoMuon) override;

  void generateCoefficients();

private:
  //[layer][wheel_ring][etaBin][x=timing, y = 1_Beta]
  std::vector<std::vector<std::vector<TH2I*> > > timigVs1_BetaHists;  //gives average 1/beta

  TH1I* betaDist = nullptr;

  edm::Service<TFileService> fileService;
};

#endif /* L1TkMuonBayes_MUTIMINGMODULEWITHSTAT_H_ */

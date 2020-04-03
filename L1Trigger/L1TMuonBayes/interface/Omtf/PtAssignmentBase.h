/*
 * PtAssignmentBase.h
 *
 *  Created on: Mar 16, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF_PTASSIGNMENTBASE_H_
#define INTERFACE_OMTF_PTASSIGNMENTBASE_H_


#include "L1Trigger/L1TMuonBayes/interface/Omtf/AlgoMuon.h"
#include "lutNN/lutNN2/interface/LutInterNetwork.h"
#include "lutNN/lutNN2/interface/ClassifierToRegression.h"


class PtAssignmentBase {
public:
  PtAssignmentBase(const OMTFConfiguration* omtfConfig): omtfConfig(omtfConfig) {};
  virtual ~PtAssignmentBase();

  virtual std::vector<float> getPts(const AlgoMuons::value_type& algoMuon) = 0;

protected:
  const OMTFConfiguration* omtfConfig = nullptr;
};


class PtAssignmentNN: public PtAssignmentBase {
public:
  PtAssignmentNN(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::string networkFile);
  virtual ~PtAssignmentNN();

  virtual std::vector<float> getPts(const AlgoMuons::value_type& algoMuon);

private:
  lutNN::LutInterNetwork network;

  std::vector<float> ptBins;

  std::vector<std::unique_ptr<lutNN::ClassifierToRegressionBase> > classifierToRegressions;

};

#endif /* INTERFACE_OMTF_PTASSIGNMENTBASE_H_ */

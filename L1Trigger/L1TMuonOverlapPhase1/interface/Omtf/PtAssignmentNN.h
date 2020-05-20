/*
 * PtAssignmentNN.h
 *
 *  Created on: May 8, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF_PTASSIGNMENTNN_H_
#define INTERFACE_OMTF_PTASSIGNMENTNN_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/PtAssignmentBase.h"
#include "lutNN/lutNN2/interface/LutInterNetwork.h"
#include "lutNN/lutNN2/interface/ClassifierToRegression.h"

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


#endif /* INTERFACE_OMTF_PTASSIGNMENTNN_H_ */

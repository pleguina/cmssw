/*
 * PtAssignmentNN.h
 *
 *  Created on: May 8, 2020
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF_PTASSIGNMENTNN_H_
#define INTERFACE_OMTF_PTASSIGNMENTNN_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/PtAssignmentBase.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/LutNetworkFixedPointRegression2Outputs.h"

//#include "lutNN/lutNN2/interface/ClassifierToRegression.h"

class PtAssignmentNNRegression: public PtAssignmentBase {
public:
  PtAssignmentNNRegression(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::string networkFile);
  ~PtAssignmentNNRegression() override;

  std::vector<float> getPts(const AlgoMuons::value_type& algoMuon) override;

private:
  static const int input_I = 10;
  static const int input_F = 4;
  static const std::size_t networkInputSize = 18;

  static const int layer1_neurons = 16;
  static const int layer1_lut_I = 3;
  static const int layer1_lut_F = 13;

  static const int layer1_output_I = 4;
  //4 bits are for the count of the noHit layers which goes to the input of the layer2
  static const int layer2_input_I = 8;

  static const int layer2_neurons = 9;
  static const int layer2_lut_I = 5;
  static const int layer2_lut_F = 11;

  static const int layer3_input_I = 5;

  static const int layer3_0_inputCnt = 8;
  static const int layer3_0_lut_I = 5;
  static const int layer3_0_lut_F = 11;
  static const int output0_I = 8;

  static const int layer3_1_inputCnt = 1;
  static const int layer3_1_lut_I = 4;
  static const int layer3_1_lut_F = 11;
  static const int output1_I = 8;

  lutNN::LutNetworkFixedPointRegression2Outputs<input_I,  input_F,  networkInputSize,
                       layer1_lut_I, layer1_lut_F, layer1_neurons,        //layer1_lutSize = 2 ^ input_I
                       layer1_output_I,
                       layer2_input_I,
                       layer2_lut_I, layer2_lut_F, layer2_neurons,
                       layer3_input_I,
                       layer3_0_inputCnt, layer3_0_lut_I, layer3_0_lut_F, output0_I,
                       layer3_1_inputCnt, layer3_1_lut_I, layer3_1_lut_F, output1_I> lutNetworkFP;

  std::vector<float> ptBins;

  //std::vector<std::unique_ptr<lutNN::ClassifierToRegressionBase> > classifierToRegressions;

};


#endif /* INTERFACE_OMTF_PTASSIGNMENTNN_H_ */

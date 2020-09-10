/*
 * ProcConfigurationBase.cc
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"

ProcConfigurationBase::ProcConfigurationBase(): cscLctCentralBx_(CSCConstants::LCT_CENTRAL_BX) {

}

ProcConfigurationBase::~ProcConfigurationBase() {

}


int ProcConfigurationBase::foldPhi(int phi) const {
  int phiBins = nPhiBins();
  if(phi > phiBins/2)
    return (phi - phiBins );
  else if(phi < -phiBins /2)
    return (phi + phiBins );

  return phi;
}


void ProcConfigurationBase::configureFromEdmParameterSet(const edm::ParameterSet& edmParameterSet) {
  if(edmParameterSet.exists("rpcMaxClusterSize") )
    setRpcMaxClusterSize(edmParameterSet.getParameter<int>("rpcMaxClusterSize"));

  if(edmParameterSet.exists("rpcMaxClusterCnt") )
    setRpcMaxClusterCnt(edmParameterSet.getParameter<int>("rpcMaxClusterCnt"));

  if(edmParameterSet.exists("rpcDropAllClustersIfMoreThanMax") )
    setRpcDropAllClustersIfMoreThanMax(edmParameterSet.getParameter<bool>("rpcDropAllClustersIfMoreThanMax"));

  if(edmParameterSet.exists("lctCentralBx")) {
    cscLctCentralBx_  = edmParameterSet.getParameter<int>("lctCentralBx");
  }

  if(edmParameterSet.exists("minDtPhiQuality")) {
    minDtPhiQuality  = edmParameterSet.getParameter<int>("minDtPhiQuality");
  }

  if(edmParameterSet.exists("minDtPhiBQuality")) {
    minDtPhiBQuality  = edmParameterSet.getParameter<int>("minDtPhiBQuality");
  }
}

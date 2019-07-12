/*
 * ProcConfigurationBase.h
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef MUONBAYES_PROCCONFIGURATIONBASE_H_
#define MUONBAYES_PROCCONFIGURATIONBASE_H_

class ProcConfigurationBase {
public:
  ProcConfigurationBase();
  virtual ~ProcConfigurationBase();

  virtual unsigned int nPhiBins() const = 0;

  virtual double hwPtToGev(unsigned int hwPt) const = 0;

  ///uGMT pt scale conversion: [0GeV, 0.5GeV) = 1 [0.5GeV, 1 Gev) = 2
  virtual int ptGevToHw(double ptGev) const = 0;

  //processor internal phi scale
  virtual int getProcScalePhi(double phiRad, double procPhiZeroRad = 0) const = 0;

  virtual int etaToHwEta(double eta) const = 0;

  //returns address for the  pdf LUTs
  virtual unsigned int ptHwToPtBin(int ptHw) const {
    return 0;
  }

  virtual unsigned int ptGeVToPtBin(float ptGeV) const {
    return 0;
  }

  //returns address for the  pdf LUTs
  virtual unsigned int etaHwToEtaBin(int etaHw) const {
    return 0;
  }

  virtual int foldPhi(int phi) const;

  virtual unsigned int nLayers() const = 0;

  virtual bool isBendingLayer(unsigned int iLayer) const = 0;

  virtual unsigned int getBxToProcess() const {
    return 1;
  }

  bool getRpcDropAllClustersIfMoreThanMax() const {
    return rpcDropAllClustersIfMoreThanMax;
  }

  void setRpcDropAllClustersIfMoreThanMax(bool rpcDropAllClustersIfMoreThanMax =
      true) {
    this->rpcDropAllClustersIfMoreThanMax = rpcDropAllClustersIfMoreThanMax;
  }

  unsigned int getRpcMaxClusterCnt() const {
    return rpcMaxClusterCnt;
  }

  void setRpcMaxClusterCnt(unsigned int rpcMaxClusterCnt = 2) {
    this->rpcMaxClusterCnt = rpcMaxClusterCnt;
  }

  unsigned int getRpcMaxClusterSize() const {
    return rpcMaxClusterSize;
  }

  void setRpcMaxClusterSize(unsigned int rpcMaxClusterSize = 4) {
    this->rpcMaxClusterSize = rpcMaxClusterSize;
  }

protected:
  //parameters of the RpcClusterization
  unsigned int rpcMaxClusterSize = 4;
  unsigned int rpcMaxClusterCnt = 2;

  bool rpcDropAllClustersIfMoreThanMax = true; // if true no  cluster is return if there is more clusters then maxClusterCnt (counted regardless of the size)
};

#endif /* INTERFACE_PROCCONFIGURATIONBASE_H_ */

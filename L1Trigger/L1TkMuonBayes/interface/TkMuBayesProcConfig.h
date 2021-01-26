/*
 * TkMuBayesProcConfig.h
 *
 *  Created on: Jan 30, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef L1TkMuonBayes_TkMuBayesProcConfig_H_
#define L1TkMuonBayes_TkMuBayesProcConfig_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"

#include <cmath>
#include <memory>
#include <vector>

class TkMuBayesProcConfig : public ProcConfigurationBase {
public:
  TkMuBayesProcConfig();

  ~TkMuBayesProcConfig() override {}

  unsigned int nRefLayers() const { return refLayers; };

  double phiGmtUnit = 2. * M_PI / 576;  //TODO from the phase 1 uGMR interface note, correct when new binning defined
  //phi in radians
  virtual int phiToGlobalHwPhi(double phi) const { return std::floor(phi / phiGmtUnit); }

  //phi in radians
  virtual double hwPhiToGlobalPhi(int phi) const { return phi * phiGmtUnit; }

  unsigned int nPhiBins() const override { return phiBins; }

  int getProcScalePhi(double phiRad, double procPhiZeroRad = 0) const override;

  virtual float getProcScalePhiToRad(int phiHw) const;

  double etaUnit =
      0.010875;  //TODO, unit as in the phase 1 upgrade, from the uGMT interface note, should be defined somewhere globally
  ///center of eta bin
  virtual double hwEtaToEta(int hwEta) const { return (hwEta * etaUnit); }

  int etaToHwEta(double eta) const override { return round(eta / etaUnit); }

  double ptUnit = 0.5;  // GeV/unit
  ///uGMT pt scale2 conversion
  double hwPtToGev(int hwPt) const override { return (hwPt - 1.) * ptUnit; }

  ///uGMT pt scale conversion: [0GeV, 0.5GeV) = 1 [0.5GeV, 1 Gev) = 2
  int ptGevToHw(double ptGev) const override { return (ptGev / ptUnit + 1); }

  virtual bool isBarrelLayer(unsigned int layer) const { return !isEndcapLayer(layer); }

  virtual bool isEndcapLayer(unsigned int layer) const;

  virtual bool isEtaLayer(unsigned int layer) const { return !isPhiLayer(layer); }

  virtual bool isPhiLayer(unsigned int layer) const;

  bool isBendingLayer(unsigned int iLayer) const override;

  unsigned int nLayers() const override { return layers; }

  virtual unsigned int nPhiLayers() const { return phiLayers; }

  ///logicLayer - number of the reference layer,
  virtual unsigned int logLayerToRefLayar(unsigned int logicLayer, unsigned int etaBin) const;

  //number of bins in the pdf LUTs
  virtual unsigned int nPtBins() const { return ptBins; }

  //returns address for the  pdf LUTs
  unsigned int ptHwToPtBin(int ptHw) const override;

  unsigned int ptGeVToPtBin(float ptGeV) const override;

  //number of bins in the pdf LUTs - only for endcap layers
  virtual unsigned int nEtaBins() const { return etaBins; }

  //returns address for the  pdf LUTs
  unsigned int etaHwToEtaBin(int etaHw) const override;

  virtual unsigned int nMinFiredLayers() const { return minFiredLayers; }

  virtual unsigned int nMaxMuStubsPerLayer() const { return maxMuStubsPerLayer; }

  virtual float minPdfVal() const { return 0.001; };

  virtual int pdfMaxLogValue() const { return pdfMaxVal; }

  const std::vector<int>& getPtHwBins() const { return ptHwBins; }

  //mode 0 - ptHw, 1 - GeV
  std::string ptBinString(unsigned int ptBin, int mode) const;

  unsigned int getBxToProcess() const override { return bxToProcess; }

  void setBxToProcess(unsigned int bxToProcess = 4) { this->bxToProcess = bxToProcess; }

private:
  unsigned int layers = 30;

  unsigned int refLayers = 1;

  unsigned int phiLayers = 23;

  //phiBins in 2Pi
  unsigned int phiBins = 5400;

  unsigned int ptBins = 64;

  unsigned int etaBins = 16;

  unsigned int bxToProcess =
      4;  //max number of BXes are of the muon stubs that are included in the matching, min should be 1 (the bx=0)

  //minimum number of filers layer to consider the ttRack matching to  stubs
  unsigned int minFiredLayers = 2;

  unsigned int maxMuStubsPerLayer = 200;  //TODO change to the value reasonable for the firmware

  unsigned int pdfMaxVal = 1023;

  //upper edges (ptHw) of the pt bins, the ptHwBins[ptBins] = 'inf'
  std::vector<int> ptHwBins;

  void buildPtHwBins();
};

typedef std::shared_ptr<TkMuBayesProcConfig> TkMuBayesProcConfigPtr;

#endif /* L1TkMuonBayes_TkMuBayesProcConfig_H_ */

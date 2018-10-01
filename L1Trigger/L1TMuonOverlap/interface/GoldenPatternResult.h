#ifndef OMTF_GOLDENPATTERNRESULTS_H
#define OMTF_GOLDENPATTERNRESULTS_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"

//result for one refHit of one GoldenPattern
class GoldenPatternResult {
public:
  typedef std::vector<unsigned int> vector1D;
  struct LayerResult {
    omtfPdfValueType pdfVal = 0;
    bool valid = false;
    int pdfBin = 0; //hit deltaPhi, bin=0 is reserved for no valid hit, see GoldenPatternBase::process1Layer1RefLayer
    int hit = 0;

    LayerResult(omtfPdfValueType pdfVal, bool valid, int pdfBin, int hit): pdfVal(pdfVal), valid(valid), pdfBin(pdfBin), hit(hit) {};
  };

private:
  bool valid;

  //number of the layer from which the reference hit originated
  int refLayer;

  ///Pdf value found for each layer
  ///First index: layer number
  std::vector<omtfPdfValueType> pdfValues;

  ///pdfBins for each for each layer
  ///First index: layer number
  vector1D hitPdfBins;

  vector1D hits;

  ///phi at the 2nd muon station (propagated refHitPhi)
  unsigned int phi;

  ///eta at the 2nd muon station
  unsigned int eta;

  ///Sum of pdfValues
  //omtfPdfValueType
  double pdfSum;

  ///Number of fired layers - excluding banding layers
  unsigned int firedLayerCnt;

  ///bits representing fired logicLayers (including banding layers),
  unsigned int firedLayerBits;

  ///phi of the reference hits
  unsigned int refHitPhi;

  static int finalizeFunction;

  double gpProbability1 = 0;

  double gpProbability2 = 0;
public:
  void init(const OMTFConfiguration* omtfConfig);

  void reset();

  bool isValid() const {
    return valid;
  }

  void setValid(bool valid) {
    this->valid = valid;
  }

  void set(int refLayer, unsigned int phi, unsigned int eta, unsigned int refHitPhi,
      unsigned int iLayer, LayerResult layerResult);

  int getRefLayer() const {
    return this->refLayer;
  }

  void setRefLayer(int refLayer) {
    this->refLayer = refLayer;
  }

  unsigned int getEta() const {
    return eta;
  }

  void setEta(unsigned int eta) {
    this->eta = eta;
  }

  unsigned int getFiredLayerBits() const {
    return firedLayerBits;
  }

  void setFiredLayerBits(unsigned int firedLayerBits) {
    this->firedLayerBits = firedLayerBits;
  }

  unsigned int getFiredLayerCnt() const {
    return firedLayerCnt;
  }

  void setFiredLayerCnt(unsigned int firedLayerCnt) {
    this->firedLayerCnt = firedLayerCnt;
  }

  /*
   * pdfValue from each layer
   */
  const std::vector<omtfPdfValueType>& getPdfValues() const {
    return pdfValues;
  }

  /*
   * sum of the pdfValues in layers
   * if finalise2() it is product of the pdfValues
   */
  omtfPdfValueType getPdfSum() const {
    return pdfSum;
  }

  const vector1D& getHitPdfBins()  {
    return hitPdfBins;
  }

  unsigned int getPhi() const {
    return phi;
  }

  void setPhi(unsigned int phi) {
    this->phi = phi;
  }

  unsigned int getRefHitPhi() const {
    return refHitPhi;
  }

  void setRefHitPhi(unsigned int refHitPhi) {
    this->refHitPhi = refHitPhi;
  }

  bool isLayerFired(unsigned int iLayer) const {
    return firedLayerBits & (1<<iLayer);
  }

  GoldenPatternResult():  valid(0), refLayer(-2), firedLayerCnt(0), myOmtfConfig(0)  {
  };

  //dont use this in the pattern construction, since the myOmtfConfig is null then
  GoldenPatternResult(const OMTFConfiguration* omtfConfig);

  void set();

  void finalise() {
    if(finalizeFunction == 1)
      finalise1();
    else if(finalizeFunction == 2)
      finalise2();
    else
      finalise0();
  }

  //version for the "normal" patterns, i.e. without pdfSum threshold
  void finalise0();

  //version for the patterns with pdfSum threshold
  void finalise1();

  //multiplication of PDF values instead of sum
  void finalise2();

  //bool empty() const;

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternResult & aResult);

  static void setFinalizeFunction(int finalizeFunction_) {
    finalizeFunction = finalizeFunction_;
  }

  double getGpProbability1() const {
    return gpProbability1;
  }

  void setGpProbability1(double probability1 = 0) {
    this->gpProbability1 = probability1;
  }

  double getGpProbability2() const {
    return gpProbability2;
  }

  void setGpProbability2(double probability2 = 0) {
    this->gpProbability2 = probability2;
  }

private:

  const OMTFConfiguration *myOmtfConfig;

};


#endif //OMTF_GOLDENPATTERNRESULTS_H

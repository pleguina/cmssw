#ifndef OMTF_GoldenPattern_H
#define OMTF_GoldenPattern_H

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h>
#include <vector>
#include <ostream>
#undef BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"

class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPattern : public GoldenPatternBase {

public:
  typedef boost::multi_array<short, 3> PdfArrayType;
  typedef boost::multi_array<short, 3> MeanDistPhiArrayType;
  //
  // GoldenPatterns methods
  //
  GoldenPattern(const Key & aKey, unsigned int nLayers, unsigned int nRefLayers, unsigned int nPdfAddrBits): GoldenPatternBase(aKey),
    pdfAllRef(boost::extents[nLayers][nRefLayers][1<<nPdfAddrBits]),
    meanDistPhi(boost::extents[nLayers][nRefLayers][2])
  {
    reset();
  }

  GoldenPattern(const Key & aKey, const OMTFConfiguration* omtfConfig): GoldenPatternBase(aKey, omtfConfig),
      pdfAllRef(boost::extents[omtfConfig->nLayers()][omtfConfig->nRefLayers()][omtfConfig->nPdfBins()]),
      meanDistPhi(boost::extents[omtfConfig->nLayers()][omtfConfig->nRefLayers()][2]) {
    reset();
  }

  virtual ~GoldenPattern() {};

  virtual void setMeanDistPhi(const MeanDistPhiArrayType& aMeanDistPhi) { meanDistPhi = aMeanDistPhi; }

  virtual const MeanDistPhiArrayType& getMeanDistPhi() const { return meanDistPhi; }

  virtual const PdfArrayType& getPdf() const {return pdfAllRef;}

  virtual void setPdf(PdfArrayType& aPdf){  pdfAllRef = aPdf; }

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const;

  virtual int pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const {return pdfAllRef[iLayer][iRefLayer][iBin];}


  virtual void setMeanDistPhiValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int paramIndex = 0) {
    meanDistPhi[iLayer][iRefLayer][paramIndex] = value;
  }

  virtual void setPdfValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) {
    pdfAllRef[iLayer][iRefLayer][iBin] = value;
  }

  friend std::ostream & operator << (std::ostream &out, const GoldenPattern & aPattern);

  ///Add a single count to the relevant pdf bin in three dimensions
  virtual void addCount(unsigned int iRefLayer,
      unsigned int iLayer,
      const int refPhi,
      const OMTFinput::vector1D & layerHits,
      int refLayerPhiB = 0);

  ///Reset contents of all data vectors, keeping the vectors size
  virtual void reset();
/*  virtual void reset(unsigned int nLayers, unsigned int nRefLayers, unsigned int nPdfAddrBits);

  virtual void reset(const OMTFConfiguration* omtfConfig) {
    reset(omtfConfig->nLayers(), omtfConfig->nRefLayers(), omtfConfig->nPdfAddrBits());
  }*/

  ///Normalise event counts in mean dist phi, and pdf vectors to get
  ///the real values of meand dist phi and probability.
  ///The pdf width is passed to this method, since the width stored in
  ///configuration is extended during the pattern making phase.
  virtual void normalise(unsigned int nPdfAddrBits);

  ///Propagate phi from given reference layer to MB2 or ME2
  ///ME2 is used if eta of reference hit is larger than 1.1
  ///expressed in ingerer MicroGMT scale: 1.1/2.61*240 = 101
  virtual int propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer);

  ///Check if the GP has any counts in any of referecne layers;
  virtual bool hasCounts();

protected:
  ///Distributions for all reference layers
  ///First index: measurement layer number
  ///Second index: refLayer number
  ///Third index: pdf bin number within layer 
  PdfArrayType pdfAllRef;

  ///Mean positions in each layer
  ///First index: measurement layer number 
  ///Second index: refLayer number
  ///Third index: index = 0 - a0, index = 1 - a1 for the linear fit meanDistPhi = a0 + a1 * phi_b
  MeanDistPhiArrayType meanDistPhi;

  ///Vector holding number of counts.
  ///Used for making the patterns
  vector2D meanDistPhiCounts;

};

class GoldenPatternWithThresh : public GoldenPattern {
private:
  std::vector<unsigned int> thresholds;
public:
  //
  // GoldenPatterns methods
  //
  GoldenPatternWithThresh(const Key & aKey, unsigned int nLayers, unsigned int nRefLayers, unsigned int nPdfAddrBits):
    GoldenPattern(aKey, nLayers, nRefLayers, nPdfAddrBits) {}

  GoldenPatternWithThresh(const Key & aKey, const OMTFConfiguration* omtfConfig) : GoldenPattern(aKey, omtfConfig) {

  }

  virtual ~GoldenPatternWithThresh() {};

  virtual void reset(const OMTFConfiguration* omtfConfig);

  unsigned int getTreshold(unsigned int iRefLayer) const {
    return thresholds.at(iRefLayer);
  }

  void setThresholds(std::vector<unsigned int>& tresholds) {
    this->thresholds = tresholds;
  }
};
//////////////////////////////////
//////////////////////////////////
#endif 

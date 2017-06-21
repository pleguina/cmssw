#ifndef OMTF_GoldenPattern_H
#define OMTF_GoldenPattern_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/IGoldenPattern.h"

class OMTFConfigMaker;
class OMTFProcessor;
class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPattern : public IGoldenPattern {

public:
  //
  // GoldenPatterns methods
  //
  GoldenPattern(const Key & aKey) : IGoldenPattern(aKey) {}

  GoldenPattern(const Key & aKey, const OMTFConfiguration * omtfConfig) : IGoldenPattern(aKey, omtfConfig) {}

  virtual ~GoldenPattern() {};

  virtual void setMeanDistPhi(const vector2D & aMeanDistPhi){ meanDistPhi = aMeanDistPhi; }

  virtual const vector2D & getMeanDistPhi() const {return meanDistPhi;}

  virtual const vector3D & getPdf() const {return pdfAllRef;}

  virtual void setPdf(const vector3D & aPdf){  pdfAllRef = aPdf; }

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const { return meanDistPhi[iLayer][iRefLayer];}

  virtual int pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const {return pdfAllRef[iLayer][iRefLayer][iBin];}

  ///Process single measurement layer with a single ref layer
  ///Method should be thread safe
  virtual IGoldenPattern::layerResult process1Layer1RefLayer(unsigned int iRefLayer,
      unsigned int iLayer,
      const int refPhi,
      const OMTFinput::vector1D & layerHits,
      int refLayerPhiB = 0);

  friend std::ostream & operator << (std::ostream &out, const GoldenPattern & aPattern);

  ///Add a single count to the relevant pdf bin in three dimensions
  virtual void addCount(unsigned int iRefLayer,
      unsigned int iLayer,
      const int refPhi,
      const OMTFinput::vector1D & layerHits,
      int refLayerPhiB = 0);

  ///Reset contents of all data vectors, keeping the vectors size
  virtual void reset();

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

private:
  ///Distributions for all reference layers
  ///First index: measurement layer number
  ///Second index: refLayer number
  ///Third index: pdf bin number within layer 
  vector3D pdfAllRef;

  ///Mean positions in each layer
  ///First index: measurement layer number 
  ///Second index: refLayer number
  vector2D meanDistPhi;

  ///Vector holding number of counts.
  ///Used for making the patterns
  vector2D meanDistPhiCounts;
};
//////////////////////////////////
//////////////////////////////////
#endif 

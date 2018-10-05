#ifndef OMTF_GoldenPatternPdf4D_H
#define OMTF_GoldenPatternPdf4D_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"

class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPatternPdf4D : public GoldenPatternBase {

public:
  typedef std::vector<vector1D> vector2D;
  typedef std::vector<vector2D> vector3D;
  typedef std::vector<vector3D> vector4D;

  //typedef boost::multi_array<omtfPdfValueType, 4> pdfArrayType;
  typedef vector4D pdfArrayType;

  typedef boost::multi_array<short, 3> meanDistPhiArrayType;

  //
  // GoldenPatternPdf4Ds methods
  //
  GoldenPatternPdf4D(const Key & aKey, const OMTFConfiguration * omtfConfig);

  virtual ~GoldenPatternPdf4D() {};

  //virtual void setMeanDistPhi(const vector2D & aMeanDistPhi){ meanDistPhi = aMeanDistPhi; }

  //virtual const vector2D & getMeanDistPhi() const {return meanDistPhi;}

  virtual const pdfArrayType & getPdf() const {return pdfAllRef;}

  //virtual void setPdf(const vector4D & aPdf){  pdfAllRef = aPdf; }

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const {
    //return meanDistPhi[iLayer][iRefLayer];
    return 0;
  }

  virtual omtfPdfValueType pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const {
    return pdfAllRef[iLayer][iRefLayer][refLayerPhiB][iBin];
  }

  virtual void setMeanDistPhiValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int refLayerPhiB = 0) {
    //meanDistPhi[iLayer][iRefLayer] = value;
  }

  virtual void setPdfValue(omtfPdfValueType value, unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) {
    pdfAllRef[iLayer][iRefLayer][refLayerPhiB][iBin] = value;
  }

  virtual int getDistPhiBitShift(unsigned int iLayer, unsigned int iRefLayer) const  { return 0;};

  virtual void setDistPhiBitShift(int value, unsigned int iLayer, unsigned int iRefLayer) {} ;

  //TODO remove it from this class
  virtual int propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer);

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternPdf4D & aPattern);

  ///Add a single count to the relevant pdf bin in three dimensions
  virtual void addCount(unsigned int iRefLayer,
      unsigned int iLayer,
      const int refPhi,
      const OMTFinput::vector1D & layerHits,
      int refLayerPhiB = 0);

  ///Reset contents of all data vectors, keeping the vectors size
  virtual void reset();

private:
  ///Distributions for all reference layers
  ///First index: measurement layer number
  ///Second index: refLayer number
  ///Third index: pdf bin number of the phiB of the reference layer - if it is DT layer
  ///Fourth index: pdf bin number within layer
  pdfArrayType pdfAllRef;

  ///Mean positions in each layer
  ///First index: measurement layer number 
  ///Second index: refLayer number
  //meanDistPhiArrayType meanDistPhi;

  ///Vector holding number of counts.
  ///Used for making the patterns
  //vector2D meanDistPhiCounts;
};
//////////////////////////////////
//////////////////////////////////
#endif 

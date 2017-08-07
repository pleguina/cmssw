#ifndef OMTF_GoldenPatternPdf4D_H
#define OMTF_GoldenPatternPdf4D_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/IGoldenPattern.h"

class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPatternPdf4D : public IGoldenPattern {

public:
  //
  // GoldenPatternPdf4Ds methods
  //
  GoldenPatternPdf4D(const Key & aKey, const OMTFConfiguration * omtfConfig);

  virtual ~GoldenPatternPdf4D() {};

  virtual Key key() const {return theKey;}

  virtual void setMeanDistPhi(const vector2D & aMeanDistPhi){ meanDistPhi = aMeanDistPhi; }

  //virtual const vector2D & getMeanDistPhi() const {return meanDistPhi;}

  virtual const vector4D & getPdf() const {return pdfAllRef;}

  virtual void setPdf(const vector4D & aPdf){  pdfAllRef = aPdf; }

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const { return meanDistPhi[iLayer][iRefLayer];}

  virtual int pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const {
    return pdfAllRef[iLayer][iRefLayer][refLayerPhiB][iBin];
  }

  virtual void setMeanDistPhiValue(int value, unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) {
    //meanDistPhi[iLayer][iRefLayer] = value;
  }

  virtual void setPdfValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) {
    pdfAllRef[iLayer][iRefLayer][refLayerPhiB][iBin] = value;
  }



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
  vector4D pdfAllRef;

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

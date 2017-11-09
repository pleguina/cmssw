#ifndef OMTF_GoldenPatternParametrised_H
#define OMTF_GoldenPatternParametrised_H

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h>
#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"

class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPatternParametrised : public GoldenPattern {

public:
  GoldenPatternParametrised(const Key & aKey, const OMTFConfiguration* omtfConfig);

  GoldenPatternParametrised(const GoldenPattern* gp);

  virtual ~GoldenPatternParametrised() {};

/*  virtual void setMeanDistPhi(const vector3D& aMeanDistPhi) { meanDistPhi = aMeanDistPhi; }

  virtual const vector3D & getMeanDistPhi() const {return meanDistPhi;}

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const;*/

  //virtual int pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const;

/*  virtual void setMeanDistPhiValue(int value, unsigned int iLayer, unsigned int iRefLayer, int paramIndex = 0) {
    meanDistPhi[iLayer][iRefLayer][paramIndex] = value;
  }*/

/*  virtual void setPdfValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) {
    //pdfAllRef[iLayer][iRefLayer][refLayerPhiB][iBin] = value;
  }*/

  virtual int propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer);

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternParametrised & aPattern);

  ///Reset contents of all data vectors, keeping the vectors size
  virtual void reset(const OMTFConfiguration* omtfConfig);

  typedef double floatType;

  float polynomianValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin) const;
  void generateLuts();


private:
  ///Mean positions in each layer
  ///First index: measurement layer number 
  ///Second index: refLayer number
  //vector2D meanDistPhi;

  ///Mean positions in each layer
  ///First index: measurement layer number
  ///Second index: refLayer number
  ///Third index: index = 0 - a0, index = 1 - a1 for the linear fit meanDistPhi = a0 + a1 * phi_b
  //vector3D meanDistPhi;

  ///Distributions for all reference layers
  ///First index: measurement layer number
  ///Second index: refLayer number
  ///Third index: pdf parameter index
  std::vector<std::vector<std::vector<floatType> > > pdfParamsAllRef;

};
//////////////////////////////////
//////////////////////////////////
#endif 

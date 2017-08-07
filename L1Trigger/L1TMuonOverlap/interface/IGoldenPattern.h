#ifndef OMTF_IGoldenPattern_H
#define OMTF_IGoldenPattern_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h"

class OMTFConfigMaker;
class OMTFProcessor;
class OMTFConfiguration;
//////////////////////////////////
// Key
//////////////////////////////////
struct Key {

Key(int iEta=99, unsigned int iPt=0, int iCharge= 0, unsigned int iNumber=999): 
  theEtaCode(iEta), thePt(iPt), theCharge(iCharge), theNumber(iNumber) {}
    
  inline bool operator< (const Key & o) const {return (theNumber < o.theNumber);}
   
  bool operator==(const Key& o) const {
    //return theNumber==o.theNumber;
    return theEtaCode==o.theEtaCode && thePt==o.thePt && theCharge==o.theCharge && theNumber==o.theNumber;
  }
  
  friend std::ostream & operator << (std::ostream &out, const Key & o) {
    out << "Key_"<<o.theNumber<<": (eta="<<o.theEtaCode<<", pt="<<o.thePt<<", charge="<<o.theCharge<<")";
    return out;
  }

  unsigned int number() const {return theNumber;}

  int theEtaCode;
  unsigned int thePt; //hardware pt, ptInGeV = (thePt-1) * 0.5GeV, where ptInGeV denotes the lover edge of the pt range cover by this pattern
  int          theCharge;
  unsigned int theNumber;

  ///in GeV
  double ptRangeFrom() const;
  ///in GeV
  double ptRangeTo() const;
 };
//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class IGoldenPattern {

 public:

  typedef std::vector<int> vector1D;
  typedef std::vector<vector1D> vector2D;
  typedef std::vector<vector2D> vector3D;
  typedef std::vector<vector3D> vector4D;
  typedef std::pair<int,bool> layerResult;

  //
  // IGoldenPatterns methods
  //
  IGoldenPattern(const Key & aKey) : theKey(aKey), myOmtfConfig(0) {}

  IGoldenPattern(const Key& aKey, const OMTFConfiguration * omtfConfig) : theKey(aKey), myOmtfConfig(omtfConfig),
      results(myOmtfConfig->nTestRefHits(), GoldenPatternResult(omtfConfig)) {}

  virtual ~IGoldenPattern() {}
  
  virtual Key key() const {return theKey;}

  //void setMeanDistPhi(const vector2D & aMeanDistPhi) { meanDistPhi = aMeanDistPhi; }

  //const vector2D & getMeanDistPhi() const {return meanDistPhi;}

  //const vector3D & getPdf() const {return pdfAllRef;}

  //void setPdf(const vector3D & aPdf){  pdfAllRef = aPdf; }

  //virtual const vector2D & getMeanDistPhi() const = 0;

  virtual int meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB = 0) const = 0;

  virtual int pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) const = 0;

  virtual void setMeanDistPhiValue(int value, unsigned int iLayer, unsigned int iRefLayer, int paramIndex = 0) = 0;

  virtual void setPdfValue(int value, unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB = 0) = 0;

  /*
  ///Process single measurement layer with a single ref layer
  ///Method should be thread safe
  virtual IGoldenPattern::layerResult process1Layer1RefLayer(unsigned int iRefLayer,
						    unsigned int iLayer,
						    const int refPhi,
						    const OMTFinput::vector1D & layerHits,
						    int refLayerPhiB = 0) = 0;
    
  ///Add a single count to the relevant pdf bin in three dimensions
  virtual void addCount(unsigned int iRefLayer,
		unsigned int iLayer,
		const int refPhi,
		const OMTFinput::vector1D & layerHits,
		int refLayerPhiB = 0) = 0;

  ///Reset contents of all data vectors, keeping the vectors size
  virtual void reset() = 0;

  ///Normalise event counts in mean dist phi, and pdf vectors to get
  ///the real values of meand dist phi and probability.
  ///The pdf width is passed to this method, since the width stored in
  ///configuration is extended during the pattern making phase.
  virtual void normalise(unsigned int nPdfAddrBits) = 0;

  ///Propagate phi from given reference layer to MB2 or ME2
  ///ME2 is used if eta of reference hit is larger than 1.1
  ///expressed in ingerer MicroGMT scale: 1.1/2.61*240 = 101
  virtual int propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer) = 0;

  ///Check if the GP has any counts in any of referecne layers;
  virtual bool hasCounts() = 0;

  ///last step of the event processing, before sorting and ghost busting
  virtual void finalise() = 0;
  */

  std::vector<GoldenPatternResult>& getResults() {
    return results;
  }

 protected:

  ///Pattern kinematical identification (iEta,iPt,iCharge)
  Key theKey;

  const OMTFConfiguration* myOmtfConfig;

  std::vector<GoldenPatternResult> results;

};
//////////////////////////////////
//////////////////////////////////
#endif 

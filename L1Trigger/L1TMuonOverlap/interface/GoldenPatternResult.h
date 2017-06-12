#ifndef OMTF_GOLDENPATTERNRESULTS_H
#define OMTF_GOLDENPATTERNRESULTS_H

#include <vector>
#include <ostream>

class OMTFConfiguration;

//result for one refHit of one GoldenPattern
class GoldenPatternResult {
public:
  typedef std::vector<unsigned int> vector1D;

private:
  //number of the layer from which the reference hit originated
  int refLayer;

  ///Pdf weight found for each layer
  ///First index: layer number
  vector1D pdfWeights;

  ///phi at the 2nd muon station (propagated refHitPhi)
  unsigned int phi;

  ///eta at the 2nd muon station
  unsigned int eta;

  ///Sum of pdf weights
  unsigned int pdfWeightSum;

  ///Number of fired layers - excluding banding layers
  unsigned int firedLayerCnt;

  ///bits representing fired logicLayers (including banding layers),
  unsigned int firedLayerBits;

  ///phi of the reference hits
  unsigned int refHitPhi;

public:
  void reset();

  bool isValid() const {
    return (refLayer >= 0);
  }

  void set(int refLayer, unsigned int phi, unsigned int eta, unsigned int refHitPhi,
      unsigned int iLayer, unsigned int pdfVal);

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

  const vector1D& getPdfWeights() const {
    return pdfWeights;
  }

  void setPdfWeights(const vector1D& pdfWeigts) {
    this->pdfWeights = pdfWeigts;
  }

  unsigned int getPdfWeigtSum() const {
    return pdfWeightSum;
  }

  void setPdfWeigtSum(unsigned int pdfWeigtSum) {
    this->pdfWeightSum = pdfWeigtSum;
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


  GoldenPatternResult(): refLayer(-2), pdfWeights(8, 0), myOmtfConfig(0)  {
  };

  void configure(const OMTFConfiguration * omtfConfig);

  /*  const GoldenPatternResults::vector2D & getResults() const {return results;}

  const GoldenPatternResults::vector1D & getSummaryVals() const {return results1D;}

  const GoldenPatternResults::vector1D & getSummaryHits() const {return hits1D;}

  const GoldenPatternResults::vector1D & getRefPhis() const {return refPhi1D;}

  const GoldenPatternResults::vector1D & getRefEtas() const {return refEta1D;}

  const GoldenPatternResults::vector1D & getHitsWord() const { return hitsBits;}

  const GoldenPatternResults::vector1D & getRefPhiRHits() const {return refPhiRHit1D;}*/

 /* void setRefPhiRHits(unsigned int iRefLayer, int iRefPhiRHit);

  void addResult(unsigned int iRefLayer,
      unsigned int iLayer,
      unsigned int val,
      int iRefPhi, int iRefEta);
        void clear();
        */

  void set();

  void finalise();



  //bool empty() const;

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternResult & aResult);

private:
  /*  ///Pdf weight found for each layer
  ///First index: layer number
  ///Second index: ref layer number
  vector2D results; 

  ///Reference phi for each reference layer
  vector1D refPhi1D; 

  ///Reference phi for each reference layer
  vector1D refEta1D; 

  ///Sum of pdf weights for each reference layer
  vector1D results1D; 

  ///Number of hits for each reference layer
  vector1D hits1D; 

  ///Words representing nimber of hit layers for each reference layer
  vector1D hitsBits;

  ///Reference phi for each reference layer - the input value
  vector1D refPhiRHit1D; */

  const OMTFConfiguration *myOmtfConfig;

};


#endif //OMTF_GOLDENPATTERNRESULTS_H

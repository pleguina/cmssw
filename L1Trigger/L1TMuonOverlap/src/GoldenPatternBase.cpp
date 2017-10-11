/*
 * GoldenPatternBase.cpp
 *
 *  Created on: Oct 3, 2017
 *      Author: kbunkow
 */


#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h"

////////////////////////////////////////////////////
////////////////////////////////////////////////////
GoldenPatternResult::layerResult GoldenPatternBase::process1Layer1RefLayer(unsigned int iRefLayer,
    unsigned int iLayer,
    const int phiRefHit,
    const OMTFinput::vector1D & layerHits,
    int refLayerPhiB){

  GoldenPatternResult::layerResult aResult; //0, 0

  int phiMean = meanDistPhiValue(iLayer, iRefLayer, refLayerPhiB);
  int phiDist = 1<<(myOmtfConfig->nPdfAddrBits()); //"infinite" value for the beginning
  ///Select hit closest to the mean of probability
  ///distribution in given layer
  for(auto itHit: layerHits){
    if(itHit >= (int)myOmtfConfig->nPhiBins())
      continue;  //empty itHits are marked with nPhiBins() in OMTFProcessor::restrictInput

    if(abs(itHit-phiMean-phiRefHit) < abs(phiDist))
      phiDist = itHit-phiMean-phiRefHit;
  }

  ///Check if phiDist is within pdf range -63 +63
  ///in firmware here the arithmetic "value and sign" is used, therefore the range is -63 +63, and not -64 +63
  if(abs(phiDist) > ( (1<<(myOmtfConfig->nPdfAddrBits()-1)) -1) ) {
    //return aResult;
    return GoldenPatternResult::layerResult(pdfValue(iLayer, iRefLayer, 0), false);
    //in the version with thresholds we use the bin 0 to store the pdf value returned when there was no hit.
    //in the version without thresholds, the value in the bin 0 should be 0
  }

  ///Shift phidist, so 0 is at the middle of the range
  phiDist += 1<<(myOmtfConfig->nPdfAddrBits()-1);

  int pdfVal = pdfValue(iLayer, iRefLayer, phiDist);
  if(pdfVal <= 0)
    return GoldenPatternResult::layerResult(pdfValue(iLayer, iRefLayer, 0), false);
  return GoldenPatternResult::layerResult(pdfVal, true);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternBase::finalise() {
  for(auto& result : getResults()) {
    result.finalise();
  }
}

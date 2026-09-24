/*
 * GoldenPatternBase.cpp
 *
 *  Created on: Oct 3, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternBase.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iomanip>

std::ostream& operator<<(std::ostream& out, const Key& o) {
  out << "Key_" << std::setw(2) << o.theNumber << " hwNum " << std::setw(2) << o.getHwPatternNumber() << " group "
      << std::setw(2) << o.theGroup << ":" << o.theIndexInGroup << " : (eta=" << o.theEtaCode << ", pt=" << std::setw(3)
      << o.thePt << ", charge=" << setw(2) << o.theCharge << ")";
  return out;
}

GoldenPatternBase::GoldenPatternBase(const Key& aKey) : theKey(aKey), myOmtfConfig(nullptr) {}

GoldenPatternBase::GoldenPatternBase(const Key& aKey, const OMTFConfiguration* omtfConfig)
    : theKey(aKey),
      myOmtfConfig(omtfConfig),
      results(boost::extents[myOmtfConfig->processorCnt()][myOmtfConfig->nTestRefHits()]) {
  for (unsigned int iProc = 0; iProc < results.size(); iProc++) {
    for (unsigned int iTestRefHit = 0; iTestRefHit < results[iProc].size(); iTestRefHit++) {
      results[iProc][iTestRefHit].init(omtfConfig);
    }
  }
}

void GoldenPatternBase::setConfig(const OMTFConfiguration* omtfConfig) {
  myOmtfConfig = omtfConfig;
  results.resize(boost::extents[myOmtfConfig->processorCnt()][myOmtfConfig->nTestRefHits()]);
  for (unsigned int iProc = 0; iProc < results.size(); iProc++) {
    for (unsigned int iTestRefHit = 0; iTestRefHit < results[iProc].size(); iTestRefHit++) {
      results[iProc][iTestRefHit].init(omtfConfig);
    }
  }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
StubResult GoldenPatternBase::process1Layer1RefLayer(unsigned int iRefLayer,
                                                     unsigned int iLayer,
                                                     MuonStubPtrs1D layerStubs,
                                                     const std::vector<int>& extrapolatedPhi,
                                                     const MuonStubPtr& refStub) {
  //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<key()<<this->getDistPhiBitShift(iLayer, iRefLayer)<<std::endl;

  int phiMean = this->meanDistPhiValue(iLayer, iRefLayer, refStub->phiBHw);
  int phiDistMin = myOmtfConfig->nPhiBins();
  int deltaPhiSel = myOmtfConfig->nPhiBins();

  ///Select hit closest to the mean of probability
  ///distribution in given layer
  MuonStubPtr selectedStub;

  int phiRefHit = 0;
  if (refStub)
    phiRefHit = refStub->phiHw;

  if (this->myOmtfConfig->isBendingLayer(iLayer)) {
    phiRefHit = 0;  //phi ref hit for the bending layer set to 0, since it should not be included in the phiDist
  }

  for (size_t iStub = 0; iStub < layerStubs.size(); iStub++) {
    const auto& stub = layerStubs[iStub];
    if (!stub)  //empty pointer
      continue;

    int hitPhi = stub->phiHw;
    if (this->myOmtfConfig->isBendingLayer(iLayer)) {
      //rejecting phiB of the low quality DT stubs is done in the OMTFInputMaker
      hitPhi = stub->phiBHw;
    }

    if (hitPhi >= (int)myOmtfConfig->nPhiBins())  //TODO is this needed now? the empty hit will be empty stub
      continue;  //empty itHits are marked with nPhiBins() in OMTFProcessor::restrictInput

    int deltaPhi = hitPhi - extrapolatedPhi[iStub] - phiRefHit;
    int phiDist = deltaPhi - phiMean;
    if ((theKey.theNumber == 3 || theKey.theNumber == 1) ) {
      std::cout << "DBG693 " << __FUNCTION__ << ":" << __LINE__ << " " << theKey
                 << " iRefLayer " << iRefLayer << " iLayer " << iLayer
                 << " iStub " << iStub
                 << " hitPhi " << hitPhi << " extrapolatedPhi " << extrapolatedPhi[iStub]
                 << " phiMean " << phiMean << " phiRefHit " << phiRefHit
                 << " deltaPhi " << deltaPhi << " phiDist " << phiDist
                 << " refStub->phiBHw " << (refStub ? refStub->phiBHw : -99999) << std::endl;
    }

    //firmware works on the sign-value, shift must be done on std::abs(phiDist)
    int sign = phiDist < 0 ? -1 : 1;
    phiDist = std::abs(phiDist) >> this->getDistPhiBitShift(iLayer, iRefLayer);
    phiDist *= sign;
    //if the shift is done here, it means that the phiMean in the xml should be the same as without shift
    //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) std::cout<<__FUNCTION__<<":"<<__LINE__<<" phiDist "<<phiDist<<std::endl;
    if (std::abs(phiDist) < std::abs(phiDistMin)) {
      phiDistMin = phiDist;
      selectedStub = stub;
      deltaPhiSel = deltaPhi;
    }
  }

  bool dbgOn = (theKey.theNumber == 3 || theKey.theNumber == 1);

  if (!selectedStub) {
    PdfValueType pdfVal = 0;
    if (this->myOmtfConfig->isNoHitValueInPdf())
      pdfVal = this->pdfValue(iLayer, iRefLayer, 0);
    if (dbgOn)
      std::cout << "DBG693RES " << theKey << " iRefLayer " << iRefLayer << " iLayer " << iLayer
                 << " NOSTUB pdfVal " << pdfVal << " valid 0" << std::endl;
    return StubResult(pdfVal, false, myOmtfConfig->nPhiBins(), myOmtfConfig->nPhiBins(), iLayer, selectedStub);
  }

  int pdfMiddle = 1 << (myOmtfConfig->nPdfAddrBits() - 1);

  if (iLayer == 8 && dbgOn) {
    std::cout << "DBG693 " << __FUNCTION__ << ":" << __LINE__ << " " << theKey
               << " iRefLayer " << iRefLayer << " iLayer " << iLayer << " SELECTED " << *selectedStub
               << " phiDistMin " << phiDistMin << " phiMean " << phiMean
               << " shift " << this->getDistPhiBitShift(iLayer, iRefLayer) << std::endl;
  }

  ///Check if phiDistMin is within pdf range -63 +63
  ///in firmware here the arithmetic "value and sign" is used, therefore the range is -63 +63, and not -64 +63
  if (std::abs(phiDistMin) > ((1 << (myOmtfConfig->nPdfAddrBits() - 1)) - 1)) {
    if (dbgOn)
      std::cout << "DBG693RES " << theKey << " iRefLayer " << iRefLayer << " iLayer " << iLayer
                 << " OOR phiDistMin " << phiDistMin << " valid 0" << std::endl;
    return StubResult(0, false, phiDistMin + pdfMiddle, deltaPhiSel, iLayer, selectedStub);

    //in some algorithms versions with thresholds we use the bin 0 to store the pdf value returned when there was no hit.
    //in the version without thresholds, the value in the bin 0 should be 0
  }

  ///Shift phidist, so 0 is at the middle of the range
  phiDistMin += pdfMiddle;
  //if (this->getDistPhiBitShift(iLayer, iRefLayer) != 0) LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" phiDistMin "<<phiDistMin<<std::endl;
  PdfValueType pdfVal = this->pdfValue(iLayer, iRefLayer, phiDistMin);
  if (pdfVal <= 0) {
    if (dbgOn)
      std::cout << "DBG693RES " << theKey << " iRefLayer " << iRefLayer << " iLayer " << iLayer
                 << " ZEROPDF phiDistMin " << phiDistMin << " valid 0" << std::endl;
    return StubResult(0, false, phiDistMin, deltaPhiSel, iLayer, selectedStub);
  }
  if (dbgOn)
    std::cout << "DBG693RES " << theKey << " iRefLayer " << iRefLayer << " iLayer " << iLayer
               << " FIRED phiDistMin " << phiDistMin << " pdfVal " << pdfVal << " valid 1" << std::endl;
  return StubResult(pdfVal, true, phiDistMin, deltaPhiSel, iLayer, selectedStub);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternBase::finalise(unsigned int procIndx) {
  for (auto& result : getResults()[procIndx]) {
    result.finalise();
  }
}

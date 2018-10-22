/*
 * PdfGeneratorProcessor.h
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PDFGENERATORPROCESSOR_H_
#define OMTF_PDFGENERATORPROCESSOR_H_

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"

template <class GoldenPatternType>
class PdfGeneratorProcessor: public ProcessorBase<GoldenPatternType>  {
public:
  PdfGeneratorProcessor(OMTFConfiguration* omtfConfig, const L1TMuonOverlapParams* omtfPatterns): ProcessorBase<GoldenPatternType>(omtfConfig, omtfPatterns) {};

  PdfGeneratorProcessor(OMTFConfiguration* omtfConfig, const typename ProcessorBase<GoldenPatternType>::GoldenPatternVec& gps): ProcessorBase<GoldenPatternType>(omtfConfig, gps) {};


  virtual ~PdfGeneratorProcessor() {};


  ///Fill counts for a GoldenPattern of this
  ///processor unit. Pattern key is selcted according
  ///to the SimTrack parameters.
  virtual void fillCounts(unsigned int iProcessor,
      const OMTFinput & aInput,
      const SimTrack* aSimMuon);

  ///Average patterns. Use same meanDistPhi[0] for two
  ///patterns neighboring in pt code.
  //not sutable for patterns with non-zero meanDistPhi[1]
  virtual void averagePatterns();
private:
 ///Shift pdf indexes by differecne between averaged and
 ///original meanDistPhi
/*  virtual void shiftGP(GoldenPatternType *aGP,
        const GoldenPattern::meanDistPhiArrayType& meanDistPhiNew,
        const GoldenPattern::meanDistPhiArrayType& meanDistPhiOld);*/
};

#endif /* OMTF_PDFGENERATORPROCESSOR_H_ */

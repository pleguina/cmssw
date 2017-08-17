/*
 * PdfGeneratorProcessor.h
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PDFGENERATORPROCESSOR_H_
#define OMTF_PDFGENERATORPROCESSOR_H_

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"

class PdfGeneratorProcessor: public ProcessorBase<GoldenPattern>  {
public:
  PdfGeneratorProcessor():ProcessorBase<GoldenPattern>() {};

  virtual ~PdfGeneratorProcessor() {};

  ///Average patterns. Use same meanDistPhi[0] for two
  ///patterns neighboring in pt code.
  //not sutable for patterns with non-zero meanDistPhi[1]
  virtual void averagePatterns();

private:
 ///Shift pdf indexes by differecne between averaged and
 ///original meanDistPhi
  virtual void shiftGP(GoldenPattern *aGP,
        const GoldenPattern::vector3D & meanDistPhiNew,
        const GoldenPattern::vector3D & meanDistPhiOld);
};

#endif /* OMTF_PDFGENERATORPROCESSOR_H_ */

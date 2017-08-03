/*
 * PatternsGeneratorProcessor.h
 *
 *  Created on: Aug 3, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PATTERNSGENERATORPROCESSOR_H_
#define OMTF_PATTERNSGENERATORPROCESSOR_H_

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"

class PatternsGeneratorProcessor: public ProcessorBase<GoldenPattern> {
public:
  PatternsGeneratorProcessor();
  virtual ~PatternsGeneratorProcessor();

  void generatePatterns();
};

#endif /* OMTF_PATTERNSGENERATORPROCESSOR_H_ */

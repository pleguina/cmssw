/*
 * GoldenPatternPdfGen.h
 *
 *  Created on: Oct 23, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_GOLDENPATTERNPDFGEN_H_
#define OMTF_GOLDENPATTERNPDFGEN_H_

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"

class GoldenPatternPdfGen : public GoldenPattern {
public:
  GoldenPatternPdfGen(const Key & aKey, const OMTFConfiguration* omtfConfig);
  virtual ~GoldenPatternPdfGen();

  ///Add a single count to the relevant pdf bin in three dimensions
  virtual void addCount(unsigned int iRefLayer,
      unsigned int iLayer,
      const int refPhi,
      const OMTFinput::vector1D & layerHits,
      int refLayerPhiB = 0);

  ///Normalise event counts in mean dist phi, and pdf vectors to get
  ///the real values of meand dist phi and probability.
  ///The pdf width is passed to this method, since the width stored in
  ///configuration is extended during the pattern making phase.
  virtual void normalise(unsigned int nPdfAddrBits);

  ///Check if the GP has any counts in any of referecne layers;
  virtual bool hasCounts();

protected:
  ///Vector holding number of counts.
  ///Used for making the patterns
  boost::multi_array<int, 2> meanDistPhiCounts;

};

#endif /* OMTF_GOLDENPATTERNPDFGEN_H_ */

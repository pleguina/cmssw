#ifndef OMTF_GoldenPatternWithStat_H
#define OMTF_GoldenPatternWithStat_H

#include <vector>
#include <ostream>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"

class OMTFConfiguration;

//////////////////////////////////
// Golden Pattern
//////////////////////////////////

class GoldenPatternWithStat : public GoldenPattern {
public:
  static const unsigned int STAT_BINS = 4;
  typedef boost::multi_array<float, 4> StatArrayType;
  //GoldenPatternWithStat(const Key & aKey) : GoldenPattern(aKey) {}

  GoldenPatternWithStat(const Key & aKey, unsigned int nLayers, unsigned int nRefLayers, unsigned int nPdfAddrBits):
    GoldenPattern(aKey, nLayers, nRefLayers, nPdfAddrBits),
    statisitics(boost::extents[nLayers][nRefLayers][1<<nPdfAddrBits][STAT_BINS] ) {};

  GoldenPatternWithStat(const Key & aKey, const OMTFConfiguration* omtfConfig): GoldenPattern(aKey, omtfConfig),
    statisitics(boost::extents[omtfConfig->nLayers()][omtfConfig->nRefLayers()][omtfConfig->nPdfBins()][STAT_BINS]) {};

  virtual ~GoldenPatternWithStat() {};

  virtual void updateStat(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, unsigned int what, double value);

  //virtual void updatePdfs(double learingRate);

  friend std::ostream & operator << (std::ostream &out, const GoldenPatternWithStat & aPattern);

  friend class PatternOptimizer;

private:
  StatArrayType statisitics;
};
//////////////////////////////////
//////////////////////////////////
#endif 

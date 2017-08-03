/*
 * ProcessorBase.h
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PROCESSORBASE_H_
#define OMTF_PROCESSORBASE_H_

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/IGoldenPattern.h"

class L1TMuonOverlapParams;
class SimTrack;

template <class GoldenPatternType>
class ProcessorBase {
public:
  ProcessorBase():myOmtfConfig(0)  {
  };

  //virtual ~ProcessorBase();
  virtual ~ProcessorBase() {
    for(auto it: theGPs) delete it;
  }

  ///Just sets the myOmtfConfig
  virtual void setConfigurataion(const OMTFConfiguration* omtfParams) {
    myOmtfConfig = omtfParams;
  }

  ///Return vector of GoldenPatterns
  virtual const std::vector<GoldenPatternType*> & getPatterns() const  {
    return theGPs;
  };

  ///Fill GP vec with patterns from CondFormats object
  virtual bool configure(const OMTFConfiguration * omtfParams, const L1TMuonOverlapParams* omtfPatterns);

  ///Add GoldenPattern to pattern vec.
  ///If GP key already exists in map, a new entry is ignored
  virtual bool addGP(GoldenPatternType *aGP);


  ///Fill counts for a GoldenPattern of this
  ///processor unit. Pattern key is selcted according
  ///to the SimTrack parameters.
  virtual void fillCounts(unsigned int iProcessor,
      const OMTFinput & aInput,
      const SimTrack* aSimMuon);

protected:
  ///vector holding Golden Patterns
  std::vector<GoldenPatternType*> theGPs;

  ///Configuration of the algorithm. This object
  ///does not contain the patterns data.
  const OMTFConfiguration* myOmtfConfig;

  ///Reset all configuration parameters
  virtual void resetConfiguration();

  ///Remove hits whis are outside input range
  ///for given processor and cone
  virtual OMTFinput::vector1D restrictInput(unsigned int iProcessor,
            unsigned int iCone,
            unsigned int iLayer,
            const OMTFinput::vector1D & layerHits);
};

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.hxx"
#endif /* OMTF_PROCESSORBASE_H_ */

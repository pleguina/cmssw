/*
 * TrackMeasureModule.h
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF2_PTMODULELUT2D_H_
#define INTERFACE_OMTF2_PTMODULELUT2D_H_

#include "L1Trigger/L1TMuonBayes/interface/MuonStubsInput.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/AlgoMuon2.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/Lut2d.h"

class IPtModule {
public:
  virtual ~IPtModule() {};

  virtual AlgoMuon2Ptr processStubs(const MuonStubsInput& muonStubs) = 0;

  virtual void saveAsRootHists() = 0;
};

class PtModuleLut2D: public IPtModule {
public:
  PtModuleLut2D(const ProcConfigurationBase* config);

  virtual ~PtModuleLut2D();

  virtual AlgoMuon2Ptr processStubs(const MuonStubsInput& muonStubs);

  virtual void saveAsRootHists();

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
      ar & BOOST_SERIALIZATION_NVP(lut2ds);
  }

  void init();

  const std::vector<std::shared_ptr<Lut2dInter> >& getLut2ds() const {
    return lut2ds;
  }

  //each bit in the interpolate switch on (if 1) or off (if 0) the interpolation in calculation of each output value from the LUT
  void setInterpolate(const boost::dynamic_bitset<>& interpolate) ;


  //if weightDisplacement = true, the weights from LUTs are included in the displacement calculation, otherwise the weight = 1 are used
  void setWeightDisplacement(bool weightDisplacement = true) {
    this->weightDisplacement = weightDisplacement;
  }

  //if weightPt = true, the weights from LUTs are included in the pt calculation, otherwise the weight = 1 are used
  void setWeightPt(bool weightPt = true) {
    this->weightPt = weightPt;
  }

protected:

  const ProcConfigurationBase* config;

  std::vector< std::shared_ptr<Lut2dInter> > lut2ds;

  bool weightPt = true;
  bool weightDisplacement = true;
};

#endif /* INTERFACE_OMTF2_PTMODULELUT2D_H_ */

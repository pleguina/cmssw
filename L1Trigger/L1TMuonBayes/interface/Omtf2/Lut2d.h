/*
 * Lut2d.h
 *
 *  Created on: May 14, 2019
 *      Author: Karol Bunkowski, kbunkow@cern.ch
 */

#include "L1Trigger/L1TMuonBayes/interface/MuonStubsInput.h"

//#include <boost/multi_array.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/dynamic_bitset.hpp>

#ifndef INTERFACE_OMTF2_LUT2D_H_
#define INTERFACE_OMTF2_LUT2D_H_

class Lut2d {
public:
  struct LutConfig {
    static const unsigned int NOT_A_LAYER = 0xffff;

    struct AddressConfig {
      unsigned int layerA = NOT_A_LAYER;
      unsigned int layerB = NOT_A_LAYER;

      unsigned int bitShift = 0;

      unsigned int addrBits = 5;
      unsigned int restBits = 5;

      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
          ar & BOOST_SERIALIZATION_NVP(layerA);
          ar & BOOST_SERIALIZATION_NVP(layerB);

          ar & BOOST_SERIALIZATION_NVP(bitShift);

          ar & BOOST_SERIALIZATION_NVP(addrBits);
          ar & BOOST_SERIALIZATION_NVP(restBits);
      }

      unsigned int getBinsCnt() const {
        return (1<<addrBits);
      }
    };

    AddressConfig addr1;
    AddressConfig addr2;

    unsigned int outValCnt = 2;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(addr1);
        ar & BOOST_SERIALIZATION_NVP(addr2);
        ar & BOOST_SERIALIZATION_NVP(outValCnt);
    }
  };

  Lut2d();

  Lut2d(const LutConfig& lutConfig);

  virtual ~Lut2d();

  const LutConfig& getLutConfig() const {
    return lc;
  }

  virtual std::vector<float> getValues(const MuonStubsInput& muonStubs) const = 0;

  std::string getName();

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
      ar & BOOST_SERIALIZATION_NVP(lc);
  }

protected:

  LutConfig lc;
};

class Lut2dInter: public Lut2d {
public:
 typedef float lutValueType;
 typedef std::vector<std::vector<std::vector<lutValueType> > > lutValuesType;

  Lut2dInter();

  Lut2dInter(const LutConfig& lutConfig);

  virtual ~Lut2dInter();

  virtual std::vector<lutValueType> getValues(unsigned int addr1, unsigned int rest1, unsigned int addr2, unsigned int rest2) const;

  virtual std::vector<lutValueType> getValues(const MuonStubsInput& muonStubs) const;

  //adds incremet to the value at valIndx corresponding to a given stub
  virtual bool updateValue(const MuonStubsInput& muonStubs, unsigned int valIndx, float incremet);

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
      //ar & boost::serialization::base_object<Lut2d>(*this);

      ar & boost::serialization::make_nvp("Lut2d", boost::serialization::base_object<Lut2d>(*this));

      ar & BOOST_SERIALIZATION_NVP(lutValues);
  }

  friend std::ostream & operator<< (std::ostream &out, const Lut2dInter& lut);

  const lutValuesType& getLutValues() const {
    return lutValues;
  }

  lutValuesType& getLutValues() {
    return lutValues;
  }

  const boost::dynamic_bitset<>& getInterpolate() const {
    return interpolate;
  }

  void setInterpolate(const boost::dynamic_bitset<>& interpolate) {
    this->interpolate = interpolate;
  }

  friend class PtModuleLut2D;

protected:
  bool getInput(const MuonStubsInput& muonStubs, const LutConfig::AddressConfig& addrConf, unsigned int& addr, unsigned int& rest) const;

  //boost::multi_array<int, 3> lutValues;  would be better, but cannon be serialized with boost/serialization
  lutValuesType lutValues;

  boost::dynamic_bitset<> interpolate;
};

#endif /* INTERFACE_OMTF2_LUT2D_H_ */

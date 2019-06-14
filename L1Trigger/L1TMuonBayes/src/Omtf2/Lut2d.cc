/*
 * Lut2d.cpp
 *
 *  Created on: May 14, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonBayes/interface/Omtf2/Lut2d.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iomanip>

Lut2d::Lut2d() {

}

Lut2d::Lut2d(const LutConfig& lutConfig): lc(lutConfig) {
}

Lut2d::~Lut2d() {
  // TODO Auto-generated destructor stub
}

std::string Lut2d::getName() {
  std::ostringstream ostr;

  ostr<<"LUT2D";
  ostr<<"_layer1A_"<<getLutConfig().addr1.layerA;
  if(getLutConfig().addr1.layerB != LutConfig::NOT_A_LAYER) {
    ostr<<"_layer1B_"<<getLutConfig().addr1.layerB;
  }

  ostr<<"_layer2A_"<<getLutConfig().addr2.layerA;
  if(getLutConfig().addr2.layerB != LutConfig::NOT_A_LAYER) {
    ostr<<"_layer2B_"<<getLutConfig().addr2.layerB;
  }

  return ostr.str();
}

Lut2dInter::Lut2dInter(): Lut2d() {

}

Lut2dInter::Lut2dInter(const LutConfig& lutConfig): Lut2d(lutConfig),
    //lutValues(boost::extents[lc.outValCnt][1<<lc.addr1.addrBits][1<<lc.addr2.addrBits])
    lutValues(lc.outValCnt, std::vector<std::vector<Lut2dInter::lutValueType> >( 1<<lc.addr1.addrBits, std::vector<Lut2dInter::lutValueType>(1<<lc.addr2.addrBits, 0) ) ),
    interpolate(lc.outValCnt, 0xfffff)
{
  for(unsigned int iOut = 0; iOut < lutValues.size() -1; iOut++) { //the last value is for string the number of events - needed for debug etc
    for(unsigned int addr1 = 0; addr1 < lutValues[iOut].size(); addr1++) {
      for(unsigned int addr2 = 0; addr2 < lutValues[iOut][addr1].size(); addr2++) {
        lutValues[iOut][addr1][addr2] = 0;//(iOut<< 10 ) | (addr1 << 5) | addr2;
      }
    }
  }
}

Lut2dInter::~Lut2dInter() {

}

std::ostream & operator<< (std::ostream &out, const Lut2dInter& lut) {
  out <<"Lut2dInter:"<<
      " | addr1.layerA "<<std::setw(2)<<(lut.lc.addr1.layerA)<<
      " | addr1.layerB "<<std::setw(2)<<(lut.lc.addr1.layerB)<<
      " | addr2.layerA "<<std::setw(2)<<(lut.lc.addr2.layerA)<<
      " | addr2.layerB "<<std::setw(2)<<(lut.lc.addr2.layerB)<<" | ";
  return out;
}

std::vector<float> Lut2dInter::getValues(unsigned int addr1, unsigned int rest1, unsigned int addr2, unsigned int rest2) const {
  std::vector<float> outVals(lutValues.size(), 0);

  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" | addr1 "<<addr1<<" | rest1 "<<rest1<<" || addr2 "<<addr2<<" | rest2 "<<rest2<<" | "<<std::endl;

  for(unsigned int iOut = 0; iOut < lutValues.size(); iOut++) {
    if(addr1 >= lutValues[iOut].size() -1)
      return outVals;

    if(addr2 >= lutValues[iOut][addr1].size() -1)
      return outVals;

    lutValueType val = 0;

    unsigned int restMask = (1 << lc.addr1.restBits) -1; //restMask = 31 if restBits = 5
    if(interpolate[iOut]) {
      if( (rest1 + rest2) < restMask) {//fixme assuming here rest1Mask = rest2Mask
        val = lutValues[iOut][addr1][addr2] +
            (((lutValues[iOut][addr1+1][addr2] - lutValues[iOut][addr1][addr2]) * rest1 ) / (restMask + 1) ) + //division by (restMask + 1) i.e bin size can be done with  >>addr1BitShift, but then one have to be careful with the negative LUT values
            (((lutValues[iOut][addr1][addr2+1] - lutValues[iOut][addr1][addr2]) * rest2 ) / (restMask + 1))   ;
      }
      else {
        val = lutValues[iOut][addr1+1][addr2+1] -
            (((lutValues[iOut][addr1+1][addr2+1] - lutValues[iOut][addr1][addr2+1]) * (restMask +1 - rest1) ) / (restMask + 1)) -
            (((lutValues[iOut][addr1+1][addr2+1] - lutValues[iOut][addr1+1][addr2]) * (restMask +1 - rest2) ) / (restMask + 1))   ;
      }
    }
    else {
      int half1 = (rest1 >> (lc.addr1.restBits-1) ); //increasing addr by 1 if rest1 is greater then half
      int half2 = (rest2 >> (lc.addr2.restBits-1) );
      val = lutValues.at(iOut)[addr1 + half1][addr2 + half2];
    }
    outVals[iOut] = val;
  }
  return outVals;
}

bool Lut2dInter::getInput(const MuonStubsInput& muonStubs, const LutConfig::AddressConfig& addrConf, unsigned int& addr, unsigned int& rest) const {
  int inVal = muonStubs.getPhiHw(addrConf.layerA, 0);
  if(inVal == MuonStub::EMTPY_PHI)
    return false;

  if(addrConf.layerB != LutConfig::NOT_A_LAYER) {
    int phi2 = muonStubs.getPhiHw(addrConf.layerB, 0);
    if(phi2 == MuonStub::EMTPY_PHI)
        return false;
    inVal -= phi2; //calculation of delta phi
  }

  inVal = inVal << addrConf.bitShift; //if the max value of inVal is less then then 3/4 of the mask, then it is worth to increase the bitShift, to better utilize the LUT

  unsigned int inValMaxBits = addrConf.addrBits + addrConf.restBits;
  unsigned int mask = (1 << inValMaxBits ) -1;

  unsigned int inValS = inVal + (mask >> 1); //shifting the inVal by half of the possible value, i.e. mask/2 to assure that it is usually positive

  if( (inValS & ~mask) != 0)
    return false;

  unsigned int restMask = (1 << addrConf.restBits) -1;
  rest = inValS & restMask;
  addr = inValS >> addrConf.restBits;
  return true;
}

//at this poin the muonStubs must contain at most one stub per layer
std::vector<float> Lut2dInter::getValues(const MuonStubsInput& muonStubs) const {
  std::vector<float> outVals(lutValues.size(), 0);

  unsigned int addr1 = 0;
  unsigned int rest1 = 0;
  if(!getInput(muonStubs, lc.addr1, addr1, rest1) )
    return outVals;

  unsigned int addr2 = 0;
  unsigned int rest2 = 0;
  if(!getInput(muonStubs, lc.addr2, addr2, rest2) )
    return outVals;

  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<std::endl;
  return getValues(addr1, rest1, addr2, rest2);
}

bool Lut2dInter::updateValue(const MuonStubsInput& muonStubs, unsigned int valIndx, float incremet) {
  unsigned int addr1 = 0;
  unsigned int rest1 = 0;
  if(!getInput(muonStubs, lc.addr1, addr1, rest1) )
    return false;

  unsigned int addr2 = 0;
  unsigned int rest2 = 0;
  if(!getInput(muonStubs, lc.addr2, addr2, rest2) )
    return false;

  addr1 += (rest1 >> (lc.addr1.restBits-1) ); //increasing addr by 1 if rest1 is greater then half
  addr2 += (rest2 >> (lc.addr2.restBits-1) );

  if(addr1 >= lutValues.at(valIndx).size())
    return false;

  if(addr2 >= lutValues[valIndx][addr1].size())
    return false;

  lutValues[valIndx][addr1][addr2] += incremet;
  return true;
}

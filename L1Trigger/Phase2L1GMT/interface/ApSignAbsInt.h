/*
 * ApSingAbsInt.h
 *
 *  Created on: Nov 24, 2021
 *      Author: kbunkow
 */

#ifndef INTERFACE_APSIGNABSINT_H_
#define INTERFACE_APSIGNABSINT_H_

#include <ap_int.h>

/*
 * sing plus absolute value integer representation for the hls ap_int
 */
template <int _AP_W>
struct ApSignAbsInt {
public:
  ApSignAbsInt():
    sign(0), absValue(0) {};

  ApSignAbsInt(const ap_uint<1>& sign, const ap_uint<_AP_W>& absValue):
    sign(sign), absValue(absValue) {};

  ApSignAbsInt(uint sign, uint absValue):
    sign(sign), absValue(absValue) {};

  virtual ~ApSignAbsInt() {};

  ap_uint<_AP_W+1> toRaw() {
    return (ap_uint<_AP_W+1>(sign) << _AP_W) | absValue;
  }

  ap_int<_AP_W+1> to_signed() {
    return (sign.to_uchar() == 1 ? (ap_int<1>(-1) * absValue) : ap_int<_AP_W+1>(absValue));
  }

  auto to_int() {
    return to_signed().to_int();
  }

  auto to_string() {
    return to_signed().to_string();
  }

  ap_uint<1> sign; //0 is positive, 1 is negative
  ap_uint<_AP_W> absValue;
};

#endif /* INTERFACE_APSIGNABSINT_H_ */

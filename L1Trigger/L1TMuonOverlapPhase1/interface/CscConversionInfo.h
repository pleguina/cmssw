#ifndef L1T_OmtfP1_CSCCONVERSIONINFO_H
#define L1T_OmtfP1_CSCCONVERSIONINFO_H

struct CscConversionInfo {
  int phi = 0;
  int offset = 0;  // fixOff
  double scale = 0;
  int order = 0;
  int halfStrip = 0;  // halfStrip value used in phi calculation
};

#endif

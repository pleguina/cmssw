#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfPhase2AngleConverter.h"

namespace {
  int sgn(float val) { return (0 < val) - (val < 0); }

  int etaVal2CodePhase2(float etaVal) {
    int sign = sgn(etaVal);
    int code = (int)round(fabs(etaVal) * 115 / 1.25);
    return sign * code;
  }
}  // namespace

int OmtfPhase2AngleConverter::getProcessorPhi(int phiZero, l1t::tftype part, int dtScNum, int dtPhi) const {
  constexpr int dtPhiBins = 65536;          //65536. for [-0.5,0.5] radians
  double hsPhiPitch = 2 * M_PI / nPhiBins;  // width of phi Pitch, related to halfStrip at CSC station 2

  int sector = dtScNum + 1;  //NOTE: there is a inconsistency in DT sector numb. Thus +1 needed to get detector numb.

  double scale = 0.5 / dtPhiBins / hsPhiPitch;  //was 0.8
  int scale_coeff = lround(scale * (1 << 15));

  int ichamber = sector - 1;
  if (ichamber > 6)
    ichamber = ichamber - 12;

  int offsetGlobal = (int)nPhiBins * ichamber / 12;

  int phiConverted = ((dtPhi * scale_coeff) >> 15) + offsetGlobal - phiZero;

  return config->foldPhi(phiConverted);
}

int OmtfPhase2AngleConverter::getGlobalEta(DTChamberId dTChamberId,
                                           const L1Phase2MuDTThContainer* dtThDigis,
                                           int bxNum) const {
  // In firmware export mode use fixed mid-chamber eta values that match the RTL constants.
  // In sample production mode fall through to the theta-digi LUT lookup below.
  if (dtFixedPointEtaForFirmware_) {
    if (dTChamberId.station() == 1)
      return 92;
    else if (dTChamberId.station() == 2)
      return 79;
    else if (dTChamberId.station() == 3)
      return 75;
    return 95;
  }

  // Sample production mode: real theta-digi LUT lookup
  int dtThBins = 65536;  //65536. for [-6.3,6.3]
  float kconv = 1 / (dtThBins / 2.);

  float eta = -999;
  bool foundeta = false;
  int thetaDigiCnt = 0;
  for (const auto& thetaDigi : (*(dtThDigis->getContainer()))) {
    if (thetaDigi.whNum() == dTChamberId.wheel() && thetaDigi.stNum() == dTChamberId.station() &&
        thetaDigi.scNum() == (dTChamberId.sector() - 1) && (thetaDigi.bxNum() - 20) == bxNum) {
      float k = thetaDigi.k() * kconv;
      int sign = sgn(thetaDigi.z());
      eta = -1. * sign * log(fabs(tan(atan(1 / k) / 2.)));
      LogTrace("OMTFReconstruction") << "OmtfPhase2AngleConverter::getGlobalEta(" << dTChamberId << ") eta: " << eta
                                     << " k: " << k << " thetaDigi.k(): " << thetaDigi.k();
      thetaDigiCnt++;
      if ((dTChamberId.station() == 1 && (std::abs(eta) < 0.85 || std::abs(eta) > 1.20)) ||
          (dTChamberId.station() == 2 && (std::abs(eta) < 0.75 || std::abs(eta) > 1.04)) ||
          (dTChamberId.station() == 3 && (std::abs(eta) < 0.63 || std::abs(eta) > 0.92))) {
        foundeta = false;
      } else
        foundeta = true;
    }
  }

  // If more than 1 thetaDigi per chamber they are ambiguous - fall back to mid-chamber
  if (thetaDigiCnt > 1)
    foundeta = false;

  if (foundeta) {
    return std::abs(etaVal2CodePhase2(eta));
  } else {
    // Fall back to mid-chamber value
    if (dTChamberId.station() == 1)
      return 92;
    else if (dTChamberId.station() == 2)
      return 79;
    else if (dTChamberId.station() == 3)
      return 75;
    return 95;
  }
}

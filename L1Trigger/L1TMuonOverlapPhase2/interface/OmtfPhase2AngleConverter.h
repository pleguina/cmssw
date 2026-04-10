#ifndef OmtfPhase2AngleConverter_h
#define OmtfPhase2AngleConverter_h

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfAngleConverter.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTThContainer.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

class OmtfPhase2AngleConverter : public OmtfAngleConverter {
public:
  OmtfPhase2AngleConverter() : OmtfAngleConverter() {}
  ~OmtfPhase2AngleConverter() override = default;

  // Convert DT phi to OMTF coordinate system.
  int getProcessorPhi(int phiZero, l1t::tftype part, int dtScNum, int dtPhi) const override;

  int getGlobalEta(DTChamberId dTChamberId, const L1Phase2MuDTThContainer *dtThDigis, int bxNum) const;

  /// When true, getGlobalEta returns fixed mid-chamber values (firmware export mode).
  /// When false (default), the real theta-digi LUT lookup is used (sample production mode).
  void setDtFixedPointEtaForFirmware(bool v) { dtFixedPointEtaForFirmware_ = v; }
  bool getDtFixedPointEtaForFirmware() const { return dtFixedPointEtaForFirmware_; }

private:
  bool dtFixedPointEtaForFirmware_ = false;
};

#endif

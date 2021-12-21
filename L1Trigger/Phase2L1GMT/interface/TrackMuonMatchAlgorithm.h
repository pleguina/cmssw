#ifndef PHASE2GMT_TRACKMUONMATCHALGO
#define PHASE2GMT_TRACKMUONMATCHALGO

#include "L1Trigger/Phase2L1GMT/interface/ApSignAbsInt.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1Trigger/interface/L1TObjComparison.h"
#include "L1Trigger/Phase2L1GMT/interface/TrackConverter.h"
#include "L1Trigger/Phase2L1GMT/interface/PreTrackMatchedMuon.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"


namespace Phase2L1GMT {

  const unsigned int PHIDIVIDER = 1 << (BITSPHI - BITSSTUBCOORD);
  const unsigned int ETADIVIDER = 1 << (BITSETA - BITSSTUBETA);

  typedef struct {
    ap_int<BITSSTUBCOORD> coord1;
    ap_uint<BITSSIGMACOORD> sigma_coord1;
    ap_int<BITSSTUBCOORD> coord2;
    ap_uint<BITSSIGMACOORD> sigma_coord2;
    ap_int<BITSSTUBETA> eta;
    ap_uint<BITSSIGMAETA> sigma_eta1;
    ap_uint<BITSSIGMAETA> sigma_eta2;
    ap_uint<1> valid;
    ap_uint<1> is_barrel;
  } propagation_t;

  typedef struct {
    ap_uint<BITSMATCHQUALITY - 2> quality;
    ap_uint<BITSSTUBID> id;
    ap_uint<1> valid;
    bool isGlobal;
    l1t::RegionalMuonCandRef muRef;
    l1t::MuonStubRef stubRef;

    //ap_uint<BITSSIGMACOORD + 1> deltaCoord1;
    //ap_uint<BITSSIGMACOORD + 1> deltaCoord2;

    ApSignAbsInt<BITSSIGMACOORD> deltaCoord1;
    ApSignAbsInt<BITSSIGMACOORD> deltaCoord2;

    ApSignAbsInt<BITSSIGMAETA> deltaEta1;
    ApSignAbsInt<BITSSIGMAETA> deltaEta2;

  } match_t;

  class TrackMuonMatchAlgorithm {
  public:
    TrackMuonMatchAlgorithm(const edm::ParameterSet& iConfig) : verbose_(iConfig.getParameter<int>("verbose")) {}

    ~TrackMuonMatchAlgorithm() {}

    std::vector<PreTrackMatchedMuon> processNonant(const std::vector<ConvertedTTTrack>& convertedTracks,
                                                   const std::vector<MuonROI>& rois) {
      LogTrace("phase2L1GMT")<<"processNonant: rois.size(): "<<rois.size();
      for(auto& roi : rois) {
        LogTrace("phase2L1GMT")<<" roi ";
        for(auto& stub : roi.stubs()) {
          stub->print();
        }
      }

      std::vector<PreTrackMatchedMuon> preMuons;
      for (const auto& track : convertedTracks) {
        PreTrackMatchedMuon mu = processTrack(track, rois);
        if (mu.valid() && preMuons.size() < MATCHER_OUT_BUFFER_SIZE)
          preMuons.push_back(mu);
      }
      std::vector<PreTrackMatchedMuon> cleanedMuons = clean(preMuons);
      return cleanedMuons;
    }

    std::vector<PreTrackMatchedMuon> cleanNeighbor(std::vector<PreTrackMatchedMuon>& muons,
                                                   std::vector<PreTrackMatchedMuon>& muonsPrevious,
                                                   std::vector<PreTrackMatchedMuon>& muonsNext,
                                                   bool equality) {
      std::vector<PreTrackMatchedMuon> out;

      if (muons.empty())
        return out;

      if (verbose_ == 1) {
        printf("-----Cleaning Up Muons in the neighbours\n");
        printf("Before:\n");
      }

      LogTrace("phase2L1GMT")<<"\ncleanNeighbor";
      for (uint i = 0; i < muons.size(); ++i) {
        LogTrace("phase2L1GMT")<<"cleanNeighbor for "<<muons[i];

        ap_uint<10> mask = 0x3ff; //TODO use constant
        for (uint j = 0; j < muonsPrevious.size(); ++j) {
          mask = mask & cleanMuon(muons[i], muonsPrevious[j], equality);
        }
        for (uint j = 0; j < muonsNext.size(); ++j) {
          mask = mask & cleanMuon(muons[i], muonsNext[j], equality);
        }
        if (mask) {
          LogTrace("phase2L1GMT")<<"alive";
          out.push_back(muons[i]);
        } else {
          LogTrace("phase2L1GMT")<<"killed";
        }
      }
      return out;
    }

    std::vector<l1t::TrackerMuon> convert(std::vector<PreTrackMatchedMuon>& muons, uint maximum) {
      std::vector<l1t::TrackerMuon> out;
      for (const auto& mu : muons) {
        if (out.size() == maximum)
          break;
        l1t::TrackerMuon muon(mu.trkPtr(), mu.charge(), mu.pt(), mu.eta(), mu.phi(), mu.z0(), mu.d0(), mu.quality());
        //muon.setMuonRef(mu.muonRef());
        for (const auto& stub : mu.stubs())
          muon.addStub(stub);
        out.push_back(muon);
        if (verbose_ == 1) {
          printf("Final Muon:");
          muon.print();
        }
        //LogTrace("phase2L1GMT")<<"\mconvert ";
      }
      return out;
    }

    bool outputGT(std::vector<l1t::TrackerMuon>& muons) {
      for (auto& mu : muons) {
        wordtype word1 = 0;
        wordtype word2 = 0;

        int bstart = 0;
        bstart = wordconcat<wordtype>(word1, bstart, mu.hwPt(), BITSGTPT);
        bstart = wordconcat<wordtype>(word1, bstart, mu.hwPhi(), BITSGTPHI);
        bstart = wordconcat<wordtype>(word1, bstart, mu.hwEta(), BITSGTETA);
        bstart = wordconcat<wordtype>(word1, bstart, mu.hwZ0(), BITSGTZ0);
        bstart = wordconcat<wordtype>(word1, bstart, (mu.hwD0() >> 2), BITSGTD0);

        bstart = 0;
        bstart = wordconcat<wordtype>(word2, bstart, mu.hwCharge(), 1);
        bstart = wordconcat<wordtype>(word2, bstart, mu.hwQual(), BITSGTQUALITY);
        bstart = wordconcat<wordtype>(word2, bstart, mu.hwIso(), BITSGTISO);
        bstart = wordconcat<wordtype>(word2, bstart, mu.hwBeta(), BITSMUONBETA);

        std::array<uint64_t, 2> wordout = {{word1, word2}};
        mu.setWord(wordout);
      }
      return true;
    }

    std::vector<l1t::TrackerMuon> sort(std::vector<l1t::TrackerMuon>& muons, uint maximum) {
      if (muons.size() < 2)
        return muons;

      std::sort(muons.begin(), muons.end(), [](l1t::TrackerMuon a, l1t::TrackerMuon b) { return a.hwPt() > b.hwPt(); });
      std::vector<l1t::TrackerMuon> out;
      for (unsigned int i = 0; i < muons.size(); ++i) {
        out.push_back(muons[i]);
        if (i == (maximum - 1))
          break;
      }

      return out;
    }

  private:
    int verbose_;

    propagation_t propagate(const ConvertedTTTrack& track, uint layer) {
      ap_uint<BITSPROPCOORD> prop_coord1 = 0;
      ap_uint<BITSPROPCOORD> prop_coord2 = 0;
      ap_uint<BITSPROPSIGMACOORD_A> res0_coord1 = 0;
      ap_uint<BITSPROPSIGMACOORD_B> res1_coord1 = 0;
      ap_uint<BITSPROPSIGMACOORD_A> res0_coord2 = 0;
      ap_uint<BITSPROPSIGMACOORD_B> res1_coord2 = 0;
      ap_uint<BITSPROPSIGMAETA_A> res0_eta1 = 0;
      ap_uint<BITSPROPSIGMAETA_B> res1_eta = 0;
      ap_uint<BITSPROPSIGMAETA_A> res0_eta2 = 0;
      ap_uint<1> is_barrel = 0;

      uint reducedAbsEta = track.abseta() / 8; //looks good here
      //abseta_ is 3138 for eta 2.40676 so 12 bits, track.eta is 13 bits

      if (layer == 0) {
        prop_coord1 = lt_prop_coord1_0[reducedAbsEta];
        prop_coord2 = lt_prop_coord2_0[reducedAbsEta];
        res0_coord1 = lt_res0_coord1_0[reducedAbsEta];
        res1_coord1 = lt_res1_coord1_0[reducedAbsEta];
        res0_coord2 = lt_res0_coord2_0[reducedAbsEta];
        res1_coord2 = lt_res1_coord2_0[reducedAbsEta];
        res0_eta1 = lt_res0_eta1_0[reducedAbsEta];
        res1_eta = lt_res1_eta_0[reducedAbsEta];
        res0_eta2 = lt_res0_eta2_0[reducedAbsEta];
        is_barrel = reducedAbsEta < barrelLimit0_ ? 1 : 0; //so barrelLimit*_ are not good
      } else if (layer == 1) {
        prop_coord1 = lt_prop_coord1_1[reducedAbsEta];
        prop_coord2 = lt_prop_coord2_1[reducedAbsEta];
        res0_coord1 = lt_res0_coord1_1[reducedAbsEta];
        res1_coord1 = lt_res1_coord1_1[reducedAbsEta];
        res0_coord2 = lt_res0_coord2_1[reducedAbsEta];
        res1_coord2 = lt_res1_coord2_1[reducedAbsEta];
        res0_eta1 = lt_res0_eta1_1[reducedAbsEta];
        res1_eta = lt_res1_eta_1[reducedAbsEta];
        res0_eta2 = lt_res0_eta2_1[reducedAbsEta];
        is_barrel = reducedAbsEta < barrelLimit1_ ? 1 : 0;

      } else if (layer == 2) {
        prop_coord1 = lt_prop_coord1_2[reducedAbsEta];
        prop_coord2 = lt_prop_coord2_2[reducedAbsEta];
        res0_coord1 = lt_res0_coord1_2[reducedAbsEta];
        res1_coord1 = lt_res1_coord1_2[reducedAbsEta];
        res0_coord2 = lt_res0_coord2_2[reducedAbsEta];
        res1_coord2 = lt_res1_coord2_2[reducedAbsEta];
        res0_eta1 = lt_res0_eta1_2[reducedAbsEta];
        res1_eta = lt_res1_eta_2[reducedAbsEta];
        res0_eta2 = lt_res0_eta2_2[reducedAbsEta];
        is_barrel = reducedAbsEta < barrelLimit2_ ? 1 : 0;

      } else if (layer == 3) {
        prop_coord1 = lt_prop_coord1_3[reducedAbsEta];
        prop_coord2 = lt_prop_coord2_3[reducedAbsEta];
        res0_coord1 = lt_res0_coord1_3[reducedAbsEta];
        res1_coord1 = lt_res1_coord1_3[reducedAbsEta];
        res0_coord2 = lt_res0_coord2_3[reducedAbsEta];
        res1_coord2 = lt_res1_coord2_3[reducedAbsEta];
        res0_eta1 = lt_res0_eta1_3[reducedAbsEta];
        res1_eta = lt_res1_eta_3[reducedAbsEta];
        res0_eta2 = lt_res0_eta2_3[reducedAbsEta];
        is_barrel = reducedAbsEta < barrelLimit3_ ? 1 : 0;

      } else if (layer == 4) {
        prop_coord1 = lt_prop_coord1_4[reducedAbsEta];
        prop_coord2 = lt_prop_coord2_4[reducedAbsEta];
        res0_coord1 = lt_res0_coord1_4[reducedAbsEta];
        res1_coord1 = lt_res1_coord1_4[reducedAbsEta];
        res0_coord2 = lt_res0_coord2_4[reducedAbsEta];
        res1_coord2 = lt_res1_coord2_4[reducedAbsEta];
        res0_eta1 = lt_res0_eta1_4[reducedAbsEta];
        res1_eta = lt_res1_eta_4[reducedAbsEta];
        res0_eta2 = lt_res0_eta2_4[reducedAbsEta];
        is_barrel = 0;
      }

      propagation_t out;
      ap_int<BITSTTCURV> curvature = track.curvature();
      ap_int<BITSPHI> phi = track.phi();
      ap_int<BITSPROPCOORD + BITSTTCURV> c1kFull = prop_coord1 * curvature;
      ap_int<BITSPROPCOORD + BITSTTCURV - 10> c1k = (c1kFull) / 1024;
      ap_int<BITSPHI> coord1 = phi - c1k;

      out.coord1 = coord1 / PHIDIVIDER; //PHIDIVIDER is 2^5 = 32

      ap_int<BITSPROPCOORD + BITSTTCURV> c2kFull = prop_coord2 * curvature;

      ap_int<BITSPROPCOORD + BITSTTCURV - 10> c2k = (c2kFull) / 1024;
      if (is_barrel)
        out.coord2 = -c2k / PHIDIVIDER;
      else
        out.coord2 = (phi - c2k) / PHIDIVIDER;

      ap_int<BITSETA> eta = track.eta();
      out.eta = eta / ETADIVIDER; //ETADIVIDER = 32

      ap_uint<2 * BITSTTCURV - 2> curvature2All = curvature * curvature;
      ap_uint<BITSTTCURV2> curvature2 = curvature2All / 2;

      ap_ufixed<3, 2, AP_TRN_ZERO, AP_SAT_SYM> sigmaFactor(SIGMA_FACTOR);

      //Remember to change emulator with new k2
      ap_uint<BITSPROPSIGMACOORD_B + BITSTTCURV2> rescoord1k = (res1_coord1 * curvature2) >> 23;
      ap_ufixed<BITSSIGMACOORD, BITSSIGMACOORD, AP_TRN_ZERO, AP_SAT_SYM> sigma_coord1 = (res0_coord1 + rescoord1k) * sigmaFactor; //TODO<<<<<<<<<<<<<<
      out.sigma_coord1 = ap_uint<BITSSIGMACOORD>(sigma_coord1);

      ap_uint<BITSPROPSIGMACOORD_B + BITSTTCURV2> rescoord2k = (res1_coord2 * curvature2) >> 23;
      ap_ufixed<BITSSIGMACOORD, BITSSIGMACOORD, AP_TRN_ZERO, AP_SAT_SYM> sigma_coord2 = (res0_coord2 + rescoord2k) * sigmaFactor; //TODO<<<<<<<<<<<<<<
      out.sigma_coord2 = ap_uint<BITSSIGMACOORD>(sigma_coord2);

      ap_uint<BITSPROPSIGMAETA_B + BITSTTCURV2> resetak = (res1_eta * curvature2) >> 23;
      ap_ufixed<BITSSIGMAETA, BITSSIGMAETA, AP_TRN_ZERO, AP_SAT_SYM> sigma_eta1 = res0_eta1 + resetak;
      out.sigma_eta1 = ap_uint<BITSSIGMAETA>(sigma_eta1);
      ap_ufixed<BITSSIGMAETA, BITSSIGMAETA, AP_TRN_ZERO, AP_SAT_SYM> sigma_eta2 = res0_eta2 + resetak;
      out.sigma_eta2 = ap_uint<BITSSIGMAETA>(sigma_eta2);
      out.valid = 1;
      out.is_barrel = is_barrel;

      LogTrace("phase2L1GMT")<<"\n\nTrackMuonMatchAlgorithm::propagate: "
          <<" layer "<<layer
          <<" is_barrel "<<out.is_barrel.to_int()
          <<" coord1 "<<out.coord1.to_int()
          <<" sigma_coord1 "<<out.sigma_coord1.to_int()
          <<" coord2 "<<out.coord2.to_int()
          <<" sigma_coord2 "<<out.sigma_coord2.to_int()
          <<" eta "<< out.eta.to_int()
          <<" sigma_eta1 "<< out.sigma_eta1.to_int()
          <<" sigma_eta2 "<< out.sigma_eta2.to_int();

      return out;
    }

    ap_uint<BITSSIGMAETA + 1> deltaEta(const ap_int<BITSSTUBETA>& eta1, const ap_int<BITSSTUBETA>& eta2) {
      ap_fixed<BITSSIGMAETA + 2, BITSSIGMAETA + 2, AP_TRN_ZERO, AP_SAT_SYM> dEta = eta1 - eta2;
      if (dEta < 0)
        return ap_uint<BITSSIGMAETA + 1>(-dEta);
      else
        return ap_uint<BITSSIGMAETA + 1>(dEta);
    }

    ap_uint<BITSSIGMACOORD + 1> deltaCoord(const ap_int<BITSSTUBCOORD>& phi1, const ap_int<BITSSTUBCOORD>& phi2) {
      //dPhiRoll has the same width as  phi1 and phi2 on purpose,
      //because its overflowing provides folding of dPhi in around the +-180deg
      ap_int<BITSSTUBCOORD> dPhiRoll = phi1 - phi2;

     /* ap_int<BITSSTUBCOORD+1> dPhiRoll_test= phi1 - phi2;
      if(dPhiRoll != dPhiRoll_test) {
        edm::LogError("gmtDataDumper")<<__FUNCTION__<<":"<<__LINE__
            <<"dPhiRoll overflow! dPhiRoll: "
            <<dPhiRoll<<" = "<<(dPhiRoll * 2. /256.*180.)<<" dPhiRoll_test "<<dPhiRoll_test<<" = "<<(dPhiRoll_test * 2. /256.*180.)
            <<" phi1 "<<phi1<<" = "<<(phi1 * 2. /256.*180.)<<" deg phi2 "<<phi2<<" = "<<(phi2 * 2. /256.*180.);
      }*/

      ap_ufixed<BITSSIGMACOORD + 1, BITSSIGMACOORD + 1, AP_TRN_ZERO, AP_SAT_SYM> dPhi;
      if (dPhiRoll < 0)
        dPhi = ap_ufixed<BITSSIGMACOORD + 1, BITSSIGMACOORD + 1, AP_TRN_ZERO, AP_SAT_SYM>(-dPhiRoll);
      else
        dPhi = ap_ufixed<BITSSIGMACOORD + 1, BITSSIGMACOORD + 1, AP_TRN_ZERO, AP_SAT_SYM>(dPhiRoll);

      return ap_uint<BITSSIGMACOORD + 1>(dPhi);
    }

    match_t match(const propagation_t prop, const l1t::MuonStubRef& stub) {
      LogTrace("phase2L1GMT")<<"\nTrackMuonMatchAlgorithm::match: Matching to ";
      stub->print();

      //Matching of Coord1
      ap_uint<1> coord1Matched;
      ap_uint<BITSSIGMACOORD + 1> deltaCoord1 = deltaCoord(prop.coord1, stub->coord1());
      if (deltaCoord1 <= prop.sigma_coord1 && (stub->quality() & 0x1)) {
        coord1Matched = 1;
      } else {
        coord1Matched = 0;
      }

      LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match:"
          <<" coord1Matched "<<coord1Matched.to_int()
          <<" deltaCoord1 "<<deltaCoord1.to_int()
          <<" prop.sigma_coord1 "<<prop.sigma_coord1.to_int()
          <<" stub->quality "<<stub->quality() ;

      //Matching of Coord2
      ap_uint<1> coord2Matched;
      ap_uint<BITSSIGMACOORD + 1> deltaCoord2 = deltaCoord(prop.coord2, stub->coord2());
      if (deltaCoord2 <= prop.sigma_coord2 && (stub->quality() & 0x2)) {
        coord2Matched = 1;
      } else {
        coord2Matched = 0;
      }

      LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match "
          <<" coord2Matched "<<coord2Matched.to_int()
          <<" deltaCoord2 "<<deltaCoord2.to_int()
          <<" prop.sigma_coord2 "<<prop.sigma_coord2.to_int();
      //Matching of Eta1

      ap_uint<1> eta1Matched;

      //if we have really bad quality[Barrel no eta]
      //increase the resolution
      ap_ufixed<BITSSIGMAETA, BITSSIGMAETA, AP_TRN_ZERO, AP_SAT_SYM> prop_sigma_eta1;
      if (stub->etaQuality() == 0) {
        prop_sigma_eta1 = prop.sigma_eta1 + 6; //TODO check how saturation here impact the performance
      }
      else
        prop_sigma_eta1 = prop.sigma_eta1;

      //Eta Quality encoding:
      //0 - no ete measurement, the coarseEta* is given (middle of the chamber)
      //1 - first eta measurement is present
      //3 - first and second eta measurement are present
      ap_uint<BITSSIGMAETA + 1> deltaEta1 = deltaEta(prop.eta, stub->eta1());
      if (deltaEta1 <= prop_sigma_eta1 && (stub->etaQuality() == 0 || (stub->etaQuality() & 0x1)))
        eta1Matched = 1;
      else
        eta1Matched = 0;

      if (verbose_ == 1)
        printf("eta1 matched=%d delta=%d res=%d\n", eta1Matched.to_int(), deltaEta1.to_int(), prop_sigma_eta1.to_int());

      LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match:"
          <<" eta1Matched "<<eta1Matched.to_int()
          <<" deltaEta1 "<<deltaEta1.to_int()
          <<" prop.prop_sigma_eta1 "<<prop_sigma_eta1.to_int();
      //Matching of Eta2

      ap_uint<1> eta2Matched;

      ap_uint<BITSSIGMAETA + 1> deltaEta2 = deltaEta(prop.eta, stub->eta2());
      if (deltaEta2 <= prop.sigma_eta2 && (stub->etaQuality() & 0x2))
        eta2Matched = 1;
      else
        eta2Matched = 0;
      match_t out;
      out.id = stub->id();

      if (verbose_ == 1)
        printf("eta2 matched=%d delta=%d res=%d\n", eta2Matched.to_int(), deltaEta2.to_int(), prop.sigma_eta2.to_int());

      LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match:"
          <<" eta2Matched "<<eta2Matched.to_int()
          <<" deltaEta2 "<<deltaEta2.to_int()
          <<" prop.prop_sigma_eta2 "<<prop.sigma_eta2.to_int();

      //absValue = 0 and sign = 0 means no match
      //absValue = 0 and sign = 1 means matched with delta = 0
      out.deltaCoord1.absValue = 0;
      out.deltaCoord1.sign = 0;
      out.deltaCoord2.absValue = 0;
      out.deltaCoord2.sign = 0;

      out.deltaEta1.absValue = 0;
      out.deltaEta1.sign = 0;
      out.deltaEta2.absValue = 0;
      out.deltaEta2.sign = 0;

      unsigned int maxQuality = 1 << (BITSSIGMACOORD + 1);

      //if barrel, coord1 has to always be matched, coord2 maybe and eta1 is needed if etaQ=0 or then the one that depends on eta quality
      if (prop.is_barrel) {
        out.valid = (coord1Matched == 1 && (eta1Matched == 1 || eta2Matched == 1));
        if (out.valid == 0) {
          out.quality = 0;
        } else {
          out.quality = maxQuality - deltaCoord1;
          if (coord2Matched == 1)
            out.quality += maxQuality - deltaCoord2;
          //fixme KB. out.quality is 7 bits, deltaCoord1 is 5 bits, so it looks that 7 bits for the quality is one bit too much

          //out.deltaCoord1 overflows are not possible, because coord1Matched = deltaCoord1 <= prop.sigma_coord1

          /*
          if(prop.coord1 >= stub->coord1())
            out.deltaCoord1 = deltaCoord1Middle + deltaCoord1;
          else
            out.deltaCoord1 = deltaCoord1Middle - deltaCoord1;

          if (coord2Matched == 1) {
            if(prop.coord2 >= stub->coord2())
              out.deltaCoord2 = deltaCoord1Middle + deltaCoord2;
            else
              out.deltaCoord2 = deltaCoord1Middle - deltaCoord2;
          }*/

          out.deltaCoord1.absValue = deltaCoord1;
          out.deltaCoord1.sign = (prop.coord1 >= stub->coord1());  //sign is minus if prop.coord1 >= stub->coord1()

          if (coord2Matched == 1) {
            out.deltaCoord2.absValue = deltaCoord2;
            out.deltaCoord2.sign = (prop.coord2 >= stub->coord2() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
          }

          if(eta1Matched == 1) {
            if(stub->etaQuality() == 0) { //no eta hit, just the chamber is matched
              out.deltaEta1.absValue = (1 << out.deltaEta1.absValue.width) -1 ; //max possible value
              out.deltaEta1.sign = 1;
            }
            else {
              out.deltaEta1.absValue = deltaEta1;
              out.deltaEta1.sign = (prop.eta >= stub->eta1() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
            }
          }

          if(eta2Matched == 1) {
            out.deltaEta2.absValue = deltaEta2;
            out.deltaEta2.sign = (prop.eta >= stub->eta2() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
          }
          if(stub->type() == 0)
            LogTrace("phase2L1GMT")<<" prop.is_barrel but stub.type() == 0 !!!!!!!!!!!!!!!!!!!!!!";
          LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match barrel: "
              <<"out.deltaCoord1 "<<std::setw(3)<<out.deltaCoord1.to_int()<<" sign-abs: "<<out.deltaCoord1.toRaw().to_string()<<" sign "<<out.deltaCoord1.sign.to_int()<<" absVal "<<out.deltaCoord1.absValue.to_int()<<" "
              <<"out.deltaCoord2 "<<std::setw(3)<<out.deltaCoord2.to_int()<<" sign-abs: "<<out.deltaCoord2.toRaw().to_string()<<" sign "<<out.deltaCoord2.sign.to_int()<<" absVal "<<out.deltaCoord2.absValue.to_int()
              <<" deltaEta1 "<<out.deltaEta1.to_int() <<" deltaEta2 "<<out.deltaEta2.to_int();
        }
      }
      //if endcap each coordinate is independent except the case where phiQuality=1 and etaQuality==3
      //KB but phiQuality=1 and etaQuality==3 is not possible???
      else { //endcap
        bool match1 = (coord1Matched == 1 && eta1Matched == 1);
        bool match2 = (coord2Matched == 1 && eta2Matched == 1);
        bool match3 =
            (coord1Matched == 1 && (eta1Matched || eta2Matched) && stub->etaQuality() == 3 && stub->quality() == 1);
        //etaQuality = 3 - first and second eta measurement are present
        //stub->quality() == 1 - CSCOnlyStub <<<<
        //stub->quality() == 2 - RPCOnlyStub
        //stub->quality() == 3 - CSC and RPC stub

        LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match endcap: match1 "<< match1<<" match2 "<< match2<<" match3 "<<match3<<endl;

        out.valid = match1 || match2 || match3;
        if (out.valid == 0)
          out.quality = 0;
        else {
          out.quality = 0;
          if (match1 || match3)
            out.quality += maxQuality - deltaCoord1;
          if (match2)
            out.quality += maxQuality - deltaCoord2;

          if (match1 || match3) {
            out.deltaCoord1.absValue = deltaCoord1;
            out.deltaCoord1.sign = (prop.coord1 >= stub->coord1());  //sign is minus if prop.coord1 >= stub->coord1()

            if(eta1Matched == 1) {
              out.deltaEta1.absValue = deltaEta1;
              out.deltaEta1.sign = (prop.eta >= stub->eta1() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
            }
            else if(eta1Matched == 2) {
              //in some rare cases, if match3 the two etas are little different, and then  one can be matched, and the other no
              //this case should be removed in the L1TPhase2GMTEndcapStubProcessor
              //but for the moment we take the deltaEta2 and put into the out.deltaEta1
              out.deltaEta1.absValue = deltaEta2;
              out.deltaEta1.sign = (prop.eta >= stub->eta2() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
            }
          }

          if (match2) {
            out.deltaCoord2.absValue = deltaCoord2;
            out.deltaCoord2.sign = (prop.coord2 >= stub->coord2() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()

            out.deltaEta2.absValue = deltaEta2;
            out.deltaEta2.sign = (prop.eta >= stub->eta1() ? 1 : 0);  //sign is minus if prop.coord1 >= stub->coord2()
          }

          if(stub->type() == 1)
            LogTrace("phase2L1GMT")<<" prop.is_endcap but stub.type() == 1 !!!!!!!!!!!!!!!!!!!!!!";
          LogTrace("phase2L1GMT")<<"TrackMuonMatchAlgorithm::match endcap: "
              <<"out.deltaCoord1 "<<std::setw(3)<<out.deltaCoord1.to_int()<<" sign_abs: "<<out.deltaCoord1.toRaw().to_string()<<" "
              <<"out.deltaCoord2 "<<std::setw(3)<<out.deltaCoord2.to_int()<<" sign_abs: "<<out.deltaCoord2.toRaw().to_string()
              <<" deltaEta1 "<<out.deltaEta1.to_int() <<" deltaEta2 "<<out.deltaEta2.to_int();
        }
      }
      if (verbose_ == 1)
        printf("GlobalMatchQuality = %d\n", out.quality.to_int());
      out.stubRef = stub;
      return out;
    }

    match_t propagateAndMatch(const ConvertedTTTrack& track, const l1t::MuonStubRef& stub) {
      propagation_t prop = propagate(track, stub->tfLayer());
      return match(prop, stub);
    }

    match_t getBest(const std::vector<match_t>& matches) {
      auto* best = &(matches[0]); //using pointer to avoid copying
      for (const auto& m : matches) {
        if (m.quality > best->quality)
          best = &m;
      }

      return *best;
    }

    PreTrackMatchedMuon processTrack(const ConvertedTTTrack& track, const std::vector<MuonROI>& rois) {
      std::vector<match_t> matchInfo0;
      std::vector<match_t> matchInfo1;
      std::vector<match_t> matchInfo2;
      std::vector<match_t> matchInfo3;
      std::vector<match_t> matchInfo4;

      LogTrace("phase2L1GMT")<<"\nprocessTrack: -----------------";
      LogTrace("phase2L1GMT")<<track;

      for (const auto& roi : rois) {
        if (verbose_ == 1) {
          printf("New ROI with %d stubs \n", int(roi.stubs().size()));
        }
        for (const auto& stub : roi.stubs()) {
          match_t m = propagateAndMatch(track, stub);
          if (m.valid == 1) {
            
            if (roi.isGlobalMuon() && roi.muonRef().isNonnull()) {
              m.isGlobal = true;
              m.muRef = roi.muonRef();
            }

            if (stub->tfLayer() == 0)
              matchInfo0.push_back(m);
            else if (stub->tfLayer() == 1)
              matchInfo1.push_back(m);
            else if (stub->tfLayer() == 2)
              matchInfo2.push_back(m);
            else if (stub->tfLayer() == 3)
              matchInfo3.push_back(m);
            else if (stub->tfLayer() == 4)
              matchInfo4.push_back(m);

          }
        }
      }

      ap_ufixed<6, 6, AP_TRN_ZERO, AP_SAT_SYM> ptPenalty = ap_ufixed<6, 6, AP_TRN_ZERO, AP_SAT_SYM>(track.pt() / 32);

      ap_uint<BITSMATCHQUALITY> quality = 0;
      PreTrackMatchedMuon muon(track.curvature(), track.charge(), track.pt(), track.eta(), track.phi(), track.z0(), track.d0());

      if (!matchInfo0.empty()) {
        match_t b = getBest(matchInfo0);
        if (b.valid) {
          muon.addStub(b.stubRef);
          muon.getDeltaCoords1()[0] = b.deltaCoord1;
          muon.getDeltaCoords2()[0] = b.deltaCoord2;

          muon.getDeltaEtas1()[0] = b.deltaEta1;
          muon.getDeltaEtas2()[0] = b.deltaEta2;

          if (b.isGlobal)
            muon.addMuonRef(b.muRef);
          quality += b.quality;
        }
      }
      if (!matchInfo1.empty()) {
        match_t b = getBest(matchInfo1);
        if (b.valid) {
          muon.addStub(b.stubRef);
          muon.getDeltaCoords1()[1] = b.deltaCoord1;
          muon.getDeltaCoords2()[1] = b.deltaCoord2;

          muon.getDeltaEtas1()[1] = b.deltaEta1;
          muon.getDeltaEtas2()[1] = b.deltaEta2;

          if (b.isGlobal)
            muon.addMuonRef(b.muRef);
          quality += b.quality;
        }
      }
      if (!matchInfo2.empty()) {
        match_t b = getBest(matchInfo2);
        muon.getDeltaCoords1()[2] = b.deltaCoord1;
        muon.getDeltaCoords2()[2] = b.deltaCoord2;

        muon.getDeltaEtas1()[2] = b.deltaEta1;
        muon.getDeltaEtas2()[2] = b.deltaEta2;

        if (b.valid) {
          muon.addStub(b.stubRef);
          if (b.isGlobal)
            muon.addMuonRef(b.muRef);
          quality += b.quality;
        }
      }
      if (!matchInfo3.empty()) {
        match_t b = getBest(matchInfo3);
        muon.getDeltaCoords1()[3] = b.deltaCoord1;
        muon.getDeltaCoords2()[3] = b.deltaCoord2;

        muon.getDeltaEtas1()[3] = b.deltaEta1;
        muon.getDeltaEtas2()[3] = b.deltaEta2;

        if (b.valid) {
          muon.addStub(b.stubRef);
          if (b.isGlobal)
            muon.addMuonRef(b.muRef);
          quality += b.quality;
        }
      }
      if (!matchInfo4.empty()) {
        match_t b = getBest(matchInfo4);
        muon.getDeltaCoords1()[4] = b.deltaCoord1;
        muon.getDeltaCoords2()[4] = b.deltaCoord2;

        muon.getDeltaEtas1()[4] = b.deltaEta1;
        muon.getDeltaEtas2()[4] = b.deltaEta2;

        if (b.valid) {
          muon.addStub(b.stubRef);
          if (b.isGlobal)
            muon.addMuonRef(b.muRef);
          quality += b.quality;
        }
      }

      muon.setOfflineQuantities(track.offline_pt(), track.offline_eta(), track.offline_phi());
      muon.setTrkPtr(track.trkPtr());

      ap_uint<8> etaAddr = muon.eta() < 0 ? ap_uint<8>(-muon.eta() / 256) : ap_uint<8>((muon.eta()) / 256);
      ap_uint<8> ptAddr = muon.pt() > 4095 ? ap_uint<8>(15) : ap_uint<8>(muon.pt() / 256);
      ap_uint<8> addr = ptAddr | (etaAddr << 4);
      ap_uint<8> qualityCut = lt_tpsID[addr];

      if (quality >= (qualityCut - QUALITY_CUT_MOD)) { //TODO changing quality cut to looser!!!!!!!!!!!!!!!!!!!!!!!111
        muon.setValid(true);
        muon.setQuality(quality + ptPenalty);
      } else {
        muon.setValid(false);
        muon.setQuality(0);
        muon.resetGlobal();
      }

      if (verbose_ == 1 && !rois.empty()) {  //patterns for HLS

        printf("TPS %d", track.trkPtr()->phiSector());
        track.printWord();

        for (uint i = 0; i < 16; ++i) {
          if (rois.size() > i) {
            rois[i].printROILine();
          } else {
            printf("%08x", 0);
            printf("%016lx", 0x1ff000000000000);
            printf("%016lx", 0x1ff000000000000);
            printf("%016lx", 0x1ff000000000000);
            printf("%016lx", 0x1ff000000000000);
            printf("%016lx", 0x1ff000000000000);
          }
        }
        muon.printWord();
        printf("\n");
      }
      return muon;
    }

//original cleanMuon
    ap_uint<5> cleanMuonV0(const PreTrackMatchedMuon& mu, const PreTrackMatchedMuon& other, bool eq) {
      ap_uint<5> valid = 0;
      ap_uint<5> overlap = 0;
      if (mu.stubID0() != 511) {
        valid = valid | 0x1;
        if (mu.stubID0() == other.stubID0())
          overlap = overlap | 0x1;
      }
      if (mu.stubID1() != 511) {
        valid = valid | 0x2;
        if (mu.stubID1() == other.stubID1())
          overlap = overlap | 0x2;
      }
      if (mu.stubID2() != 511) {
        valid = valid | 0x4;
        if (mu.stubID2() == other.stubID2())
          overlap = overlap | 0x4;
      }
      if (mu.stubID3() != 511) {
        valid = valid | 0x8;
        if (mu.stubID3() == other.stubID3())
          overlap = overlap | 0x8;
      }
      if (mu.stubID4() != 511) {
        valid = valid | 0x10;
        if (mu.stubID4() == other.stubID4())
          overlap = overlap | 0x10;
      }

      if (((mu.quality() < other.quality()) && (!eq)) || ((mu.quality() <= other.quality()) && (eq)))
        return valid & (~overlap);
      else
        return valid;
    }

    /*
     * TODO
     * this method leaves littl bit more duplicates when cleanin between the nonants
     * because the duplicated ttTracks sometimes have little bit different pt, eta or phi
     * and then the deltaCoords are not exactly equal
     * TODO
     * would be better to clean the cross nonant duplicated just comparing pt eta and phi
     * all with appropriate margins that must be found from simulation
     */
    ap_uint<10> cleanMuonV1(PreTrackMatchedMuon& mu, PreTrackMatchedMuon& other, bool eq) {
      ap_uint<10> valid = 0;

      LogTrace("phase2L1GMT")<<"other    muon: "<<other;
      for(unsigned int layer = 0 ; layer < 5; layer++) {
        LogTrace("phase2L1GMT")<<"layer "<<layer<<" this  stubID "<<std::setw(3)<<mu.stubID(layer)
                <<" deltaCoords1 "<<std::setw(3)<<mu.getDeltaCoords1()[layer].to_int()<<" sign "<<mu.getDeltaCoords1()[layer].sign
                <<" deltaCoords2 "<<std::setw(3)<<mu.getDeltaCoords2()[layer].to_int()<<" sign "<<mu.getDeltaCoords2()[layer].sign;
        LogTrace("phase2L1GMT")<<"layer "<<layer<<" other stubID "<<std::setw(3)<<other.stubID(layer)
                <<" deltaCoords1 "<<std::setw(3)<<other.getDeltaCoords1()[layer].to_int()<<" sign "<<mu.getDeltaCoords1()[layer].sign
                <<" deltaCoords2 "<<std::setw(3)<<other.getDeltaCoords2()[layer].to_int()<<" sign "<<mu.getDeltaCoords2()[layer].sign<<"\n";

        if (mu.stubID(layer) != 511) {
          unsigned int coord1Idx = 2 * layer;
          unsigned int coord2Idx = coord1Idx + 1;

          //absValue = 0 and sign = 0 means no match, so then valid = 0
          if( !(mu.getDeltaCoords1()[layer].absValue == 0 && mu.getDeltaCoords1()[layer].sign == 0) )
            valid = valid | (1 << coord1Idx);

          if( !(mu.getDeltaCoords2()[layer].absValue == 0 && mu.getDeltaCoords2()[layer].sign == 0) )
            valid = valid | (1 << coord2Idx);

          if (mu.stubID(layer) == other.stubID(layer)) {
            bool otherIsValid1 = !(other.getDeltaCoords1()[layer].absValue == 0 && other.getDeltaCoords1()[layer].sign == 0);
            if( (!eq && (otherIsValid1 && mu.getDeltaCoords1()[layer].absValue >  other.getDeltaCoords1()[layer].absValue)) ||
                ( eq && (otherIsValid1 && mu.getDeltaCoords1()[layer].absValue >= other.getDeltaCoords1()[layer].absValue))   ) {
              valid = valid & ~(1 << coord1Idx);
            }

            bool otherIsValid2 = !(other.getDeltaCoords2()[layer].absValue == 0 && other.getDeltaCoords2()[layer].sign == 0);
            if( (!eq && (otherIsValid2 && mu.getDeltaCoords2()[layer].absValue >  other.getDeltaCoords2()[layer].absValue)) ||
                ( eq && (otherIsValid2 && mu.getDeltaCoords2()[layer].absValue >= other.getDeltaCoords2()[layer].absValue))    ) {
              valid = valid & ~(1 << coord2Idx);
            }
          }
        }
      }

      LogTrace("phase2L1GMT")<<"valid "<<valid.to_string();

      return valid;
    }

    ap_uint<10> cleanMuon(PreTrackMatchedMuon& mu, PreTrackMatchedMuon& other, bool eq) {
      if(CLEANMUON_VERSION == 1)
        return cleanMuonV1(mu, other, eq) ;

      return cleanMuonV0(mu, other, eq) ;
    }

    std::vector<PreTrackMatchedMuon> clean(std::vector<PreTrackMatchedMuon>& muons) {
      std::vector<PreTrackMatchedMuon> out;
      if (muons.empty())
        return out;

      LogTrace("phase2L1GMT")<<"\nclean: ";
      for (uint i = 0; i < muons.size(); ++i) {
        ap_uint<10> mask = 0x3ff;
        LogTrace("phase2L1GMT")<<"cleaning muon: "<<muons[i];
        for (uint j = 0; j < muons.size(); ++j) {
          if (i == j)
            continue;
          mask = mask & cleanMuon(muons[i], muons[j], false);
        }
        if (mask) {
          out.push_back(muons[i]);
        }
        LogTrace("phase2L1GMT")<<"final mask : "<<mask.to_string()<<" = "<<(mask == 0 ? "killed" : "alive")<<"\n";
      }
      return out;
    }
  };
}  // namespace Phase2L1GMT

#endif

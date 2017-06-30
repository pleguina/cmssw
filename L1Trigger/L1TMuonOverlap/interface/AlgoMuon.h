#ifndef AlgoMuon_H
#define AlgoMuon_H

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include <L1Trigger/L1TMuonOverlap/interface/IGoldenPattern.h>
#include <ostream>

class AlgoMuon{

 public:
  
  // AlgoMuon() : pt(-1.), eta(99.), phi(9999.), disc(-999), bx(0), q(-1), charge(99), refLayer(-1), hits(0) {} // the old one version 
  AlgoMuon() : m_disc(-999), m_phi(9999), m_eta(99), m_refLayer(-1), m_firedLayerBits(0), m_q(-1), m_bx(0), m_pt(-1), m_charge(99) {}
/*  AlgoMuon(int disc=-999, int phi=9999, int eta=99, int refLayer=-1,
              int hits=0, int q=-1, int bx=0, int pt=-1, int charge=99):
              m_disc(disc), m_phi(phi), m_eta(eta), m_refLayer(refLayer), 
              m_firedLayerBits(hits), m_q(q), m_bx(bx), m_pt(pt), m_charge(charge), 
              m_patNumb(999), m_rhitNumb(999) {}*/

  AlgoMuon(const GoldenPatternResult& gpResult, const Key& gpKey, unsigned int refHitNumber):
              m_disc(gpResult.getPdfWeigtSum()), m_phi(gpResult.getPhi()), m_eta(gpResult.getEta()),
              m_refLayer(gpResult.getRefLayer()),
              m_firedLayerBits(gpResult.getFiredLayerBits()), m_q(gpResult.getFiredLayerCnt()),
              m_bx(0), m_pt(gpKey.thePt), m_charge(gpKey.theCharge),
              m_phiRHit(gpResult.getRefHitPhi()),
              m_patNumb(gpKey.theNumber), m_rhitNumb(refHitNumber) {}

  int getDisc() const { return m_disc; }
  int getPhi()  const { return m_phi; }
  int getEta()  const { return m_eta; }
  int getRefLayer() const { return m_refLayer; }
  int getFiredLayerBits() const { return m_firedLayerBits; }
  int getQ()  const { return m_q; }
  int getBx() const { return m_bx; }
  int getPt() const { return m_pt; }
  int getCharge()   const { return m_charge; }
  int getPhiRHit()  const { return m_phiRHit; }
  unsigned int getPatternNumber() const { return m_patNumb; }
  unsigned int getRefHitNumber() const { return m_rhitNumb; }

  void setDisc(int disc) { m_disc = disc; }
  void setPhi(int phi)   { m_phi = phi; }
  void setEta(int eta)   { m_eta = eta; }
  void setRefLayer(int refLayer) { m_refLayer = refLayer; }
  void setHits(int hits) { m_firedLayerBits = hits; }
  void setQ(int q)    { m_q = q; }
  void setBx(int bx)  { m_bx = bx; }
  void setPt(int pt)  { m_pt = pt; }
  void setCharge(int charge)   { m_charge = charge; }
  void setPhiRHit(int phiRHit) { m_phiRHit = phiRHit; }
  void setPatternNumber(unsigned int aPatNum) { m_patNumb = aPatNum;}
  void setRefHitNumber(unsigned int aRefHitNum) { m_rhitNumb = aRefHitNum; }

  bool isValid() const;  

  bool operator< (const AlgoMuon & o) const;

  friend std::ostream & operator<< (std::ostream &out, const AlgoMuon &o);

 private: 

  int m_disc;
  int m_phi;
  int m_eta;
  int m_refLayer;
  int m_firedLayerBits;
  int m_q;
  int m_bx; 
  int m_pt;
  int m_charge;
  int m_phiRHit;
  // to add 
  // int m_pdf; 
  unsigned int m_patNumb;
  unsigned int m_rhitNumb;

};
#endif

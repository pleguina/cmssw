#ifndef AlgoMuon_H
#define AlgoMuon_H

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h>
#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include <ostream>

class AlgoMuon {

 public:
  AlgoMuon() {}
/*  AlgoMuon(int disc=-999, int phi=9999, int eta=99, int refLayer=-1,
              int hits=0, int q=-1, int bx=0, int pt=-1, int charge=99):
              m_disc(disc), m_phi(phi), m_eta(eta), m_refLayer(refLayer), 
              m_firedLayerBits(hits), m_q(q), m_bx(bx), m_pt(pt), m_charge(charge), 
              m_patNumb(999), m_rhitNumb(999) {}*/

  AlgoMuon(const GoldenPatternResult& gpResult, GoldenPatternBase* gp, unsigned int refHitNumber, int bx = 0):
              gpResult(gpResult), goldenPatern(gp),
              m_q(gpResult.getFiredLayerCnt()), //initial value of quality, can be altered later
              m_bx(bx), m_rhitNumb(refHitNumber) {}

/*  AlgoMuon(const GoldenPatternResult& gpResult, const Key& gpKey, unsigned int refHitNumber):
              goldenPatern(0),
              m_disc(gpResult.getPdfWeigtSum()), m_phi(gpResult.getPhi()), m_eta(gpResult.getEta()),
              m_refLayer(gpResult.getRefLayer()),
              m_firedLayerBits(gpResult.getFiredLayerBits()), m_q(gpResult.getFiredLayerCnt()),
              m_bx(0), m_pt(gpKey.thePt), m_charge(gpKey.theCharge),
              m_phiRHit(gpResult.getRefHitPhi()),
              m_patNumb(gpKey.theNumber), m_rhitNumb(refHitNumber) {}*/

  GoldenPatternBase* getGoldenPatern() const {
    return goldenPatern;
  }

  virtual ~AlgoMuon() {}

  const GoldenPatternResult& getGpResult() const {
    return gpResult;
  }

  omtfPdfValueType getDisc() const { return gpResult.getPdfSum(); }
  int getPhi()  const { return gpResult.getPhi(); }
  int getEta()  const { return gpResult.getEta(); }
  int getRefLayer() const { return gpResult.getRefLayer(); }
  int getFiredLayerBits() const { return gpResult.getFiredLayerBits(); }
  int getQ()  const { return m_q; }
  int getBx() const { return m_bx; }

  int getPt() const {
    if(goldenPatern ==  nullptr)
      return -1;
    return goldenPatern->key().thePt;
  }

  int getCharge()   const {
    if(goldenPatern ==  nullptr)
      return 0;
    return goldenPatern->key().theCharge;
  }
  int getPhiRHit()  const { return gpResult.getRefHitPhi(); }

  unsigned int getPatternNumber() const {
    if(goldenPatern ==  nullptr)
      return 0;
    return goldenPatern->key().theNumber;
  }

  unsigned int getRefHitNumber() const { return m_rhitNumb; }

  void setQ(int q)    { m_q = q; }
  void setEta(int eta)   { gpResult.setEta(eta); }

/*  void setDisc(omtfPdfValueType disc) { m_disc = disc; }
  void setPhi(int phi)   { m_phi = phi; }

  void setRefLayer(int refLayer) { m_refLayer = refLayer; }
  void setHits(int hits) { m_firedLayerBits = hits; }

  void setBx(int bx)  { m_bx = bx; }
  void setPt(int pt)  { m_pt = pt; }
  void setCharge(int charge)   { m_charge = charge; }
  void setPhiRHit(int phiRHit) { m_phiRHit = phiRHit; }
  void setPatternNumber(unsigned int aPatNum) { m_patNumb = aPatNum;}*/

  void setRefHitNumber(unsigned int aRefHitNum) { m_rhitNumb = aRefHitNum; }

  bool isValid() const;

  const bool isKilled() const {
    return killed;
  }

  void kill() {
    killed = true;
    //FIXME maybe also valid = false???
  }

  bool operator< (const AlgoMuon & o) const;

  friend std::ostream & operator<< (std::ostream &out, const AlgoMuon &o);

 private: 
  ///gpResult stored in the goldenPatern can be updated more then once in one event (many ttTracks in one event in one processor!),
  ///so gpResult cannot be a reference or pointer but a copy
  GoldenPatternResult gpResult;

  GoldenPatternBase* goldenPatern = nullptr;

  int m_q = -1;
  int m_bx  = 0;

/*  omtfPdfValueType m_disc;
  int m_phi;
  int m_eta;
  int m_refLayer;
  int m_firedLayerBits;


  int m_pt;
  int m_charge;
  int m_phiRHit;
  // to add 
  // int m_pdf; 
  unsigned int m_patNumb;*/
  unsigned int m_rhitNumb = 0;

  bool killed = false;
};

typedef std::vector<std::shared_ptr<AlgoMuon> > AlgoMuons;

#endif

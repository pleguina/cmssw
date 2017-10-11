#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include <sstream>
#include <string>
#include <vector>

using namespace std;
void drawHitDistr(int ipt, int charge, int layer, int refLayer, int color);

TCanvas* c0 = new TCanvas("histo", "histo", 800, 600);
TLegend* legend = new TLegend(0.7,0.4,0.9,0.6);

int compareDistr(){
  cout << "ROOT" << endl;
  c0->SetGrid();

  int layer = 7;
  int refLayer = 3;

  drawHitDistr(2*50.0 + 1, 1, layer, refLayer, kBlack);
  drawHitDistr(2*20  + 1, 1, layer, refLayer, kBlue);
  drawHitDistr(2*18 + 1, 1, layer, refLayer, kRed);
  drawHitDistr(2*4.5 + 1, 1, layer, refLayer, kGreen);
  drawHitDistr(2*10.0 + 1, 1, layer, refLayer, kMagenta);


  legend->Draw();

  return 0;

}

bool first = true;
void drawHitDistr(int iPt, int charge, int layer, int refLayer, int color) {
  string path = "/afs/cern.ch/work/k/kpijanow/public/CMSSW_9_2_0/src/L1Trigger/L1TMuonOverlap/test/expert/";
  std::ostringstream s1;
  s1 << path << "bendingDistr_ptCode_" << iPt << "_ch_" << charge << ".root";

  if(access( (s1.str()).c_str(), F_OK ) == -1) {
    cout<<"cannot open the file "<<s1.str()<<endl;
    return;;
  }

  TFile* file = new TFile((TString)s1.str());
  TObject* hist = file->Get(Form("ipt_%i_ch%i_layer_%i_refLayer_%i", iPt, charge, layer, refLayer) );
  if( dynamic_cast<TH1F*>(hist) ) {
    TH1F* h = dynamic_cast<TH1F*>(hist);
    Double_t scale = 1./(h->Integral());
    h->Scale(scale);
    if(first)
      h->Draw("hist");
    else
      h->Draw("histsame");

    h->SetMarkerColor(color);
    h->SetLineColor(color);

    std::ostringstream ostr;
    ostr<<" iPt "<<iPt<<" = "<<( (iPt-1)/2.)<<" GeV";
    legend->AddEntry(h, ostr.str().c_str(), "l");
  }
  else if( dynamic_cast<TH2F*>(hist) ) {
    TH2F* h = dynamic_cast<TH2F*>(hist);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    string opt = ""; //"lego"; //"cont1";//box";
    if(first)
      h->Draw("box");
    else
      h->Draw((opt+"same").c_str() );

    std::ostringstream ostr;
    ostr<<" iPt "<<iPt<<" = "<<( (iPt-1)/2.)<<" GeV";
    legend->AddEntry(h, ostr.str().c_str(), "l");
  }

  first = false;
}


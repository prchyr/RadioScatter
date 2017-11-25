#ifndef RS_EVENT_H
#define RS_EVENT_H
#include "TObject.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TVirtualFFT.h"
#include <deque>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/LorentzVector.h"




// class TH1F;
// class TTree;
// class TGraph;
// class TObject;

using namespace CLHEP;
using namespace std;

class RadioScatterEvent:public TObject
{
private:
  //default histogram for storing spectra, etc. 
  TH1F* spectrumHist=new TH1F("spectrumHist", "spectrumHist", 100, 0, 100);
  TH2F* spectrogramHist=new TH2F("spectrogramHist", "spectrogramHist", 100, 0, 100, 100, 0, 100);
  Int_t dummy;
  TRandom *ran = new TRandom();
public:
  RadioScatterEvent();


  Hep3Vector direction;
  Hep3Vector position;
  Hep3Vector tx;
  Hep3Vector rx;
  double primaryEnergy=0;
  double sampleRate;
  double nPrimaries=0;
  double txVoltage=0;
  double freq=0;
  //  double power=0;
  double totNScatterers=0;
  TH1F * eventHist=0;
  TH1F * reHist=0;
  TH1F * imHist=0;

  TGraph * eventGraph=0;


  //plotting things

  TGraph  getComplexEnvelope(double cutoff=0);
  TGraph  getLowpassFiltered(double cutoff);
  TGraph  getGraph();
  //  TGraph getSpectrum(bool dbflag=false);
  TH1F * getSpectrum(bool dbflag=false);  
  void spectrogram(Int_t binsize = 128, Int_t overlap=32);
  int plotEvent(int bins=64, int overlap=8);
  
  double chirpSlope();
  double startFreq();
  double stopFreq();
  double peakFreq();
  double bandWidth();
  double peakV();
  double rms();
  double duration();
  double power();
  int triggered(double thresh);
  
  ClassDef(RadioScatterEvent, 2);
};
#endif

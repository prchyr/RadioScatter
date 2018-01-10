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
//#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"




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
  HepLorentzVector *tx;
  HepLorentzVector *rx;
  //the energy of the primary as set in geant
  double primaryEnergy=0;

  double sampleRate;
  double nPrimaries=0;
  double txVoltage=0;
  double freq=0;
  //  double power=0;
  double totNScatterers=0;
  vector<vector<TH1F*>> eventHist;
  vector<vector<TH1F*>> reHist;
  vector<vector<TH1F*>> imHist;

  vector<vector<TGraph*>> eventGraph;

  int ntx=1;
  int nrx=1;
  //plotting things

  TGraph  getComplexEnvelope(int txindex, int rxindex,double cutoff=0);
  TGraph  getLowpassFiltered(int txindex, int rxindex,double cutoff);
  TGraph  getGraph();
  //  TGraph getSpectrum(bool dbflag=false);
  TH1F * getSpectrum(int txindex, int rxindex,bool dbflag=false);  
  void spectrogram(int txindex, int rxindex,Int_t binsize = 128, Int_t overlap=32);
  int plotEvent(int txindex, int rxindex,int bins=64, int overlap=8);
  
  double chirpSlope();
  double startFreq();
  double stopFreq();
  double peakFreq(int txindex, int rxindex);
  double bandWidth();
  double peakV(int txindex, int rxindex);
  double rms(int txindex, int rxindex);
  double duration(int txindex, int rxindex);
  double power(int txindex, int rxindex);
  double pathLength(int txindex, int rxindex);
  //energy calculated from the geant4 energy and the number of primaries
  double primaryParticleEnergy();
  int triggered(double thresh);
  
  ClassDef(RadioScatterEvent, 3);
};
#endif

/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */
#ifndef RS_EVENT_H
#define RS_EVENT_H
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TLegend.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include <deque>

//#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

//#include <gsl/gsl_linalg.h>


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
  TCanvas *ccc=0;
  TH1F* spectrumHist=new TH1F("spectrumHist", "spectrumHist", 100, 0, 100);
  TH1F* ceHist = new TH1F("cehist", "cehist", 100, 1, -1);//complex envelope
  TH2F* spectrogramHist=new TH2F("spectrogramHist", "spectrogramHist", 100, 0, 100, 100, 0, 100);
  TH3F* rxhist = new TH3F("rxhist", "rxhist", 100, 1, -1, 100, 1, -1, 100, 1, -1);
  TH3F* txhist = new TH3F("txhist", "txhist", 100, 1, -1, 100, 1, -1, 100, 1, -1);
  TH3F* vertexhist = new TH3F("vertexhist", "vertexhist", 100, 1, -1, 100, 1, -1, 100, 1, -1);
  TH3F* pointingHist = new TH3F("pointingHist", "pointingHist", 100, 1, -1, 100, 1, -1, 100, 1, -1);
  TH3F* triggeredhist = new TH3F("triggeredhist", "triggeredhist", 100, 1, -1, 100, 1, -1, 100, 1, -1);
  Int_t dummy;
  TRandom *ran = new TRandom();

  double dt[64000][3];
  HepLorentzVector source[64000];
  int POINTING_MAP_BUILT=0;
  int RXHIST_FILLED=0;  
public:
  RadioScatterEvent();


  Hep3Vector direction=Hep3Vector(0,0,1);
  Hep3Vector position=Hep3Vector(1, 1, 1);
  //HepLorentzVector *tx;
  //HepLorentzVector *rx;
  vector<HepLorentzVector> tx;
  vector<HepLorentzVector> rx;
  //the energy of the primary as set in geant
  double primaryEnergy=1000;//1000MeV (geant units) as a default

  double sampleRate;
  double nPrimaries=0;
  double txVoltage=0;
  double txPowerW=0;
  double freq=0;
  double txGain=1.;
  double rxGain=1.;
  //  double power=0;
  double totNScatterers=0;
  vector<vector<TH1F*>> eventHist;
  vector<vector<TH1F*>> reHist;
  vector<vector<TH1F*>> imHist;
  vector<TH1F*> testHist;

  vector<vector<TGraph*>> eventGraph;

  int SINE_SUBTRACT=0;
  int ntx=1;
  int nrx=1;
  //plotting things

  TH1F *  getComplexEnvelope(int txindex, int rxindex,double cutoff=0);
  TGraph  getLowpassFiltered(int txindex, int rxindex,double cutoff);
  TGraph * getGraph(int txindex, int rxindex);
  //  TGraph getSpectrum(bool dbflag=false);
  TH1F * getSpectrum(int txindex, int rxindex,bool dbflag=false);  
  void spectrogram(int txindex, int rxindex,Int_t binsize = 128, Int_t overlap=32);
  int plotEvent(int txindex, int rxindex, double noise_flag=0, int show_geom=0, int bins=64, int overlap=8);
  int plotEventNotebook(int txindex, int rxindex, int noise_flag=0, int show_geom=0, int bins=64, int overlap=8);
  int reset();
  double thermalNoiseRMS();
  double chirpSlope();
  double startFreq();
  double stopFreq();
  double sineSubtract(int txindex, int rxindex, double rangestart=0, double rangeend=240, double p0=.02, double p1=1., double p2=0.);
  int backgroundSubtract(int txindex, int rxindex, TH1F *bSubHist);
  TH1F* makeBackgroundSubtractHist(int txindex, int rxindex, TString bfile);
  double peakFreq(int txindex, int rxindex);
  double bandWidth();
  double peakV(int txindex, int rxindex);
  double effectiveCrossSection(int txindex, int rxindex);
  double rms(int txindex, int rxindex);
  double duration(int txindex, int rxindex);
  double integratedPower(int txindex, int rxindex);
  double integratedPower(int txindex, int rxindex, double tlow, double thigh, double dcoffset=0.);
  double integratedPowerAroundPeak(int txindex, int rxindex, double window=100);
  double integratedVoltage(int txindex, int rxindex);
  double integratedVoltage(int txindex, int rxindex, double tlow, double thigh, double dcoffset=0.);
  double peakPowerMW(int txindex, int rxindex);
  double peakPowerW(int txindex, int rxindex);
  double pathLengthM(int txindex, int rxindex);
  double pathLengthMM(int txindex, int rxindex);
  //energy calculated from the geant4 energy and the number of primaries
  double primaryParticleEnergy();
  int triggered(double thresh, int n_antennas=1);
  int nTriggered(double thresh);
  int trigSingle(double thresh, int ant=0);

  //pointing things

  //HepLorentzVector findSource(int debug=0);
  // HepLorentzVector pointingTest();
  int buildMap();
  HepLorentzVector pointUsingMap();



  ClassDef(RadioScatterEvent, 5);
};
#endif

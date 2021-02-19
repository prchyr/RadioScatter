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
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TUtilRadioScatter.hh"
//#include <vector>

//#include "CLHEP/Units/PhysicalConstants.h"
//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/LorentzVector.h"

//#include <gsl/gsl_linalg.h>


// class TH1F;
// class TTree;
// class TGraph;
// class TObject;
using namespace TUtilRadioScatter;
//using namespace CLHEP;
//using namespace std;
/**\brief This is the storage class for a RadioScatter object, called an event. an event is a single scatter from a cascade, as detected in all of the receivers. 

The RadioScatterEvent is a TObject, meaning that it is easily stored in a ROOT tree and is plottable from ttree->Draw(). It has many member variables that store information about the cascade, the geometry of the receiver(s) and transmitter(s), and some basic variables about the events that can be easiliy plotted. But it also stores the raw waveforms captured in each receiver, which can be used in analysis. 

*/

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
  TLorentzVector source[64000];
  int POINTING_MAP_BUILT=0;
  int RXHIST_FILLED=0;  
public:
  RadioScatterEvent();


  TVector3 direction=TVector3(0,0,1);
  TVector3 position=TVector3(1, 1, 1);
  TVector3 polarization=TVector3(0,0,1);//antenna polarization
  //TLorentzVector *tx;
  //TLorentzVector *rx;
  std::vector<TLorentzVector> tx;
  std::vector<TLorentzVector> rx;
  //the energy of the primary as set in geant
  double primaryEnergy=1*GeV;//1GeV as a default
  double targetEnergy=1*GeV;///1GeV
  double inelasticity=1.;
  double weight=1.; //used for some calculations. 
  double sampleRate;
  double nPrimaries=0;
  double txVoltage=0;
  double txPowerW=0;
  double freq=0;
  double txGain=1.;
  double rxGain=1.;
  std::vector<std::vector<double>> angRVP;///angle subtended by receiver, vertex, and cascade momentum direction
  std::vector<std::vector<double>> angTVP;///angle subtended by transmitter, vertex, and cascade momentum direction
  std::vector<std::vector<double>> angTVR;///angle subtended by transmitter, vertex, and receiver.
  std::vector<std::vector<double>> beta;//the bistatic angle
  std::vector<std::vector<double>> delta;//the aspect angle
  std::vector<std::vector<double>> doppler;//the doppler shift
  //  double power=0;
  double totNScatterers=0;
  std::vector<std::vector<TH1F*>> eventHist;
  std::vector<std::vector<TH1F*>> reHist;
  std::vector<std::vector<TH1F*>> imHist;
  std::vector<TH1F*> testHist;

  std::vector<std::vector<TGraph*>> eventGraph;

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
  int plotEvent(int txindex, int rxindex, double noise_flag=0, int show_geom=0, int bins=256, int overlap=128, int logFlag=2);
  int plotEventNotebook(int txindex, int rxindex, double noise_flag=0, int show_geom=0, int bins=64, int overlap=8, int logFlag=2);
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

  double duration(int txindex, int rxindex, double highThreshRatio=.3, double lowThreshRatio=.1);///the duration of a pulse. expressed in terms of a threshold as a ratio of the peak voltage. so here the high threshold is .3*peakV, and the low threshold is .1*peakV.
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

  //TLorentzVector findSource(int debug=0);
  // TLorentzVector pointingTest();
  int buildMap();
  TLorentzVector pointUsingMap();



  ClassDef(RadioScatterEvent, 5);
};
#endif

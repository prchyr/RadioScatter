/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */
#ifndef RS_EVENT_SUM_H
#define RS_EVENT_SUM_H
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
#include <deque>
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

//using namespace CLHEP;
using namespace std;
/**\brief This is a storage class of summary variables of a RadioScatterEvent, that are useful for plotting and calculations. 

The summary class is useful for things like effective volume calculations, where a large number of events are run, and the storage space for the individual waveform data for each event and each antenna would be prohibitive. This allows for some useful variables to be calculated for each event to be stored in this object for analysis. 

An RSEventSummary is created when setMakeSummary(1) is called from a radioscatter object

 */
class RSEventSummary:public TObject
{
public:

  RSEventSummary(int ntransmitters=1, int nreceivers=1);
  virtual~RSEventSummary();

  TVector3 direction;
  TVector3 position;
  TVector3 polarization;

  vector<TLorentzVector> tx{1};
  vector<TLorentzVector> rx{1};
  //the energy of the primary as set in geant
  double primaryEnergyG4=0;

  double sampleRate=0;
  double nPrimaries=0;
  double txVoltageV=0;
  double txPowerW=0;
  double freq=0;
  double weight=1.;
  //  double power=0;
  double totNScatterers=0;
  double primaryParticleEnergy=0;
  double inelasticity=1;
  int ntx=1;
  int nrx=1;
  //plotting things

  vector<vector<double>> peakFreq;
  vector<vector<double>> peakV;
  vector<vector<double>> effectiveCrossSection;
  vector<vector<double>> rms;
  vector<vector<double>> duration;
  vector<vector<double>> integratedPower;
  vector<vector<double>> peakPowerW;
  vector<vector<double>> pathLengthM;
  std::vector<std::vector<double>> angRVP;///angle subtended by receiver, vertex, and cascade momentum direction
  std::vector<std::vector<double>> angTVP;///angle subtended by transmitter, vertex, and cascade momentum direction
  std::vector<std::vector<double>> angTVR;///angle subtended by transmitter, vertex, and receiver.
  std::vector<std::vector<double>> beta;///the bistatic angle
  std::vector<std::vector<double>> delta;///the aspect angle
  std::vector<std::vector<double>> doppler;
  

  int triggered(double thresh, int n_antennas=1);
  int nTriggered(double thresh);
  int trigSingle(double thresh, int txx=0, int rxx=0);

  
  ClassDef(RSEventSummary, 1);
};
#endif

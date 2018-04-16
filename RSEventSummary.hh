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

class RSEventSummary:public TObject
{
public:

  RSEventSummary(int ntransmitters=1, int nreceivers=1);
  virtual~RSEventSummary();

  Hep3Vector direction;
  Hep3Vector position;

  vector<HepLorentzVector> tx{1};
  vector<HepLorentzVector> rx{1};
  //the energy of the primary as set in geant
  double primaryEnergyG4=0;

  double sampleRate=0;
  double nPrimaries=0;
  double txVoltageV=0;
  double txPowerW=0;
  double freq=0;
  //  double power=0;
  double totNScatterers=0;
  double primaryParticleEnergy=0;

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
  

  int triggered(double thresh, int n_antennas=1);
  int nTriggered(double thresh);
  int trigSingle(double thresh, int txx=0, int rxx=0);

  
  ClassDef(RSEventSummary, 1);
};
#endif

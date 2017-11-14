#ifndef EVENT_TREE_H
#define EVENT_TREE_H
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
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

class EventTree:public TObject
{
private:
  Int_t dummy;
public:
  EventTree();
  //  ~eventTree();

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
  TGraph getSpectrum(bool dbflag=false);  
  
  double chirpSlope();
  double startFreq();
  double stopFreq();
  double bandWidth();
  double peakV();
  double rms();
  double duration();
  double power();
  int triggered(double thresh);
  
  ClassDef(EventTree, 1);
};
#endif

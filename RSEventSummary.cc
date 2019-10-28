/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */
#include "RSEventSummary.hh"


ClassImp(RSEventSummary)

RSEventSummary::RSEventSummary(int ntransmitters, int nreceivers){
  ntx=ntransmitters;
  nrx=nreceivers;
  tx.resize(ntx);
  rx.resize(nrx);
  peakFreq.resize(ntx, vector<double>(nrx, 0.));
  peakV.resize(ntx, vector<double>(nrx, 0.));
  effectiveCrossSection.resize(ntx, vector<double>(nrx, 0.));
  rms.resize(ntx, vector<double>(nrx, 0.));
  duration.resize(ntx, vector<double>(nrx, 0.));
  integratedPower.resize(ntx, vector<double>(nrx, 0.));
  peakPowerW.resize(ntx, vector<double>(nrx, 0.));
  pathLengthM.resize(ntx, vector<double>(nrx, 0.));
}

RSEventSummary::~RSEventSummary(){
}

//simple thing-did any receiver trigger?
int RSEventSummary::triggered(double thresh, int n_antennas){
  int trig=0,num=0;
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      trig=peakV[i][j]>=thresh?1:0;
      if(trig==1)num++;
      if(num>=n_antennas)return 1;
    }
  }
  return trig;
}
//how many triggered?
int RSEventSummary::nTriggered(double thresh){
  int trig=0,num=0;
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      trig=peakV[i][j]>=thresh?1:0;
      if(trig==1)num++;
    }
  }
  return num;
}

int RSEventSummary::trigSingle(double thresh, int txx, int rxx){
 return peakV[txx][rxx]>=thresh?1:0;
}


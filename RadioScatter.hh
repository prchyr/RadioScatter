/*
todo:
antenna effective height
antenna effective area
(currently these are both set to 1 accross all frequencies)

 
 */
#ifndef R_Scat
#define R_Scat


//#include "errno.h"
#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "EventTree.hh"

 using namespace CLHEP;
 using namespace std;


// #ifdef __CINT__
// #pragma link C++ class Hep3Vector;
// #endif


class RadioScatter{


  EventTree event;


  TString output_file_name;
  char* pol = "horizontal";//default, also set as default in set_simulation_paramaters
  TString polarization="horizontal";
  double x_offset=0*m, z_offset = 2.*m, y_offset = 5.*m;
  double omega, omega_c, omega_e = 7.79e20, omega_0=twopi*2000*megahertz;
  double nu_col = 0;
  double tcs = .655e-24;//thompson cross section
  double n_primaries = 1;//set this based on number of events in run
  double tx_voltage = 1.;//mW
  double impedance = 50;
  double tx_gain=1.;
  double rx_gain=1.;
  double step_length=0;//default g4 steplength
  double E_i=.000038;//default electron ion pair energy
  double e_charge_cgs = 4.803e-10;//statcoloumbs
  double e_mass_cgs = 9.109e-28;//g
  //TFile *outfile;
  //  TFile * outfile=new TFile("/home/natas/Documents/physics/geant/root/time.root", "RECREATE");
  TH1F *fft_hist, *power_hist;
  TH1F *time_hist = new TH1F("eventHist", "eventHist", 100, 0, 10);
  TH1F *re_hist = new TH1F("reHist", "reHist", 100, 0, 10);
  TH1F *im_hist = new TH1F("imHist", "imHist", 100, 0, 10);
  //e^2/m_e
  double plasma_const = sqrt(4*pi*electron_charge*electron_charge/electron_mass_c2);

  //e^2/(4pi epislon0 m c^2), dividing by meters puts things in terms of meters, not mm as is default. makes distance calculations more accurate
  double cross_section=classic_electr_radius/m;

  double half_window = 300;//number of nanoseconds in 1/2 of the record window. can be changed;   

  int useAttnLengthFlag=0;
  double attnLength=460*m;
  //  ofstream of;
  
  double frequency, period, lambda, k;

  double samplerate=10, samplingperiod=.1, start_time=0, end_time=1000;

  double txp, tx_on=-999999999., tx_off=999999999.;

  int includeCW_flag=0;//whether to simulate the direct signal as well.
  
  vector<double> amplitudes, timeofarrival, phases, field, plasma; 

  //variables for refraction manipulation
  
  double k_r, c_light_r, mag1, mag2;//, tof, txphase, kx;
  //distance from the antennas to the interface, must be set by user
  double tx_interface_dist,rx_interface_dist;
  //relative index of refraction, calculated to always be >1.
  double n_rel=1.5;

  
  double testvalue;

public:


  HepLorentzVector tx;//transmitter
  HepLorentzVector rx;//reciever


  
  RadioScatter();
  
  void makeOutputFile(TString filename);
  void makeOutputTextFile(char* filename);
  void writeToTextFile();
  void makeTimeHist();
  void setTxPos(double xin, double yin, double zin);
  void setRxPos(double xin, double yin, double zin);
  void setTxPos(Hep3Vector in);
  void setRxPos(Hep3Vector in);
  void setTxFreq(double f);
  void setTxVoltage(double v);
  void setNPrimaries(double n);
  void setPrimaryEnergy(double e);
  void setPolarization(char * p);
  void setTxVals(double f, double power, double gain);
  void setRxVals(double s, double gain);
  void setSimulationParameters(double n, char* tx_rx_pol, double relative_index_of_refraction, int flag);
  void setRelativeIndexOfRefraction(double iof);
  void setCalculateUsingAttnLength(int val=0);
  void setRecordWindowLength(double nanoseconds);
  void setRxSampleRate(double rate);
  void setTxInterfaceDistX(double dist);
  void setRxInterfaceDistX(double dist);
  void setShowCWFlag(double i);
  void setTxOnTime(double on);
  void setTxOffTime(double off);
  Hep3Vector getTxPos();
  Hep3Vector getRxPos();
  double getFreq();
  double getTxGain(double angle);
  double getRxGain(double angle);
  double getTxAmplitude(HepLorentzVector point);
  double getRxAmplitude(HepLorentzVector point, Hep3Vector j1,  Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);
  double getRxAmplitude(HepLorentzVector point);
  double getTxPhase(double t_0);
  double getRxTime(HepLorentzVector point, Hep3Vector j, Hep3Vector l);
  double getRxTime(HepLorentzVector point);
  double getTxTime(HepLorentzVector point, int direct);
  double getRxPhase(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);
  double getRxPhase(HepLorentzVector point);
  double getRxPhaseRel(HepLorentzVector point, double v);
  int checkTxOn(double time);
  Hep3Vector findRefractionPlane(HepLorentzVector p1, HepLorentzVector p2);
  double findRefractionPointZ(double kx, double kz, double jx);
  double findPathLengthWithRefraction(HepLorentzVector p1, HepLorentzVector p2, double interface_dist_x); 
  TH1F* getDirectSignal(const TH1F*);
  double getDirectSignalPhase(HepLorentzVector point);
  double makeRays(HepLorentzVector point, double e, double l, double e_i);
  double makeRays(HepLorentzVector point, double e, double l, double e_i, TH1F *hist);
  double power();
  void draw();
  void printEventStats();
  void writeHist(TString filename, float num_events);
  TH1F * scaleHist(float num_events);
  int writeRun(float num_events=1., int debug=0);
  void writeEvent(TString filename, float num_events);
  void close();


  void testFunc(double testVal);


private:

  int REFRACTION_FLAG=0;
  bool RSCAT_HIST_DECLARE=false;
  int fRunCounter=0;
};

#endif

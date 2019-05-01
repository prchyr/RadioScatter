/*
todo:
antenna effective height

currently we use a variant of the RICE formula
h=sqrt(lam^2 (gain)/480pi^2) where we set the gain=480pi^2. essentially, 
we say that this is a magic antenna where the effective height is 
equal to the wavelength at all freqs. in fact, this is reasonable for 
a nice broadband antenna.
 
 */
#ifndef R_Scat
#define R_Scat


//#include "errno.h"
#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "RadioScatterEvent.hh"
#include "RSEventSummary.hh"

 using namespace CLHEP;
 using namespace std;


// #ifdef __CINT__
// #pragma link C++ class Hep3Vector;
// #endif


class RadioScatter{
public:

  RadioScatterEvent event;


  TString output_file_name;
  char* pol = "horizontal";//default, also set as default in set_simulation_paramaters
  TString polarization="horizontal";
  double x_offset=0*m, z_offset = 2.*m, y_offset = 5.*m;
  double omega, omega_c, omega_e = 7.79e20, omega_0=twopi*2000*megahertz;
  double nu_col = 0;
  double tcs = .655e-24;//thompson cross section
  double n_primaries = 1;//set this based on number of events in run
  double tx_voltage = 1.;//mW
  double zscale=1.;//the longitudinal scale factor
  double tscale=1.;//the time scale factor
  double impedance = 50;
  double tx_gain=1.;
  double rx_gain=1.;
  double step_length=1;//default g4 steplength
  double E_i=.000038;//default electron ion pair energy
  double e_charge_cgs = 4.803e-10;//statcoloumbs
  double e_mass_cgs = 9.109e-28;//g
  double k_b = 1.38e-23;//j/k
  //TFile *outfile;
  //  TFile * outfile=new TFile("/home/natas/Documents/physics/geant/root/time.root", "RECREATE");
  TH1F *fft_hist, *power_hist;
  TH1F *t_h = new TH1F("eventHist", "eventHist", 100, 0, 10);
  TH1F *re_h = new TH1F("reHist", "reHist", 100, 0, 10);
  TH1F *im_h = new TH1F("imHist", "imHist", 100, 0, 10);
  TGraph *event_gr = new TGraph();
  
  vector<vector<TH1F*>> time_hist;
  vector<vector<TH1F*>> re_hist;
  vector<vector<TH1F*>> im_hist;
  vector<vector<TGraph*>> event_graph;
  //e^2/m_e
  double plasma_const = sqrt(4*pi*electron_charge*electron_charge/electron_mass_c2);

  //e^2/(4pi epislon0 m c^2), in units of m
  double e_radius=classic_electr_radius/m;

  double half_window = 300;//number of nanoseconds in 1/2 of the record window. can be changed;   

  int useAttnLengthFlag=0;
  //  double attnLength=460;//meters
  double attnLength=1400.;//dzb length for upper 1.5km at pole
  //  ofstream of;
  
  double frequency, period, lambda, k;

  double effectiveHeight=1.;//, antennaGain=1.;
  
  double samplerate=10, samplingperiod=.1, start_time=0, end_time=1000;
  //plasma lifetime, set to the samplingperiod by default
  double lifetime=.1;

  double txp, tx_on=-999999999., tx_off=999999999.;

  int includeCW_flag=0;//whether to simulate the direct signal as well.
  
  vector<double> amplitudes, timeofarrival, phases, field, plasma; 

  //variables for refraction manipulation
  
  double k_r, c_light_r, mag1, mag2;//, tof, txphase, kx;
  //distance from the antennas to the interface, must be set by user
  vector<double> tx_interface_dist{1};
  vector<double> rx_interface_dist{1};
  //relative index of refraction, calculated to always be >1.
  double n_rel=1.5;

  
  double testvalue, maxval=0, avgval=0, num=0;

  //public:

  int ntx=1;
  int nrx=1;
  vector<HepLorentzVector> tx{1};//transmitters, allow for multiple
  vector<HepLorentzVector> rx{1};//recievers, allow for multiple

  int TX_ITERATOR=0;
  int RX_ITERATOR=0;

  int FILL_BY_EVENT=1;
  int FILL_PARTICLE_INFO=0;
  TString PARTICLE_INFO_FILENAME="";
  int MAKE_SUMMARY_FILE=0;
  
  int NPRIMARIES_SET=0;
  int PRIMARY_ENERGY_SET=0;
  int SCALE_BY_ENERGY=0;


  int REAL_DATA=0;
  RadioScatter();
  
  void makeOutputFile(TString filename);
  void makeOutputTextFile(char* filename);
  void writeToTextFile();
  void makeTimeHist();
  void setMakeSummary(double val);
  int setNTx(double n);
  int setNRx(double n);
  void setTxPos(double xin, double yin, double zin, int index=0);
  void setRxPos(double xin, double yin, double zin, int index=0);
  void setTxPos(Hep3Vector in, int index);
  void setTxPos(Hep3Vector in);
  void setRxPos(Hep3Vector in, int index);
  void setRxPos(Hep3Vector in);
  void setTxFreq(double f);
  void setTxVoltage(double v);
  void setTxPower(double p);
  void setNPrimaries(double n);
  void setAntennaGain(double gain);
  void setPrimaryEnergy(double e);
  void setPrimaryPosition(Hep3Vector p);
  void setPrimaryDirection(Hep3Vector p);
  void setTargetEnergy(double e);//not gonna work yet
  int setScaleByEnergy(double val);
  void setPlasmaLifetime(double l);
  void setPolarization(char * p);
  void setTxVals(double f, double power, double gain);
  void setRxVals(double s, double gain);
  void setSimulationParameters(double n, char* tx_rx_pol, double relative_index_of_refraction, int flag);
  void setRelativeIndexOfRefraction(double iof);
  void setCalculateUsingAttnLength(double val=0.);
  void setRecordWindowLength(double nanoseconds);
  void setRxSampleRate(double rate);
  void setTxInterfaceDistX(double dist, int index=0);
  void setRxInterfaceDistX(double dist, int index=0);
  void setShowCWFlag(double i);
  void setTxOnTime(double on);
  void setTxOffTime(double off);
  void setFillByEvent(double i);
  void setFillParticleInfo(double i);
  void setParticleInfoFilename(char * filename);
  Hep3Vector getTxPos(int index=0);
  Hep3Vector getRxPos(int index=0);
  double getFreq();
  double getTxGain(int index,double angle);
  double getRxGain(int index,double angle);
  double getTxAmplitude(int index,HepLorentzVector point);
  double getRxAmplitude(int index, HepLorentzVector point, Hep3Vector j1,  Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);
  double getRxAmplitude(int txindex, int rxindex, HepLorentzVector point);
  double getTxPhase(double t_0);
  double getRxTime(HepLorentzVector point, Hep3Vector j, Hep3Vector l);
  double getRxTime(int index,HepLorentzVector point);
  double getTxTime(int index,HepLorentzVector point, int direct);
  double getRxPhase(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);

  double getRxPhase(int txindex, int rxindex, HepLorentzVector point);
  double getRxPhaseRel(int index,HepLorentzVector point, double v);
  double getAmplitudeFromAt(double E_0, HepLorentzVector from, HepLorentzVector at);
  double getPhaseFromAt(HepLorentzVector from, HepLorentzVector at);
  double doScreening(TTree *tree, int entry);
  double plotCausalCharges(TTree *tree, int entry);
  int checkTxOn(double time);
  Hep3Vector findRefractionPlane(HepLorentzVector p1, HepLorentzVector p2);
  double findRefractionPointZ(double kx, double kz, double jx);
  double findPathLengthWithRefraction(HepLorentzVector p1, HepLorentzVector p2, double interface_dist_x); 
  TH1F* getDirectSignal(int txindex, int rxindex, const TH1F*);
  double getDirectSignalPhase(int txindex, int rxindex,HepLorentzVector point);
  double makeRays(HepLorentzVector point, double e, double l, double e_i);
  double makeRays(HepLorentzVector point, double e, double l, double e_i, TH1F *hist);
  double power();
  void draw();
  void printEventStats();
  void writeHist(TString filename, float num_events);
  vector<vector<TH1F*>> scaleHist(float num_events);
  int writeRun(float num_events=1., int debug=0);
  int writeEvent(int debug=0);
  int makeSummary(TFile *f);
  void close();


  void testFunc(double testVal);


private:

  int REFRACTION_FLAG=0;
  bool RSCAT_HIST_DECLARE=false;
  bool RSCAT_HIST_RESIZE=false;
  int fRunCounter=0;
};

#endif

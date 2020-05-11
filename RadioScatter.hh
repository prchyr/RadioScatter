/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
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
#include "TAxis.h"

#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "RadioScatterEvent.hh"
#include "RSEventSummary.hh"
#include "IceRayTracing.hh"

 using namespace CLHEP;
 using namespace std;




class RadioScatter{
public:

  RadioScatterEvent event;

  /*
    a word about units. 

    GEANT uses mm and ns as the length and time, but for RF stuff things are best defined in m. so for length calculations as they pertain to RF fields, the lengths are converted into meters. 

    velocity in geant is mm/ns, so for things like phase calculations we employ these native units. 

    we use volts for the fields.
    we use watts for power units.

    if you set a distance unit for radioscatter, always multiply by the unit, like 100*m for 100 meters. for RF units like volts and watts, just use the default RadioScatter (e.g. just provide the number)

  
  
 *************************************************
 constants and variables. function definitions below.

 ***********************************************
  some global things and unit defaults. if there is no default, it must be set.

  */
  TString output_file_name;
  const char* pol = (char*)"horizontal";//default, also set as default in set_simulation_paramaters
  TString polarization="horizontal";
  double x_offset=0*m, z_offset = 2.*m, y_offset = 5.*m;//unused

  double tcs = .655e-24;//thompson cross section
  double collisionalCrossSection = 3.e-13;//mm^-3, from NIST, converted in to mm^-3
  double n_primaries = 1;//set this based on number of events in run
  double tx_voltage = 1.;//V
  double zscale=1.;//the longitudinal scale factor
  double tscale=1.;//the time scale factor
  double impedance = 50;//ohms
  double tx_gain=3.;//transmitter gain default
  double rx_gain=3.;//receiver gain default
  double step_length=1;//default g4 steplength
  double E_i=.000038;//default electron ion pair energy
  double e_charge_cgs = 4.803e-10;//statcoloumbs
  double e_mass_cgs = 9.109e-28;//g
  double k_b = 1.38e-23;//j/k
  double c_light_mns=c_light/m;//c light in m/ns


  //4*pi e^2/m_e
  double plasma_const = 4.*pi*e_charge_cgs*e_charge_cgs/e_mass_cgs;
  //e^2/(4pi epislon0 m c^2), in units of m
  double e_radius=classic_electr_radius/m;
  double nu_col = 0;//collisional frequency
  
  double half_window = 300;//number of nanoseconds in 1/2 of the record window. can be changed;   

  int useAttnLengthFlag=0;//use attenuation length?
  double attnLength=1400.;//average length for upper 1.5km at pole. is changeable

  //misc rf constants
  double frequency=1000, period, lambda, k, omega;
  double phase0=0.;

  //these are changed based on the gain, so the defaults here are meaningless.
  double rxEffectiveHeight=1., txEffectiveHeight=1., txFactor=1., rxFactor=1.;

  //receiver stuff
  double samplerate=10, samplingperiod=.1, start_time=0, end_time=1000;
  //plasma lifetime, set to the samplingperiod by default
  double lifetime=.1;
  //timing variables for pulsed CW
  double txp, tx_on=-999999999., tx_off=999999999.;

  int includeCW_flag=0;//whether to simulate the direct signal as well.
  //misc constants
  vector<double> amplitudes, timeofarrival, phases, field, plasma; 

  //variables for refraction manipulation
  double k_r, c_light_r, mag1, mag2;//, tof, txphase, kx;
  //distance from the antennas to the interface, must be set by user
  vector<double> tx_interface_dist{1};
  vector<double> rx_interface_dist{1};
  //relative index of refraction, calculated to always be >1.
  double n_rel=1.5;

  

  //some histograms 
  TH1F *fft_hist, *power_hist;
  TH1F *t_h = new TH1F("eventHist", "eventHist", 100, 0, 10);
  TH1F *re_h = new TH1F("reHist", "reHist", 100, 0, 10);
  TH1F *im_h = new TH1F("imHist", "imHist", 100, 0, 10);
  TGraph *event_gr = new TGraph();
  
  vector<vector<TH1F*>> time_hist;
  vector<vector<TH1F*>> re_hist;
  vector<vector<TH1F*>> im_hist;
  vector<vector<TGraph*>> event_graph;
  TH1F* testHist0;
  TH1F* testHist1;
  
  //public:

  int ntx=1;
  int nrx=1;
  vector<HepLorentzVector> tx{1};//transmitters, allow for multiple
  vector<HepLorentzVector> rx{1};//recievers, allow for multiple

  //flags
  int TX_ITERATOR=0;
  int RX_ITERATOR=0;

  int FILL_BY_EVENT=1;
  int FILL_PARTICLE_INFO=0;
  TString PARTICLE_INFO_FILENAME="";
  int MAKE_SUMMARY_FILE=0;
  int TX_GAIN_SET=0;
  int RX_GAIN_SET=0;
  int TX_FREQ_SET=0;
  int RX_FREQ_SET=0;
  
  int NPRIMARIES_SET=0;
  int PRIMARY_ENERGY_SET=0;
  int SCALE_BY_ENERGY=0;


  int REAL_DATA=0;
  RadioScatter();

  /*

*******************************************
function definitions

***********************************************

this is a mandatory call that needs to come before the others. it sets the output file and makes a RadioScatterEvent object that is filled with all of the outpus. 
   */
  void makeOutputFile(TString filename);
  //unused
  void makeOutputTextFile(char* filename);
  //unused
  void writeToTextFile();

  /*
    these are the settings that can be set in the simulation. 
    These can be set at the start of any radioscatter simulation. the defaults are sensible (maybe) such that not all need to be set. 
    
    many of these can be set from within GEANT4 using the RSMessenger file in the included simulations. 

    they are described individually below
   */

  //flag saying to make a useful summary file. 
  void setMakeSummary(double val);
  //set the number of transmitters. currently only 1 is allowed, but we've included framework to add more in future. 
  int setNTx(double n);
  //set the number of receivers. this can be any arbitrary number.
  int setNRx(double n);
  //set the transmitter position. the index goes from 0 to ntx-1
  void setTxPos(double xin, double yin, double zin, int index=0);
  //set the receiver position. the index goes from 0 to nrx-1
  void setRxPos(double xin, double yin, double zin, int index=0);
  //set the tx position using an Hep3Vector object
  void setTxPos(Hep3Vector in, int index);
  //same, but with built-in indexing (e.g. each time this is called the next tx will be set up to ntx-1);
  void setTxPos(Hep3Vector in);
    //set the rx position using an Hep3Vector object
  void setRxPos(Hep3Vector in, int index);
    //same, but with built-in indexing (e.g. each time this is called the next rx will be set up to nrx);
  void setRxPos(Hep3Vector in);
  //set the transmitter frequency
  void setTxFreq(double f);
  //set the transmitter voltage. can be superceded by setting power or vice versa (V)
  void setTxVoltage(double v);
  //set the tx power (W) 
  void setTxPower(double p);
  /*
set the number of primaries. this is essentially a scaling factor used to achieve higher-energy showers than GEANT can provide in a reasonable time. for a 10PeV shower, for example, you could simulate a 10GeV primary and then set nPrimaries to 1e7, to get the equivalent density of ta 10PeV shower. to get the longitudinal profile correct, you'd need to setScaleByEnergy = 1 below. 

default is 1.
   */
  void setNPrimaries(double n);
  //set the receiver gain in dB
  void setReceiverGain(double gain);
  //same
  void setRxGain(double gain);
  //set the transmitter gain in dB
  void setTransmitterGain(double gain);
  //same
  void setTxGain(double gain);
  /*
set the energy of the primary. this is important when running over some arbitrary shower input file, as there is no way for radioscatter to know the energy of the primary. it is used for output files and stuff and for scaling factors etc.
   */
  void setPrimaryEnergy(double e);
  //set the position of the primary. useful for several calculations
  void setPrimaryPosition(Hep3Vector p);
  //set the direction of the primary. 
  void setPrimaryDirection(Hep3Vector p);
  //not used
  void setTargetEnergy(double e);
  //not used
  void setCrossSection(double val);
  /*
use this flag to scale the shower by some amount. to use it, you first must call setPrimaryEnergy to the energy of the primary, and then using setNPrimaries, this will calculate the 'target' energy. it they scales the shower accordingly in the longitudinal direction.
   */
  int setScaleByEnergy(double val);
  //set the plasma lifetime in nanoseconds
  void setPlasmaLifetime(double l);
  //set the antenna polarization. very primitive now: vertical=(0,1,0), horizontal=(0,0,1). TODO: fix this.
  void setPolarization(const char * p);
  //not really useful. set them individually instead.
  void setTxVals(double f, double power, double gain);
  //same
  void setRxVals(double s, double gain);
  //same
  void setSimulationParameters(double n, char* tx_rx_pol, double relative_index_of_refraction, int flag);
  //this is n1/n2 for n1>n2. used for refraction calculations when the tx and rx are in different media.
  void setRelativeIndexOfRefraction(double iof);
  //set to calculate the RF propagation with attenuation losses
  void setCalculateUsingAttnLength(double val=0.);
  //set the length of the save window in ns
  void setRecordWindowLength(double nanoseconds);
  //set the sample rate of the receiver in GS/s
  void setRxSampleRate(double rate);
  //only used for refraction calculations when tx and rx are in different media.
  void setTxInterfaceDistX(double dist, int index=0);
  //same
  void setRxInterfaceDistX(double dist, int index=0);
  //set this flag to show the pure CW from transmitter to receiver. default is to have it off.
  void setShowCWFlag(double i);
  //set the time the transmitter is on. useful for pulsed CW, but otherwise don't mess with it and the program assumes constant CW.
  void setTxOnTime(double on);
  //set the time the tx is off. see above
  void setTxOffTime(double off);
  //when running in GEANT4, save individual events (1) or average over all events in a run (0)
  void setFillByEvent(double i);
  //save the tuples in GEANT4 that are filled with information about the particle tracks and steps. slow and will eat big memory so be careful.
  void setFillParticleInfo(double i);
  //set the filename for the above information. 
  void setParticleInfoFilename(char * filename);
  //get the transmitter position
  Hep3Vector getTxPos(int index=0);
  //get the receiver position
  Hep3Vector getRxPos(int index=0);
  //get the transmitter frequency
  double getFreq();
  //get the transmitter gain at a certain angle (not implemented)
  double getTxGain(int index,double angle);
  //get the receiver gain at a certain angle (not implemented)
  double getRxGain(int index,double angle);
  /*
    these fuctions are all that you need to call (for each point you want to scatter radio from.) 

    just give them a 4 vector, the energy deposited in a region, one length of that region, and the ionization energy of the material. 

   */
  
  double makeRays(HepLorentzVector point, double e, double l, double e_i);
  double makeRays(HepLorentzVector point, double e, double l, double e_i, TH1F *hist);
  double makeRaysRayTrace(HepLorentzVector point, double e, double l, double e_i);
  //print out some event statistics. not used much. 
  void printEventStats();
  //scale the histogram (when averaging over several events)
  vector<vector<TH1F*>> scaleHist(float num_events);
  //write a full run.
  int writeRun(float num_events=1., int debug=0);
  //write a single event.
  int writeEvent(int debug=0);
  //make a useful summary file from the full RadioScatterEvent
  int makeSummary(TFile *f);
  //close
  void close();



private:
  void makeTimeHist();
  /*
    These functions are all called from and within makeRays, and so are delicate and shouldn't be messed with/called. hence the private. 
   */
double getTxAmplitude(int index,HepLorentzVector point);
  double getRxAmplitude(int index, HepLorentzVector point, Hep3Vector j1,  Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);
  double getRxAmplitude(int txindex, int rxindex, HepLorentzVector point);
  double getTxPhase(double t_0);
  double getRxTime(HepLorentzVector point, Hep3Vector j, Hep3Vector l);
  double getRxTime(int index,HepLorentzVector point);
  double getTxTime(int index,HepLorentzVector point, int direct);
  double getRxPhase(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2);

  double getRxPhase(int txindex, int rxindex, HepLorentzVector point);
  double getRxPhaseRayTrace(HepLorentzVector point, double tof, double tx_pr, double pr_rx);
  double getRxPhaseRel(int index,HepLorentzVector point, double v);
  double getAmplitudeFromAt(double E_0, HepLorentzVector from, HepLorentzVector at);

  int getRxInfoRayTrace(int txindex,int rxindex, HepLorentzVector point, double *rxPhase, double *rxAmplitude, double *rxTime);
  double getPhaseFromAt(HepLorentzVector from, HepLorentzVector at);
  double doScreening(TTree *tree, int entry);
  double plotCausalCharges(TTree *tree, int entry);
  int checkTxOn(double time);
  Hep3Vector findRefractionPlane(HepLorentzVector p1, HepLorentzVector p2);
  double findRefractionPointZ(double kx, double kz, double jx);
  double findPathLengthWithRefraction(HepLorentzVector p1, HepLorentzVector p2, double interface_dist_x); 
  TH1F* getDirectSignal(int txindex, int rxindex, const TH1F*);
  double getDirectSignalPhase(int txindex, int rxindex,HepLorentzVector point);


  int REFRACTION_FLAG=0;
  bool RSCAT_HIST_DECLARE=false;
  bool RSCAT_HIST_RESIZE=false;

  int fRunCounter=0;
};

#endif

/*
todo:
antenna effective height
antenna effective area
(currently these are both set to 1 accross all frequencies)
 
 
 */
#ifndef R_Scat
#define R_Scat


#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

//#include "TROOT.h"

//#include "TSystem.h"

 using namespace CLHEP;
 using namespace std;
// #ifdef __CLING__
// #pragma link C++ class test;
// #endif


//ClassImp(test);
class RadioScatter{


  typedef struct{
    Hep3Vector direction, position;
    double primary_energy, n_primaries, power, tot_n_scatterers;
    TH1F * event_hist;
    TH1F * re_hist;
    TH1F * im_hist;

  }eventTree;


  eventTree event;

  int REFRACTION_FLAG=0;
  TString output_file_name;
  char* pol = "horizontal";//default, also set as default in set_simulation_paramaters
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
  TH1F *time_hist = new TH1F("time_hist", "time_hist", 100, 0, 10);
  TH1F *re_hist = new TH1F("re_hist", "re_hist", 100, 0, 10);
  TH1F *im_hist = new TH1F("im_hist", "im_hist", 100, 0, 10);
  //e^2/m_e
  double plasma_const = sqrt(4*pi*electron_charge*electron_charge/electron_mass_c2);

  //e^2/(4pi epislon0 m c^2)
  double cross_section=classic_electr_radius;

  
  

  //  ofstream of;
  
  double frequency, period, lambda, k;

  double samplerate, samplingperiod, start_time=0, end_time=1000;

  double txp, tx_on=0., tx_off=10000.;

  int includeCW_flag;//whether to simulate the direct signal as well.
  
  vector<double> amplitudes, timeofarrival, phases, field, plasma; 

  //variables for refraction manipulation
  
    double k_r, c_light_r, mag1, mag2, tof, txphase, kx;
  //distance from the antennas to the interface, must be set by user
  double tx_interface_dist,rx_interface_dist;
  //relative index of refraction, calculated to always be >1.
  double n_rel;

  
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
  void setNprimaries(double n);
  void setPolarization(char * p);
  void setTxVals(double f, double power, double gain);
  void setRxVals(double s, double gain);
  void setSimulationParameters(double n, char* tx_rx_pol, double relative_index_of_refraction, int flag);
  void setRelativeIndexOfRefraction(double iof);
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
  void writeEvent(TString filename, float num_events);
  void close();


  void testFunc(double testVal);
   
};



inline RadioScatter::RadioScatter(){
}

inline void RadioScatter::makeOutputFile(TString filename){
  output_file_name = filename;
  TFile *f = new TFile(filename, "recreate");
  TTree *t = new TTree("tree", "received events");
  t->Branch("event_hist",  &time_hist);
  t->Branch("primary_energy", &event.primary_energy);
  t->Branch("n_primaries", &event.n_primaries);
  t->Branch("position", &event.position);
  t->Branch("direction", &event.direction);
  t->Write();
  f->Close();
}
//just makes the output file stream for a text file
inline void RadioScatter::makeOutputTextFile(char* filename){
  //  ofstream f(filename);
  //  of=f;
}


inline void RadioScatter::makeTimeHist(){
  Hep3Vector dist = tx.vect()+rx.vect();
  //cout<<dist.mag()<<endl;
  start_time = 0.;
  end_time = (dist.mag()/c_light) *3.;
  //  outfile = new TFile("/home/natas/Documents/physics/geant/root/time.root", "RECREATE");
  time_hist->SetBins(samplerate*(end_time-start_time), start_time, end_time);
  re_hist->SetBins(samplerate*(end_time-start_time), start_time, end_time);
  im_hist->SetBins(samplerate*(end_time-start_time), start_time, end_time);
  //outfile.SetOption("RECREATE");
  //time_hist->SetBins(32000, 0, 16000);
  //  time_hist->SetBins(32000, 0, 3200);
  //  time_hist = hist;
  //  cout<<"hist initialized"<<endl;
 }


inline    void RadioScatter::setTxPos(double xin, double yin, double zin){
    tx.setX(xin);
    tx.setY(yin);
    tx.setZ(zin);
    //    makeTimeHist();
  }
inline   void RadioScatter::setRxPos(double xin, double yin, double zin){
    rx.setX(xin);
    rx.setY(yin);
    rx.setZ(zin);
  }
inline    void RadioScatter::setTxPos(Hep3Vector in){
  tx.setX(in.x());
  tx.setY(in.y());
  tx.setZ(in.z());
  }
inline   void RadioScatter::setRxPos(Hep3Vector in){
  rx.setX(in.x());
  rx.setY(in.y());
  rx.setZ(in.z());
  }
inline  void RadioScatter::setTxVals(double f, double power=1., double gain=1.){
    frequency = f*megahertz;
    omega = frequency*twopi;
    period = 1./omega;
    lambda = c_light/frequency;
    k = omega/c_light;
    tx_gain = gain;
    tx_voltage = power;
  }
inline void RadioScatter::setTxFreq(double f){
  frequency = f*megahertz;
  omega = frequency*twopi;
  period = 1./omega;
  lambda = c_light/frequency;
  k = omega/c_light;
  cout<<"tx frequency: "<<f<<endl;
}
inline void RadioScatter::setTxVoltage(double v){
  tx_voltage = v;
  event.power=v;
}

inline void RadioScatter::setTxOnTime(double on){
  tx_on = on;
}
inline void RadioScatter::setTxOffTime(double off){
  tx_off = off;
}
inline void RadioScatter::setNprimaries(double n){
  n_primaries=n;
  event.n_primaries=n_primaries;
  cout<<"n primaries: "<<n_primaries<<endl;
}
inline void RadioScatter::setPolarization(char * p){
  pol=p;
  cout<<"polarization: "<<p<<endl;
}
inline  void RadioScatter::setRxVals(double s=1., double gain=1.){
    samplerate = s*nanosecond;
    samplingperiod = 1./samplerate;
    rx_gain = gain;
    //  n_primaries=n;
  }

inline void RadioScatter::setShowCWFlag(double i){
  includeCW_flag=(int)i;
}

inline  void RadioScatter::setSimulationParameters(double n=1., char* tx_rx_pol="horizontal", double relative_index_of_refraction=1.5, int flag=0){
  n_primaries=n;
  pol=tx_rx_pol;
  n_rel = relative_index_of_refraction;
  k_r = omega/(c_light/n_rel);
  c_light_r = c_light/n_rel;
  includeCW_flag=flag;
}
inline void RadioScatter::setRelativeIndexOfRefraction(double iof){
  n_rel=iof;
  k_r = omega/(c_light/n_rel);
  c_light_r = c_light/n_rel;
}
inline void RadioScatter::setRxSampleRate(double rate){
  samplerate=rate*nanosecond;
  samplingperiod = 1./samplerate;
}
inline void RadioScatter::setTxInterfaceDistX(double dist=0){
    tx_interface_dist = abs(dist);
    REFRACTION_FLAG=1;
}
inline void RadioScatter::setRxInterfaceDistX(double dist=0){
  rx_interface_dist = abs(dist);
}
inline Hep3Vector RadioScatter::getTxPos(){
  return tx.vect();
}
inline Hep3Vector RadioScatter::getRxPos(){
  return rx.vect();
}
inline double RadioScatter::getFreq(){
  return frequency;
}

  //can input gain pattern later
inline  double RadioScatter::getTxGain(double angle=0.){
    double gain = tx_gain;
    return gain;
  }
inline  double RadioScatter::getRxGain(double angle=0.){
    double gain = rx_gain;
    return gain;
  }
inline  double RadioScatter::getTxAmplitude(HepLorentzVector point){
    double gain = getTxGain(point.vect().theta());
    double power = tx_voltage*gain;
    return power;
  }


/*
use the calculated refraction vectors (from makeRays()) to sort out the correct amplitude.
 */

inline double RadioScatter::getRxAmplitude(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2){
  double dist = j1.mag()+j2.mag()+l1.mag()+l2.mag();
  //refraction things:


  double alpha1 = atan(j1.z()/j1.x());
  double alpha_prime1 = atan(l1.z()/l1.x());
  double alpha2 = atan(j2.z()/j2.x());
  double alpha_prime2 = atan(l2.z()/l2.x());

  //E for polarization parallel to plane of incidence
  double E1_para = (2.*cos(alpha1))/(cos(alpha1)+(n_rel*cos(alpha_prime1)));
  double E1_perp = (2.*n_rel*cos(alpha1))/(n_rel*n_rel*cos(alpha1)+(n_rel*cos(alpha_prime1)));

  //E for polarization perpendicular to plane of incidence
  double E2_para = (2.*cos(alpha2))/(cos(alpha2)+(n_rel*cos(alpha_prime2)));
  double E2_perp = (2.*n_rel*cos(alpha2))/(n_rel*n_rel*cos(alpha2)+(n_rel*cos(alpha_prime2)));

  //find angle between plane of incidence and polarization vector
  double theta = atan(point.y()/point.z());

  if(pol=="vertical")theta = theta+(pi/2.);

  double amp1 = sqrt(pow(E1_para*cos(theta), 2)+pow(E1_perp*sin(theta), 2));
  double amp2 = sqrt(pow(E2_para*cos(theta), 2)+pow(E2_perp*sin(theta), 2));

  double angle_dependence = (l1.unit()).dot(l2.unit());
  double amplitude = (tx_voltage/dist)*amp1*amp2*angle_dependence;
  return amplitude;
}

//non-refracted amplitude
inline double RadioScatter::getRxAmplitude(HepLorentzVector point){
  double dist = (tx.vect()-point.vect()).mag()+(rx.vect()-point.vect()).mag();
  //refraction things:
  Hep3Vector one=tx.vect()-point.vect();
  Hep3Vector two=point.vect()-rx.vect();

  //  double amplitude = (tx_voltage/dist)*one.unit().dot(two.unit());
  double amplitude = (tx_voltage/dist);

  return amplitude;
}


inline   double RadioScatter::getTxPhase(double t_0){
  //    HepLorentzVector tx_pr = tx-point;
  double phase = omega*t_0;
  return phase;
}

inline  double RadioScatter::getRxTime(HepLorentzVector point, Hep3Vector l, Hep3Vector j){
  double dist = j.mag()+l.mag();
  double time = point.t()+(dist/c_light);
    return time;
  }

inline  double RadioScatter::getRxTime(HepLorentzVector point){
  //  double dist = findPathLengthWithRefraction(rx,point, rx_interface_dist);
 
  double dist = (rx.vect()-point.vect()).mag();
  double time = point.t()+(dist/c_light);

  return time;
  }

inline  double RadioScatter::getTxTime(HepLorentzVector point, int direct=0){
  //    Hep3Vector dist = point.vect()-tx.vect();
  double dist=0;
  if(direct==0){
    if(REFRACTION_FLAG==1){
      dist = findPathLengthWithRefraction(tx,point, tx_interface_dist);
    }
    else{
      dist= (point.vect()-tx.vect()).mag();
    } 
  }
  else{
    dist = (point.vect()-tx.vect()).mag();
  }
    double time = point.t()-(dist/c_light);

    return time;
  }

inline int RadioScatter::checkTxOn(double time){
  if(time>tx_on&&time<tx_off){
    return 1;
  }
  else{
    return 0;
  }
}

inline double RadioScatter::getDirectSignalPhase(HepLorentzVector point){
  double dist = (tx-rx).vect().mag();
  double t_0 = point.t()-(dist/c_light);
  return t_0*omega;
}
/*
using the refraction vectors (calculated in makeRays), calculate the wave vectors (including with modified k value in the medium) and the 
total time-of-flight and phase. the values n_rel, k_r, and c_light_r are the relative index of refraction, and 
the modified wavenumber and speed of light (for the medium) respectively  
*/
inline double RadioScatter::getRxPhase(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2){
  
  //debug:check that snell's law is satisfied for the found paths.
  // cout<<(j1.z()/j1.mag())/(l1.z()/l1.mag())<<endl;
  // cout<<(j2.z()/j2.mag())/(l2.z()/l2.mag())<<endl<<endl;
  
  //calculate the time of flight using correct values for velocity
  double  tof = j1.mag()/c_light + l1.mag()/c_light_r + j2.mag()/c_light + l2.mag()/c_light_r;

  double txphase = getTxPhase(point.t()-(j1.mag()/c_light + l1.mag()/c_light_r));//find phase at retarded time
  //wave vector calculation, with correct phase velocity
  Hep3Vector  kvec1 = k*j1;
  Hep3Vector  kvec2 = k_r*l1;
  Hep3Vector  kvec3 = k_r*l2;
  Hep3Vector  kvec4 = k*j2;
  Hep3Vector  ktot = kvec1+kvec2+kvec3+kvec4;
  double  kx = kvec1.mag()+kvec2.mag()+kvec3.mag()+kvec4.mag();

  return ((kx) - omega*tof + txphase);
  }

//non-refraction phase finder
inline   double RadioScatter::getRxPhase(HepLorentzVector point){
    double rxtime = getRxTime(point);//find advanced time
    double txtime = getTxTime(point);//find retarted time
    double txphase = getTxPhase(txtime);//find phase at retarded time
    double tof = abs(rxtime-txtime);//time of flight
    HepLorentzVector tx_pr=tx-point, pr_rx = point-rx;//make vectors
    //wave number addition
    Hep3Vector kvec1 = k*tx_pr.vect();
    Hep3Vector kvec2 = k*pr_rx.vect();
    Hep3Vector ktot = kvec1+kvec2;
    double kx = ktot.mag();
    //calculate compton effects
    double inv_omega_c = (1/omega)+(1/omega_e)*(1-cos(tx_pr.vect().unit().angle(pr_rx.vect().unit())));
    omega_c = 1/inv_omega_c;
    //    cout<<txtime<<" "<<txphase<<" "<<rxtime<<endl;
    return ((kx) - omega*tof + txphase);
  }



  //calculate relativistic doppler shift
// inline  double RadioScatter::getRxPhaseRel(HepLorentzVector point, double v){
//     double time = getRxTime(point);
//     HepLorentzVector tx_pr=tx-point, pr_rx = point-rx;
//     double gamma = 1./(sqrt(1.-pow(v, 2)));
//     double omega_prime = gamma*omega*(1.+v*cos(tx_pr.vect().unit().angle(pr_rx.vect().unit())));
//     double k_rel = omega_prime/c_light;
//     Hep3Vector kvec1 = k_rel*tx_pr.vect();
//     Hep3Vector kvec2 = k_rel*pr_rx.vect();
//     Hep3Vector ktot = kvec1+kvec2;
//     double kx = ktot.mag();
  
//     return ((kx) - (omega_prime*time));
//   }


/*
simulate the direct signal that would be seen at the receiver for CW
 */
inline TH1F * RadioScatter::getDirectSignal(const TH1F *in){
  TH1F *outhist=(TH1F*)in;
  HepLorentzVector point=rx;
  double rx_amp, rx_ph, amp;
  rx_amp = tx_voltage*m/(tx.vect()-rx.vect()).mag();
  //cout<<"rx amp: "<<tx_voltage<<" "<<tx_voltage*m<<" "<<(tx.vect()-rx.vect()).mag()<<" "<<rx_amp;
  int size = in->GetNbinsX();
  //cout<<size<<endl; 
  for(int i=0;i<size;i++){
    point.setT((double)i*samplingperiod);

    rx_ph = getDirectSignalPhase(point);
    amp =rx_amp*cos(rx_ph);
    //    cout<<amp<<endl;
    //outhist->Fill(i, amp);
    //cout<<"asdfasfd"<<endl;
    rx.setT(i*samplingperiod);
    if(checkTxOn(getTxTime(rx, 1))==1){
      outhist->Fill(i*samplingperiod, (rx_amp*cos(rx_ph)));
    }
  }

  return outhist;
}

/*
the main function. q1 and q2 are the direct path vectors between tx-ip and rx-ip respectively. 

j1, j2 are the vectors from tx-refraction point and rx-refraction point respectively, for their corresponding interfaces.
l1, l2 are the vectors from these refraction points to the interaction point.

calculate the phase, the amplitude, and the prefactors for cross-section, 

 */

inline  double RadioScatter::makeRays(HepLorentzVector point, double e, double l, double e_i){
#ifndef RSCAT_HIST_DECLARE
  makeTimeHist();
#endif
#define RSCAT_HIST_DECLARE

  double rx_time, rx_amplitude, rx_phase;
  // int dumb=1;
  // if(dumb==1){
  if(checkTxOn(getTxTime(point))==1){
    
    if(REFRACTION_FLAG==1){
      Hep3Vector  q1 = findRefractionPlane(tx, point);//make a plane where the refraction will happen
      Hep3Vector j1;
      j1.setZ(findRefractionPointZ(q1.x(), q1.z(), tx_interface_dist));//find the refraction point in this plane on interface 
      Hep3Vector  q2 = findRefractionPlane(rx, point);
      Hep3Vector  j2;
      j2.setZ(findRefractionPointZ(q2.x(), q2.z(), rx_interface_dist));
      j1.setX(tx_interface_dist);
      j1.setY(0);
      j2.setX(rx_interface_dist);
      j2.setY(0);
      Hep3Vector l1;
      l1.set(q1.x()-tx_interface_dist, 0., q1.z()-j1.z());
      Hep3Vector l2;
      l2.set(q2.x()-rx_interface_dist, 0., q2.z()-j2.z());
      
      
      
      //get the reflected signal amplitude and phase
      rx_phase = getRxPhase(point, j1, j2, l1, l2);
      rx_amplitude = getRxAmplitude(point, j1, j2, l1, l2);
      
      
      
      //get the time ray would arrive at rx, to fill histogram
      rx_time = getRxTime(point, j2, l2);
    }
    else{
      rx_phase=getRxPhase(point);
      rx_amplitude=getRxAmplitude(point);
      rx_time=getRxTime(point);
    }
      
    step_length=l;//to make our density approximation
    double n = e/e_i;//edeposited/ionization E
    double n_e =1;
    //calculate plasma freq and collison freq
    if(step_length!=0){
      //electron number density, using step length cube
      n_e = n*n_primaries/pow(step_length, 3);
      //collision frequency (approximation from Cravens), multiplied by 3 for e/i, e/e, e/n (BAD APPROXIMATION FIX PLZ)
      nu_col = 3.*54.*n_e/pow(e_i/k_Boltzmann, 1.5);
      //plasma frequency, TODO
      //    cout<<nu_col<<endl;
      omega_0=plasma_const*sqrt(n_e);
    }
    event.tot_n_scatterers+=n;//track total number of scatterers
    //for each ionization e scatterer
    //  double filter = exp(-pow(omega, 2)/pow(omega_0, 2));
    double filter=1.;//TODO
    //the full scattering amplitude pre-factor  
    double prefactor = filter*n*n_primaries*cross_section*rx_amplitude*omega/(pow(omega, 2)+pow(nu_col, 2));
    //now calculate real and imaginary e fields
    //double E= prefactor*omega*cos(rx_phase)-prefactor*nu_col*sin(rx_phase);
    double E_real= prefactor*omega*cos(rx_phase)-prefactor*nu_col*sin(rx_phase);
      
    double E_imag = prefactor*-nu_col*cos(rx_phase)+prefactor*omega*sin(rx_phase);
      
      
    //make sure phase term is analytic
    if(rx_phase>0&&rx_phase<999&&rx_amplitude>0&&rx_amplitude<999){
      time_hist->Fill(rx_time, E_real);
      re_hist->Fill(rx_time, E_real);
      im_hist->Fill(rx_time, E_imag);
    }
      
    return E_real;
  }
  else{
    return 0;
  }
}


//calculate the point of refraction in the z direction. n_rel is the relative index of refraction n2/n1
inline double RadioScatter::findRefractionPointZ(double kx, double kz, double jx){
  double  v =n_rel;
  double jz=kz/2. + sqrt(pow(kz,2) - (2*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) + (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))))/2. - sqrt(2*pow(kz,2) - (4*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) - (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) - pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))) + (8*pow(kz,3) + (16*pow(jx,2)*kz*pow(v,2))/(-1 + pow(v,2)) - (8*kz*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(-1 + pow(v,2)))/(4.*sqrt(pow(kz,2) - (2*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) + (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))))))/2;
return jz;
}

/*
analytic solution above requires that the tx and rx are in a plane. this function assumes 
this plane is in the plane containing the x-axis, interaction point, and either tx/rx. 

function rotates the k vector to be in the kz plane and then sends this plane to the above function.
might be completely wrong. but resultant vectors satisfy snell's law.
 */

inline Hep3Vector RadioScatter::findRefractionPlane(HepLorentzVector p1, HepLorentzVector p2){
  Hep3Vector k;
  k.setX(abs(p2.x()-p1.x()));

  k.setY(abs(p2.y()-p1.y()));
  k.setZ(abs(p2.z()-p1.z()));
  // cout<<endl<<k.x()<<" "<<k.y()<<" "<<k.z()<<endl;
  // cout<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<endl;
  // cout<<tx.x()<<" "<<tx.y()<<" "<<tx.z()<<endl;
  // cout<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<endl<<endl;
  double theta = atan(k.y()/k.z());
  k.rotateX(theta);
  return k;
}
/*
find the path length including refraction
 */
inline   double RadioScatter::findPathLengthWithRefraction(HepLorentzVector p1, HepLorentzVector p2, double interface_dist_x){
  Hep3Vector k = findRefractionPlane(p1,p2);
  double jz = findRefractionPointZ(k.x(), k.z(), interface_dist_x);
  Hep3Vector j(interface_dist_x, 0., jz);//from source point to refraction point
  Hep3Vector l(k.x()-interface_dist_x, 0., k.z()-jz);//from refraction point to interaction point
  double mag = j.mag()+l.mag();
  //  cout<<mag-k.mag()<<endl;
  return mag;
  
}

inline   double RadioScatter::power(){
  double p;
  //  cout<<time_hist->GetNbinsX()<<endl;
  //time_hist.Draw();
  for(int i=0;i<time_hist->GetNbinsX();i++){
    p+=pow(time_hist->GetBinContent(i), 2);
  }
  //  cout<<p<<endl;
  return p;
  }
//broken
inline  void RadioScatter::draw(){
    time_hist->Draw("l");
  }
inline  void RadioScatter::writeHist(TString filename, float num_events=1.){
  TFile *t = new TFile(filename, "recreate");
  time_hist->Scale(1./num_events);
  re_hist->Scale(1./num_events);
  im_hist->Scale(1./num_events);
  
  TH1F *outhist=0;
  if(includeCW_flag==1){
    //    cout<<"before"<<endl;
    outhist=getDirectSignal((const TH1F*)time_hist);
    //cout<<"after"<<endl;
      //      cout <<getDirectSignal();
    }
  else{
    outhist=time_hist;
  }

  outhist->Write();
  // time_hist.BufferEmpty();
  t->Close();
  // time_hist->Clear();
  //outfile.Write();
  //outfile.Close();
  //return outhist;
  cout<<filename<<endl;
}

inline  TH1F* RadioScatter::scaleHist(float num_events=1.){
  time_hist->Scale(1./num_events);
  re_hist->Scale(1./num_events);
  im_hist->Scale(1./num_events);
  TH1F *outhist=0;
  if(includeCW_flag==1){
    //    cout<<"before"<<endl;
    outhist=getDirectSignal((const TH1F*)time_hist);
    //cout<<"after"<<endl;
      //      cout <<getDirectSignal();
    }
  else{
    outhist=time_hist;
  }
  return outhist;
}

inline void RadioScatter::writeToTextFile(){
  //  of<<
}

inline void RadioScatter::printEventStats(){
  //  event.tot_n_scatterers = event.tot_n_scatterers*n_primaries;
  cout<<setprecision(10);
  cout<<"total number of scatterers: "<<  event.tot_n_scatterers<<endl;
}
inline void RadioScatter::writeEvent(TString filename, float num_events=1.){
  TFile *f = new TFile(filename, "recreate");


  //event.primary_energy = e;
  event.n_primaries = n_primaries;
  //event.position = pos;
  // event.direction = dir;
  event.event_hist = scaleHist(num_events);
  event.re_hist = re_hist;
  event.im_hist = im_hist;
  event.tot_n_scatterers = event.tot_n_scatterers*n_primaries;

  // TTree *t = (TTree*)f->Get("tree");
  // t->SetBranchAddress("event_hist", &time_hist);
  // t->SetBranchAddress("primary_energy", &event.primary_energy);
  // t->SetBranchAddress("n_primaries", &event.n_primaries);
  // t->SetBranchAddress("position", &event.position);  
  TTree *t = new TTree("tree", "treee");
  //  t->Branch("event", &event);
  t->Branch("event_hist", &event.event_hist);
  t->Branch("primary_energy", &event.primary_energy);
  t->Branch("n_primaries", &event.n_primaries);
  //cout<<"here"<<endl;
  f->Append(event.event_hist);
  f->Append(event.re_hist);
  f->Append(event.im_hist);
  f->cd();
  t->Fill();
  cout<<filename<<endl;
  f->Write();
  time_hist->Reset();
  f->Close();
}

inline  void RadioScatter::close(){
    delete(time_hist);
    delete(fft_hist);
    delete(power_hist);
  }

inline void RadioScatter::testFunc(double testVal){
  testvalue = testVal;
}

#endif

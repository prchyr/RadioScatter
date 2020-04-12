/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */

#include "RadioScatter.hh"
 using namespace CLHEP;
 using namespace std;


RadioScatter::RadioScatter(){
  setNTx(1);
  setNRx(1);
}


//mandatory call.
 void RadioScatter::makeOutputFile(TString filename){
  output_file_name = filename;
  TFile *f = new TFile(filename, "recreate");
  TTree *t = new TTree("tree", "received events");
  t->Branch("event", &event);
  f->Write();

  fRunCounter=0;
}
//just makes the output file stream for a text file
//not using now cause txt files are dumb
 void RadioScatter::makeOutputTextFile(char* filename){
  //  ofstream f(filename);
  //  of=f;
}

//called automatically, makes the output hist and figures out where the return pulse will arrive in time
 void RadioScatter::makeTimeHist(){
   //build the transmitter/receiver histograms
   if(RSCAT_HIST_DECLARE==false){
     //these histograms are 'test' histograms which can be declared to do whatever you want. kinda not fully implemented yet though. 
     auto testHistN=2;
     for(int i=0;i<testHistN;i++){
       TString ii=TString::Itoa(i, 10);
       TH1F * temphist=new TH1F("test"+ii, "test"+ii, 100, -10, 20);
       event.testHist.push_back(temphist);
     }
     
     for(int i=0;i<ntx;i++){
       time_hist.push_back(vector<TH1F*>());
       re_hist.push_back(vector<TH1F*>());
       im_hist.push_back(vector<TH1F*>());
       event_graph.push_back(vector<TGraph*>());
       for(int j=0;j<nrx;j++){
	 TString ii=TString::Itoa(i, 10);
	 TString jj=TString::Itoa(j, 10);
	 TH1F * newthist = new TH1F("t"+ii+jj, "t"+ii+jj,10,0, 10);
	 TH1F * newrhist = new TH1F("r"+ii+jj, "r"+ii+jj,10,0, 10);
	 TH1F * newihist = new TH1F("i"+ii+jj, "i"+ii+jj,10,0, 10);
	 TGraph *tgr=0;
	 time_hist[i].push_back(newthist);
	 re_hist[i].push_back(newrhist);
	 im_hist[i].push_back(newihist);
	 event_graph[i].push_back(tgr);
       }
     }
     cout<<"histograms and graphs initialized"<<endl;

   }

   for(int i=0;i<ntx;i++){
     for(int j=0;j<nrx;j++){
       Hep3Vector dv = rx[j].vect()-event.position;
       double dist = dv.mag();
       double time =(dist/c_light)-half_window; 
       
       time<0?start_time=0:start_time=time;
       end_time = start_time+(2*half_window);
       time_hist[i][j]->Reset();
       re_hist[i][j]->Reset();
       im_hist[i][j]->Reset();
       time_hist[i][j]->SetBins(samplerate*(end_time-start_time), start_time, end_time);
       re_hist[i][j]->SetBins(samplerate*(end_time-start_time), start_time, end_time);
       im_hist[i][j]->SetBins(samplerate*(end_time-start_time), start_time, end_time);
     }
   }
 }

void RadioScatter::setMakeSummary(double val){
  MAKE_SUMMARY_FILE=(int)val;
}
void RadioScatter::setRecordWindowLength(double nanoseconds){
  half_window= nanoseconds/2;
}

void RadioScatter::setCalculateUsingAttnLength(double val){
  useAttnLengthFlag=(int)val;
  if((int)val==1){
    cout<<"calculating using attenuation length. (Barwick et al)"<<endl;
  }
}

int RadioScatter::setNTx(double n){
  ntx=(int)n;
  HepLorentzVector vec;
  if(tx.size()==ntx)return 1;
  for(int i=0;i<ntx;i++){
    tx.push_back(vec);
    tx_interface_dist.push_back(0);
  }
  return 1;
}
int RadioScatter::setNRx(double n){
  nrx=(int)n;
  HepLorentzVector vec;
  if(rx.size()==nrx)return 1;
  for(int i=0;i<nrx;i++){
    rx.push_back(vec);
    rx_interface_dist.push_back(0);
  }
  return 1;
}
void RadioScatter::setTxPos(double xin, double yin, double zin, int index){
    tx[index].setX(xin);
    tx[index].setY(yin);
    tx[index].setZ(zin);
    RSCAT_HIST_RESIZE=false;//resize the output histogram
  }
void RadioScatter::setRxPos(double xin, double yin, double zin, int index){
    rx[index].setX(xin);
    rx[index].setY(yin);
    rx[index].setZ(zin);
    RSCAT_HIST_RESIZE=false;
  }
void RadioScatter::setTxPos(Hep3Vector in){
  setTxPos(in, TX_ITERATOR);
  TX_ITERATOR+=1;
}
void RadioScatter::setTxPos(Hep3Vector in, int index){
  tx[index].setX(in.x());
  tx[index].setY(in.y());
  tx[index].setZ(in.z());
  RSCAT_HIST_RESIZE=false;
  }
void RadioScatter::setRxPos(Hep3Vector in){
  setRxPos(in, RX_ITERATOR);
  RX_ITERATOR+=1;
}
void RadioScatter::setRxPos(Hep3Vector in, int index){
  rx[index].setX(in.x());
  rx[index].setY(in.y());
  rx[index].setZ(in.z());
  RSCAT_HIST_RESIZE=false;
  }
  void RadioScatter::setTxVals(double f, double power=1., double gain=1.){
    frequency = f*megahertz;
    omega = frequency*twopi;
    period = 1./omega;
    lambda = c_light/frequency;
    k = omega/c_light;
    tx_gain = gain;
    tx_voltage = sqrt(power*50.);
  }
 void RadioScatter::setTxFreq(double f){
   TX_FREQ_SET=1;
   frequency = f*megahertz;
  omega = frequency*twopi;
  period = 1./omega;
  lambda = c_light/f;
  k = omega/c_light;
  
  if(TX_GAIN_SET!=1){
    //from RICE but with right placement of gain
    txEffectiveHeight=sqrt(lambda*lambda/(tx_gain*480.*pi*pi));
    //from wikipedia antenna factor but with right placement of gain
    txEffectiveHeight=lambda/(9.73*sqrt(tx_gain));
    //first principles losses due to spherical spreading (Friis)         
    rxFactor=sqrt(rx_gain*lambda/(4.*pi));
    txFactor=sqrt(tx_gain*lambda/(4.*pi));  

  }

  attnLength = 1450000;//barwick et all south pole measurement
  cout<<"tx frequency: "<<f<<endl;
  cout<<"omega: "<<omega<<endl;
  cout<<"lambda: "<<lambda<<endl;

  cout<<"attn length (mm): "<<attnLength<<endl;
}
 void RadioScatter::setTxVoltage(double v){
  tx_voltage = v;
  event.txVoltage=v;
  event.txPowerW=pow(v, 2)/50.;
}
void RadioScatter::setTxPower(double p){//W
  tx_voltage = sqrt(p*50.);
  event.txVoltage=tx_voltage;
  event.txPowerW=p;
}
void RadioScatter::setRxGain(double gain){
  setReceiverGain(gain);
}
void RadioScatter::setReceiverGain(double gain){
  double gg=pow(10., gain/10.);
  rx_gain=gg;
  if(RX_GAIN_SET!=1){
    RX_GAIN_SET==1;

    cout<<" receiver gain: "<<gain<<" dB, ("<<gg<<" linear)"<<endl;
    //from rice paper
    rxEffectiveHeight=sqrt(lambda*lambda*rx_gain/(480.*pi*pi));
    //from wikipedia antenna factor
    rxEffectiveHeight=lambda*sqrt(rx_gain)/9.73;
    //first principles losses due to spherical spreading (Friis)         
    rxFactor=sqrt(rx_gain*lambda/(4.*pi));
    
    cout<<"effective height: "<<rxEffectiveHeight<<endl;
  }
}
void RadioScatter::setTxGain(double gain){
  setTransmitterGain(gain);
}
void RadioScatter::setTransmitterGain(double gain){
  if(TX_FREQ_SET==0){
    cout<<"WARNING:antenna gain set before transmit frequency, effective height will be incorrect. please set tx freq before setting antenna gain"<<endl;
  }
  double gg=pow(10., gain/10.);
  tx_gain=gg;
  if(TX_GAIN_SET!=1){
    TX_GAIN_SET==1;

    cout<<"transmitter voltage: "<<tx_voltage<<" transmitter gain: "<<gain<<" dB, ("<<gg<<" linear)"<<endl;

    //from RICE but with the right placement of the gain
    txEffectiveHeight=sqrt(lambda*lambda/(tx_gain*480.*pi*pi));
    //from wikipedia antenna factor but with the right placement of gain
    txEffectiveHeight=lambda/(9.73*sqrt(tx_gain));
    //first principles losses due to spherical spreading (Friis)         
    txFactor=sqrt(tx_gain*lambda/(4.*pi));
    cout<<"effective height: "<<txEffectiveHeight<<endl;
  }
}

 void RadioScatter::setTxOnTime(double on){
  tx_on = on;
}
 void RadioScatter::setTxOffTime(double off){
  tx_off = off;
}
 void RadioScatter::setNPrimaries(double n){
  n_primaries=n;
  event.nPrimaries=n_primaries;

  cout<<"n primaries: "<<n_primaries<<endl;
  NPRIMARIES_SET=1;

}
void RadioScatter::setWeight(double val){
  event.weight=val;
}
void RadioScatter::setCrossSection(double val){
  collisionalCrossSection=val;
}

void RadioScatter::setTargetEnergy(double e){
  setNPrimaries(e-event.primaryEnergy);
  NPRIMARIES_SET=1;
}
int RadioScatter::setScaleByEnergy(double val){
  if((int)val==1){
    if(PRIMARY_ENERGY_SET==0){
      cout<<"primary energy has not been set in RadioScatter, scaling will not work. please call setPrimaryEnergy to tell RadioScatter what energy your primary is!"<<endl;
      return 0;
    }
    if(NPRIMARIES_SET==0){
      cout<<"number of primaries/target energy has not been set in RadioScatter, scaling will not work. please call setNPrimaries or setTargetEnergy to tell RadioScatter how to scale!"<<endl;
      return 0;
    }
    double current_energy = event.primaryEnergy/1000.;//gev
    double target_energy = (current_energy+event.nPrimaries);//gev

    zscale = (3000.*log10(event.nPrimaries)+6000.)/((log10(event.primaryEnergy/1000.)*3000.)+6000);

    tscale = (10.*log10(event.nPrimaries)+22.)/((log10(event.primaryEnergy/1000.)*10.)+22);
    cout<<"scaling activated. zscale="<<zscale<<" , tscale="<<tscale<<endl;

    SCALE_BY_ENERGY=1;
    return 1;
  }
}
void RadioScatter::setPrimaryEnergy(double e){
  event.primaryEnergy=e;
  PRIMARY_ENERGY_SET=1;
}
void RadioScatter::setPrimaryDirection(Hep3Vector p){
  event.direction=p;
  RSCAT_HIST_RESIZE=false;
}
void RadioScatter::setPrimaryPosition(Hep3Vector p){
  event.position=p;
  RSCAT_HIST_RESIZE=false;
}
void RadioScatter::setPolarization(const char * p){
  pol=p;
  cout<<"polarization: "<<pol<<endl;
  polarization=TString(pol);
}

void RadioScatter::setPlasmaLifetime(double l){
  lifetime=l;
}
void RadioScatter::setRxVals(double s=1., double gain=1.){
    samplerate = s*nanosecond;
    samplingperiod = 1./samplerate;
    rx_gain = gain;
  }

 void RadioScatter::setShowCWFlag(double i){
  includeCW_flag=(int)i;
}

void RadioScatter::setSimulationParameters(double n=1., char* tx_rx_pol=(char*)"horizontal", double relative_index_of_refraction=1.5, int flag=0){
  n_primaries=n;
  pol=tx_rx_pol;
  n_rel = relative_index_of_refraction;
  k_r = omega/(c_light/n_rel);
  c_light_r = c_light/n_rel;
  includeCW_flag=flag;
}
 void RadioScatter::setRelativeIndexOfRefraction(double iof){
  n_rel=iof;
  k_r = omega/(c_light/n_rel);
  c_light_r = c_light/n_rel;
}
 void RadioScatter::setRxSampleRate(double rate){
  samplerate=rate*nanosecond;
  samplingperiod = 1./samplerate;
  event.sampleRate=samplerate;
}
void RadioScatter::setTxInterfaceDistX(double dist, int index){
    tx_interface_dist[index] = abs(dist);
    REFRACTION_FLAG=1;
}
void RadioScatter::setRxInterfaceDistX(double dist, int index){
  rx_interface_dist[index] = abs(dist);
 }
void RadioScatter::setFillByEvent(double i){
  FILL_BY_EVENT=(int)i;
}
void RadioScatter::setFillParticleInfo(double i){
  FILL_PARTICLE_INFO=(int)i;
}
void RadioScatter::setParticleInfoFilename(char * filename){
  PARTICLE_INFO_FILENAME=TString(filename);
}
Hep3Vector RadioScatter::getTxPos(int index){
  return tx[index].vect();
}
 Hep3Vector RadioScatter::getRxPos(int index){
  return rx[index].vect();
}
 double RadioScatter::getFreq(){
  return frequency;
}

  //can input gain pattern later
  double RadioScatter::getTxGain(int index,double angle=0.){
    double gain = tx_gain;
    return gain;
  }
  double RadioScatter::getRxGain(int index,double angle=0.){
    double gain = rx_gain;
    return gain;
  }
  double RadioScatter::getTxAmplitude(int index,HepLorentzVector point){
    double gain = getTxGain(point.vect().theta());
    double amplitude = tx_voltage/txEffectiveHeight;
    return amplitude;
  }


/*
use the calculated refraction vectors (from makeRays()) to sort out the correct amplitude.
 */

 double RadioScatter::getRxAmplitude(int index,HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2){
   double dist = ((j1.mag()+j2.mag())/m)*((l1.mag()+l2.mag())/m);
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
  double angle_dependence=1.;

  
  Hep3Vector nhat((point-rx[index]).vect().unit());

  //Hep3Vector vert(0,1.,0), horiz(0,0,1.);
  Hep3Vector vert(1.,1.,1.), horiz(1.,1.,1.); 
  //  l1plane.setTheta(0);
  //  l1plane.setPhi(0);
  
  if(polarization=="vertical"){
    //refraction angle change
    theta = theta+(pi/2.);

    //angle_dependence = vert.cross(nhat).mag();
    angle_dependence = nhat.cross(nhat.cross(vert)).mag();

  }
  else{

    //angle_dependence = horiz.cross(nhat).mag();
    angle_dependence = nhat.cross(nhat.cross(horiz)).mag();

  }
  double amp1 = sqrt(pow(E1_para*cos(theta), 2)+pow(E1_perp*sin(theta), 2));
  double amp2 = sqrt(pow(E2_para*cos(theta), 2)+pow(E2_perp*sin(theta), 2));


  //double amplitude = ((tx_voltage/txEffectiveHeight)/dist)*amp1*amp2*angle_dependence;
  double amplitude = ((tx_voltage*txFactor)/dist)*amp1*amp2*angle_dependence;

  if(useAttnLengthFlag==1){
    double attn_dist = (j1.mag()+j2.mag())+(l1.mag()+l2.mag());
    amplitude=amplitude*exp(-attn_dist/attnLength);
  }
  return amplitude;
}

//non-refracted amplitude
double RadioScatter::getRxAmplitude(int txindex,int rxindex, HepLorentzVector point){
  double dist = ((tx[txindex].vect()-point.vect()).mag()/m)*((rx[rxindex].vect()-point.vect()).mag()/m);//here we've used the product of the distances as the radiated amplitude E~(E_0/R_1)/R_2. 

  Hep3Vector one=tx[txindex].vect()-point.vect();
  Hep3Vector two=point.vect()-rx[rxindex].vect();

  double angle_dependence=1.;
  Hep3Vector nhat(two.unit());
  Hep3Vector vert(0,1.,0), horiz(0,0,1.); 

  if(polarization=="vertical"){
    Hep3Vector pol=-one.unit().cross(one.unit().cross(vert));//for a dipole rad pattern

    angle_dependence = nhat.cross(nhat.cross(pol)).mag();
  }
  else{
    Hep3Vector pol=-one.unit().cross(one.unit().cross(horiz));

    angle_dependence = nhat.cross(nhat.cross(pol)).mag();
  }

  //double amplitude = ((tx_voltage/txEffectiveHeight)/dist)*angle_dependence;
  double amplitude = ((tx_voltage*txFactor)/dist)*angle_dependence;
  if(useAttnLengthFlag==1){
    double attn_dist = ((tx[txindex].vect()-point.vect()).mag())+((rx[rxindex].vect()-point.vect()).mag());//here the overall attenuation is just calculated over the full path length. 
    amplitude=amplitude*exp(-attn_dist/attnLength);
  }

  return amplitude;
}

double RadioScatter::getAmplitudeFromAt(double E_0,HepLorentzVector from, HepLorentzVector at){
  double dist=((tx[0].vect()-from.vect()).mag()/m)*((from.vect()-at.vect()).mag()/m);

  Hep3Vector one=tx[0].vect()-from.vect();
  Hep3Vector two=from.vect()-at.vect();

  double angle_dependence=1.;
  Hep3Vector nhat(two.unit());
  Hep3Vector vert(0,1.,0), horiz(0,0,1.);

  if(polarization=="vertical"){

    angle_dependence = nhat.cross(nhat.cross(vert)).mag();
  }
  else{

    angle_dependence = nhat.cross(nhat.cross(horiz)).mag();
  }

  double amplitude = (tx_voltage*m*m/dist)*angle_dependence;
  if(useAttnLengthFlag==1){
    amplitude=amplitude*exp(-dist/attnLength);
  }
  return amplitude;
}

double RadioScatter::getPhaseFromAt(HepLorentzVector from, HepLorentzVector at){
  double txtime = getTxTime(0,from, 0);//find retarted time
  double txphase = getTxPhase(txtime);//find phase at retarded time
  //time of full flight
  //  double tof = abs(rxtime-txtime);//time of flight
  //time of flight for zero lifetime(phase is fixed at interaction point)
  //  double tof = point.t()-(point.vect()-tx[index].vect()).mag()/c_light;
  double tof=from.t()-txtime;
  HepLorentzVector tx_pr=tx[0]-from, pr_rx = from-at;//make vectors
  //wave number addition
  Hep3Vector kvec1 = k*tx_pr.vect();
  Hep3Vector kvec2 = k*pr_rx.vect();
  Hep3Vector ktot = kvec1+kvec2;
  double kx = ktot.mag();
  //calculate compton effects
  //  double inv_omega_c = (1/omega)+(1/omega_e)*(1-cos(tx_pr.vect().unit().angle(pr_rx[index].vect().unit())));
  //omega_c = 1/inv_omega_c;
  //    cout<<txtime<<" "<<txphase<<" "<<rxtime<<endl;
  return ((kx) - omega*tof + txphase);
}


double RadioScatter::getTxPhase(double t_0){
  //    HepLorentzVector tx_pr = tx-point;
  double phase = omega*t_0 + phase0;
  return phase;
}

  double RadioScatter::getRxTime(HepLorentzVector point, Hep3Vector l, Hep3Vector j){
  double dist = j.mag()+l.mag();
  double time = point.t()+(dist/c_light);
    return time;
  }

  double RadioScatter::getRxTime(int index,HepLorentzVector point){

    Hep3Vector sep(rx[index].vect()-point.vect());
  double dist = sep.mag();
  double time = point.t()+(dist/c_light);

  return time;
  }

  double RadioScatter::getTxTime(int index,HepLorentzVector point, int direct=0){
  Hep3Vector distvec = point.vect()-tx[index].vect();
  double dist=0;
  if(direct==0){
    if(REFRACTION_FLAG==1){
      dist = findPathLengthWithRefraction(tx[index],point, tx_interface_dist[index]);
    }
    else{
      dist= distvec.mag();
    } 
  }
  else{
    dist = distvec.mag();
  }
    double time = point.t()-(dist/c_light);

    return time;
  }

 int RadioScatter::checkTxOn(double time){
  if(time>tx_on&&time<tx_off){
    return 1;
  }
  else{
    return 0;
  }
}

double RadioScatter::getDirectSignalPhase(int txindex, int rxindex, HepLorentzVector point){
  double dist = (tx[txindex]-rx[rxindex]).vect().mag();
  double t_0 = point.t()-(dist/c_light);
  return t_0*omega;
}
/*
using the refraction vectors (calculated in makeRays), calculate the wave vectors (including with modified k value in the medium) and the 
total time-of-flight and phase. the values n_rel, k_r, and c_light_r are the relative index of refraction, and 
the modified wavenumber and speed of light (for the medium) respectively  
*/
double RadioScatter::getRxPhase(HepLorentzVector point, Hep3Vector j1, Hep3Vector j2, Hep3Vector l1, Hep3Vector l2){
  
  //debug:check that snell's law is satisfied for the found paths.
  // cout<<(j1.z()/j1.mag())/(l1.z()/l1.mag())<<endl;
  // cout<<(j2.z()/j2.mag())/(l2.z()/l2.mag())<<endl<<endl;
  
  //calculate the time of flight using correct values for velocity
  //double  tof = j1.mag()/c_light + l1.mag()/c_light_r + j2.mag()/c_light + l2.mag()/c_light_r;

   //this is possibly incorrect. for a lifetime of zero, it may be correct to stop at the reflection point
  //and propagate the phase at the point of scattering. so the signal is a delta function with a fixed
  //phase (polarity).
  double  tof = j1.mag()/c_light + l1.mag()/c_light_r;
  
  double txphase = getTxPhase(point.t()-tof);//find phase at retarded time
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
double RadioScatter::getRxPhase(int txindex, int rxindex, HepLorentzVector point){
  double rxtime = getRxTime(rxindex,point);//find advanced time
  double txtime = getTxTime(txindex,point);//find retarted time
  double txphase = getTxPhase(txtime);//find phase at retarded time
  //time of full flight
  //  double tof = abs(rxtime-txtime);//time of flight
  //time of flight for zero lifetime(phase is fixed at interaction point)
  //  double tof = point.t()-(point.vect()-tx[index].vect()).mag()/c_light;
  double tof=point.t()-txtime;
  HepLorentzVector tx_pr=tx[txindex]-point, pr_rx = point-rx[rxindex];//make vectors
  //wave number addition
  Hep3Vector kvec1 = k*tx_pr.vect();
  Hep3Vector kvec2 = k*pr_rx.vect();
  Hep3Vector ktot = kvec1+kvec2;
  double kx = ktot.mag();
  //calculate compton effects
  //  double inv_omega_c = (1/omega)+(1/omega_e)*(1-cos(tx_pr.vect().unit().angle(pr_rx[index].vect().unit())));
  //omega_c = 1/inv_omega_c;
  //    cout<<txtime<<" "<<txphase<<" "<<rxtime<<endl;
  return ((kx) - omega*tof + txphase);
}




/*
simulate the direct signal that would be seen at the receiver for CW
 */
TH1F * RadioScatter::getDirectSignal(int txindex, int rxindex, const TH1F *in){
  TH1F *outhist=(TH1F*)in;
  HepLorentzVector point=rx[rxindex];
  double rx_amp, rx_ph, amp;
  Hep3Vector dist_vec = tx[txindex].vect()-rx[rxindex].vect();
  double dist = dist_vec.mag()/m;//in meters
  rx_amp = tx_voltage*txFactor/dist;
  //  cout<<"rx amp: "<<tx_voltage<<" "<<tx_voltage*m*m<<" "<<(tx[index].vect()-rx[index].vect()).mag()<<" "<<rx_amp;
  int size = in->GetNbinsX();
  int start = in->GetXaxis()->GetXmin();
  int end = in->GetXaxis()->GetXmax();

  for(int i=0;i<size;i++){
    point.setT(in->GetBinCenter(i));

    rx_ph = getDirectSignalPhase(txindex, rxindex, point);
    amp =rx_amp*rxFactor*cos(rx_ph);

    rx[rxindex].setT(in->GetBinCenter(i));
    if(checkTxOn(getTxTime(txindex, rx[rxindex], 1))==1){

      outhist->Fill(rx[rxindex].t(), (rx_amp*cos(rx_ph)));
    }
  }

  return outhist;
}

/*
the main function.
lifetime and screeing are included. 

for the case of refraction at a boundary:
 q1 and q2 are the direct path vectors between tx-ip and rx-ip respectively. 

j1, j2 are the vectors from tx-refraction point and rx-refraction point respectively, for their corresponding interfaces.
l1, l2 are the vectors from these refraction points to the interaction point.

calculate the phase, the amplitude, and the prefactors for cross-section, 

 */


double RadioScatter::makeRays(HepLorentzVector point, double e, double l, double e_i){
  if(RSCAT_HIST_RESIZE==false){
    makeTimeHist();
    RSCAT_HIST_RESIZE=true;
    RSCAT_HIST_DECLARE=true;
    
  }

  double rx_time, rx_amplitude, rx_phase, point_time, t_step=0.;
  double zz=point.z()*zscale;
  double tt=point.t()*tscale;

  //cout<<point.z()<<" ";
  point.setZ(zz);
  //  cout<<zscale<<" "<<point.z()<<endl;
  point.setT(tt);


  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){

      //would RF from the transmitter reached this point?
      if(checkTxOn(getTxTime(i,point))!=1)return 0;

      //calculate plasma freq and collison freq
      
      step_length=l;//to make our density approximation
      double n = e/e_i;//edeposited/ionization E
      double n_e =1.;
      
      if(step_length==0)return 0;
      
      
      //electron number density, using step length cube
      n_e = n*n_primaries/pow(step_length, 3);


      if(n_e==0)return 0;



      //rad scat as published
      double N_ice=3.2e19;//per mm^3;
      //the collisional cross section is some number x10^-16cm^-3.
      //NIST has a plot that depends upon the incident ionization energy
      //a value of 1e-16 is for 15eV ionization electron energy,
      //we use 3 for good measure.
      nu_col = sqrt(k_Boltzmann*(300)*kelvin/electron_mass_c2)*collisionalCrossSection*(N_ice);

      event.totNScatterers+=n;//track total number of scatterers

      //the full scattering amplitude pre-factor  
      //double prefactor = -rxEffectiveHeight*n*n_primaries*e_radius*omega/(pow(omega, 2)+pow(nu_col, 2));
      double prefactor = -rxFactor*n*n_primaries*e_radius*omega/(pow(omega, 2)+pow(nu_col, 2));

      //x position of charge w/r/t shower axis
      Hep3Vector vec=(point.vect()-event.position);
      double x = vec.mag()*abs(sin(vec.unit().angle(event.direction)));
      
      //only consider charges within 1 m (10 moliere radii, to accelerate sim). uses mm. so does c_light in alpha below, which is unitless
      double x_0=(abs(x)>1000?0:(1000-abs(x)));
      //plasma frequency
      double omega_p=sqrt(plasma_const*n_e)*1e-9;//in ns^-1

      //the screening term. as derived in paper
      double alpha= ((omega_p*omega_p)/(2.*c_light))*(nu_col/(omega*omega + nu_col*nu_col))*x_0;


      double attn_factor = exp(-alpha);

      prefactor=prefactor*attn_factor;
      
      HepLorentzVector point_temp=point;      
      //are we calculating in a refraction region?
      if(REFRACTION_FLAG==1){
	Hep3Vector  q1 = findRefractionPlane(tx[i], point);//make a plane where the refraction will happen
	Hep3Vector j1;
	j1.setZ(findRefractionPointZ(q1.x(), q1.z(), tx_interface_dist[i]));//find the refraction point in this plane on interface 
	Hep3Vector  q2 = findRefractionPlane(rx[j], point);
	Hep3Vector  j2;
	j2.setZ(findRefractionPointZ(q2.x(), q2.z(), rx_interface_dist[j]));
	j1.setX(tx_interface_dist[i]);
	j1.setY(0);
	j2.setX(rx_interface_dist[j]);
	j2.setY(0);
	Hep3Vector l1;
	l1.set(q1.x()-tx_interface_dist[i], 0., q1.z()-j1.z());
	Hep3Vector l2;
	l2.set(q2.x()-rx_interface_dist[j], 0., q2.z()-j2.z());


	point_time=point_temp.t();
	double point_time_end=point_time+lifetime;
	while(point_time<point_time_end){
	  //get the reflected signal amplitude and phase
	  rx_phase = getRxPhase(point_temp, j1, j2, l1, l2);
	  rx_amplitude = getRxAmplitude(j, point_temp, j1, j2, l1, l2);
 	  rx_time = getRxTime(point_temp, j2, l2);

	  double E_real= prefactor*rx_amplitude*(omega*cos(rx_phase)+nu_col*sin(rx_phase));
      	  double E_imag = prefactor*rx_amplitude*(-nu_col*cos(rx_phase)+omega*sin(rx_phase));
      
	  if(abs(E_real)<tx_voltage){//simple sanity check, probably not needed      
	    time_hist[i][j]->Fill(rx_time, E_real);
	    re_hist[i][j]->Fill(rx_time, E_real);
	    im_hist[i][j]->Fill(rx_time, E_imag);
	  }
	  point_time+=samplingperiod;
	  point_temp.setT(point_time);


	}
	
      }
      //assuming transmitter and interaction and receiver are in a medium with same refractive index.    
      else{
	point_time=point_temp.t();
	double point_time_end=point_time+lifetime;
	while(point_time<point_time_end){
	  rx_phase=getRxPhase(i,j,point_temp);
	  rx_amplitude=getRxAmplitude(i,j,point_temp);
	  rx_time=getRxTime(j,point_temp);

	  
	  double E_real= prefactor*rx_amplitude*(omega*cos(rx_phase)+nu_col*sin(rx_phase));
	  double E_imag = prefactor*rx_amplitude*(-nu_col*cos(rx_phase)+omega*sin(rx_phase));

	  if(abs(E_real)<tx_voltage){//simple sanity check
	    time_hist[i][j]->Fill(rx_time, E_real);
	    re_hist[i][j]->Fill(rx_time, E_real);
	    im_hist[i][j]->Fill(rx_time, E_imag);
	  }
	  point_time+=samplingperiod;
	  point_temp.setT(point_time);

	}
      }
    }
    return 1;
  }
}



//calculate screeing effects. at the moment, cannot be run in 'realtime' inside of GEANT4.
double RadioScatter::doScreening(TTree * tree, int entry){
  HepLorentzVector  point(0,0,0,0);
  HepLorentzVector * test=0;
  double n_e;
  tree->SetBranchAddress("point", &test);
  tree->SetBranchAddress("n_e", &n_e);
  tree->GetEntry(entry);
  point=*test;//assign "point" to this electron
  int entries = tree->GetEntries();
  double E_eff=0;
    for(int i=entry-40000;i<entry+40000;i++){
  //  for(int i=0;i<entries;i++){
    if(i>entries)break;
    if(i<0)continue;
    tree->GetEntry(i);
    //    cout<<i<<endl;
    double tof=(test->vect()-point.vect()).mag()/c_light;
    double t_a0=test->t()+tof;
    double t_a1=t_a0+lifetime;
    double t_b0=point.t();
    double t_b1=t_b0+lifetime;
    if(t_a0>t_b1||t_a1<t_b0)continue;
      double amp = getAmplitudeFromAt(tx_voltage, *test, point);
      double phase = getPhaseFromAt(*test, point);
      double ee = -(e_radius*amp*cos(phase)*n_e*n_primaries);//negative from polarity flip

	E_eff+=ee;

  }
  return E_eff;
}


// double RadioScatter::plotCausalCharges(TTree * tree, int entry){
//   HepLorentzVector  point(0,0,0,0);
//   HepLorentzVector * test=0;
//   double n_e;
//   tree->SetBranchAddress("point", &test);
//   tree->SetBranchAddress("n_e", &n_e);
//   tree->GetEntry(entry);
//   //TH3F *hist = new TH3F("asdf", "asdf", 40, 1, -1, 40, 1, -1, 40, 1, -1);
//   vector<double> xx, yy, zz;
//   point=*test;//assign "point" to this electron
//   int entries = tree->GetEntries();
//   double E_eff=0;
//     for(int i=0;i<entries;i++){
//     tree->GetEntry(i);
//     //    cout<<i<<endl;
//     double tof=(test->vect()-point.vect()).mag()/c_light;
//     double t_a0=test->t()+tof;
//     double t_a1=t_a0+lifetime;
//     double t_b0=point.t();
//     double t_b1=t_b0+lifetime;
//     if(t_a0>t_b1||t_a1<t_b0)continue;
//     xx.push_back(test->x());
//     yy.push_back(test->y());
//     zz.push_back(test->z());
//     double amp = getAmplitudeFromAt(tx_voltage, *test, point);
//     double phase = getPhaseFromAt(*test, point);
//     double ee = -(e_radius*amp*cos(phase)*n_e*n_primaries);//negative from polarity flip
//     //	if(ee!=0&&!isnan(ee)&&!isinf(ee))E_eff+=ee;
//     E_eff+=ee;
//     //	cout<<i-entry<<endl;
//     cout<<"amp: "<<amp<<" phase: "<<phase<<" E: "<<E_eff<<" "<<ee<<endl;
//     //cout<<point.x()<<" "<<test->x()<<" "<<endl;

// 	//    cout<<E_eff<<endl;
//   }
//     cout<<"ssdf"<<xx.size()<<endl;
//     TGraph2D *gr=new TGraph2D(xx.size(), &xx[0], &yy[0], &zz[0]);
//     gr->Draw("ap");
//     return E_eff;
// }










//this is hideous, but works....

//calculate the point of refraction in the z direction. n_rel is the relative index of refraction n2/n1
 double RadioScatter::findRefractionPointZ(double kx, double kz, double jx){
  double  v =n_rel;
  double jz=kz/2. + sqrt(pow(kz,2) - (2*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) + (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))))/2. - sqrt(2*pow(kz,2) - (4*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) - (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) - pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))) + (8*pow(kz,3) + (16*pow(jx,2)*kz*pow(v,2))/(-1 + pow(v,2)) - (8*kz*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(-1 + pow(v,2)))/(4.*sqrt(pow(kz,2) - (2*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)))/(3.*(-1 + pow(v,2))) + (pow(2,0.3333333333333333)*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),2))/(3.*(-1 + pow(v,2))*pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3) + sqrt(-4*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),6) + pow(108*pow(jx,4)*pow(kz,2)*pow(v,4)*(-1 + pow(v,2)) + 108*pow(jx,2)*pow(kz,4)*pow(v,2)*pow(-1 + pow(v,2),2) - 108*pow(jx,2)*pow(kz,2)*pow(v,2)*(-1 + pow(v,2))*(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2)) + 2*pow(-pow(jx,2) + 2*jx*kx - pow(kx,2) - pow(kz,2) + pow(jx,2)*pow(v,2) + pow(kz,2)*pow(v,2),3),2)),0.3333333333333333)/(3.*pow(2,0.3333333333333333)*(-1 + pow(v,2))))))/2;
return jz;
}

/*
analytic solution above requires that the tx and rx are in a plane. this function assumes 
this plane is in the plane containing the x-axis, interaction point, and either tx/rx[index]. 

function rotates the k vector to be in the kz plane and then sends this plane to the above function.
might be completely wrong. but resultant vectors satisfy snell's law.
 */

 Hep3Vector RadioScatter::findRefractionPlane(HepLorentzVector p1, HepLorentzVector p2){
  Hep3Vector k;
  k.setX(abs(p2.x()-p1.x()));

  k.setY(abs(p2.y()-p1.y()));
  k.setZ(abs(p2.z()-p1.z()));
  // cout<<endl<<k.x()<<" "<<k.y()<<" "<<k.z()<<endl;
  // cout<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<endl;
  // cout<<tx[index].x()<<" "<<tx[index].y()<<" "<<tx[index].z()<<endl;
  // cout<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<endl<<endl;
  double theta = atan(k.y()/k.z());
  k.rotateX(theta);
  return k;
}
/*
find the path length including refraction
 */
   double RadioScatter::findPathLengthWithRefraction(HepLorentzVector p1, HepLorentzVector p2, double interface_dist_x){
  Hep3Vector k = findRefractionPlane(p1,p2);
  double jz = findRefractionPointZ(k.x(), k.z(), interface_dist_x);
  Hep3Vector j(interface_dist_x, 0., jz);//from source point to refraction point
  Hep3Vector l(k.x()-interface_dist_x, 0., k.z()-jz);//from refraction point to interaction point
  double mag = j.mag()+l.mag();
  //  cout<<mag-k.mag()<<endl;
  return mag;
  
}

vector<vector<TH1F*>> RadioScatter::scaleHist(float num_events=1.){
  vector<vector<TH1F*>> outhist(ntx, vector<TH1F*>(nrx));
  for(int i=0;i<ntx;i++){
    for(int j=0;j<nrx;j++){
      time_hist[i][j]->Scale(1./num_events);

      re_hist[i][j]->Scale(1./num_events);
      im_hist[i][j]->Scale(1./num_events);

      if(includeCW_flag==1){
	outhist[i][j]=getDirectSignal(i, j, (const TH1F*)time_hist[i][j]);

      }
      else{
	outhist[i][j]=time_hist[i][j];
      }
    }
  }
  
  return outhist;
}

 void RadioScatter::writeToTextFile(){
  //  of<<
}

 void RadioScatter::printEventStats(){
  //  event.totNScatterers = event.totNScatterers*n_primaries;
  cout<<setprecision(10);
  cout<<"total number of scatterers: "<<  event.totNScatterers<<endl;
}


int RadioScatter::writeRun(float num_events, int debug){
  //this is a stupid check for multi-threaded mode,
  //will only write the run if there has indeed been a run
  if(event.totNScatterers==0){
    return 0;
  }
  fRunCounter++;
  TFile *f = (TFile *)gROOT->GetFile(output_file_name);
  //  TFile *f = ((TFile *)(gROOT->GetListOfFiles()->At(0)));
  TTree *t = (TTree*)f->Get("tree");
  event.nPrimaries = n_primaries;

  event.ntx=ntx;
  event.nrx=nrx;
  for(int i=0;i<ntx;i++){
    event.tx.push_back(tx[i]);
  }
  for(int i=0;i<nrx;i++){
    event.rx.push_back(rx[i]);
  }

  event.txGain=tx_gain;
  event.rxGain=rx_gain;
  event.eventHist = scaleHist(num_events);
  event.reHist = re_hist;
  event.imHist = im_hist;
  //event.testHist0=(TH1F*)testHist0->Clone();
  //event.testHist1=(TH1F*)testHist1->Clone();
  event.freq=frequency;
  //total number of electrons per shower * total primaries in the bunch * the number of events in the run. 
  event.totNScatterers = event.totNScatterers*n_primaries/num_events;
  vector<double> xvals, yvals;
  
  f->cd();
  
  
  t->Fill();
  if(debug==1){
    cout<<"The RadioScatter root file: "<<endl<<f->GetName()<<endl<<" has been filled. This is run number "<<fRunCounter<<"."<<endl;
  }
  
  cout<<"Run total N scatterers:"<<event.totNScatterers<<endl; 
  
  event.reset();
  
  return 1;
  
}

int RadioScatter::writeEvent(int debug){
  //this is a stupid check for multi-threaded mode,
  //will only write the run if there has indeed been a run
  
  fRunCounter++;
  
  TFile *f = ((TFile *)(gROOT->GetFile(output_file_name)));
  TTree *t = (TTree*)f->Get("tree");
  event.nPrimaries = n_primaries;
  
  event.ntx=ntx;
  event.nrx=nrx;
  event.txGain=tx_gain;
  event.rxGain=rx_gain;
  for(int i=0;i<ntx;i++){
    event.tx.push_back(tx[i]);
  }
  for(int i=0;i<nrx;i++){
    event.rx.push_back(rx[i]);
  }
  
  event.eventHist = time_hist;
  event.reHist = re_hist;
  event.imHist = im_hist;
  //  *event.testHist0=*testHist0;
  //*event.testHist1=*testHist1;
  event.freq=frequency;
  //total number of electrons per shower * total primaries in the bunch * the number of events in the run.
  
  event.totNScatterers = event.totNScatterers*n_primaries;

  f->cd();
  

  t->Fill();

  if(debug==1){
    cout<<"The RadioScatter root file: "<<endl<<f->GetName()<<endl<<" has been filled. This is run number "<<fRunCounter<<"."<<endl;
    
  }
  
  cout<<"Event total N scatterers:"<<event.totNScatterers<<endl; 

  event.reset();

  return 1;

}

int RadioScatter::makeSummary(TFile *f){
  RadioScatterEvent *rs = new RadioScatterEvent();
  TTree *intree = (TTree*)f->Get("tree");
  intree->SetBranchAddress("event", &rs);
  intree->GetEntry(0);
  TString fname = f->GetName();
  fname.ReplaceAll(".root", "_summary.root");
  TFile *outfile=new TFile(fname, "RECREATE"); 
  RSEventSummary *rss = new RSEventSummary(rs->ntx, rs->nrx);
  TTree *outtree= new TTree("sumTree", "tree of the things");
  outtree->Branch("summary", &rss);
  
  int entries = intree->GetEntries();
  for(int i=0;i<entries;i++){
    intree->GetEntry(i);
    rss->position=rs->position;
    rss->direction=rs->direction;
    rss->primaryParticleEnergy=rs->primaryParticleEnergy();
    rss->primaryEnergyG4=rs->primaryEnergy;
    rss->sampleRate=rs->sampleRate;
    rss->txVoltageV=rs->txVoltage;
    rss->txPowerW=rs->txPowerW;
    rss->freq=rs->freq;
    rss->totNScatterers=rs->totNScatterers;
    rss->tx=rs->tx;
    rss->rx=rs->rx;
    rss->ntx=rs->ntx;
    rss->nrx=rs->nrx;
    rss->weight=rs->weight;
    
    for(int j=0;j<rs->ntx;j++){
      for(int k=0;k<rs->nrx;k++){
	rss->peakV[j][k]=rs->peakV(j,k);
	rss->peakFreq[j][k]=rs->peakFreq(j,k);
	rss->duration[j][k]=rs->duration(j,k);
	rss->peakPowerW[j][k]=rs->peakPowerW(j,k);
	rss->integratedPower[j][k]=rs->integratedPower(j,k);
	rss->pathLengthM[j][k]=rs->pathLengthM(j,k);
	rss->rms[j][k]=rs->rms(j,k);
	rss->effectiveCrossSection[j][k]=rs->effectiveCrossSection(j,k);
      }
    }
    outfile->cd();
    outtree->Fill();
  }
  outfile->Write();
  cout<<"the summary file:"<<endl<<outfile->GetName()<<endl<<"has been written."<<endl;
  outfile->Close();
  return 0;
}

    
  void RadioScatter::close(){
    TFile *f=    ((TFile *)(gROOT->GetFile(output_file_name)));
    TString fname = f->GetName();

    f->Write();

    cout<<"The RadioScatter root file: "<<endl<<fname<<endl<<"has been written."<<endl;
    f->Close();

    f=TFile::Open(fname);
    if(MAKE_SUMMARY_FILE==1){
      makeSummary(f);
    }
    f->Close();

  }





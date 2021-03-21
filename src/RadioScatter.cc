/*
this is radioscatter. copyright s. prohira 

released under GPL3.
 
 */

#include "RadioScatter.hh"
// using namespace CLHEP;


/**
default constructor.
*/
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
       time_hist.push_back(std::vector<TH1F*>());
       re_hist.push_back(std::vector<TH1F*>());
       im_hist.push_back(std::vector<TH1F*>());
       event_graph.push_back(std::vector<TGraph*>());
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
     std::cout<<"histograms and graphs initialized"<<std::endl;

   }

   for(int i=0;i<ntx;i++){
     for(int j=0;j<nrx;j++){
       TVector3 dv = rx[j].Vect()-event.position;
       double dist = abs(dv.Mag());
       double time =(dist/c_light_r)-half_window; 
       
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


   if(RSCAT_POSITION_SET==true && RSCAT_DIRECTION_SET==true){
     //initialize the variables for calculating frequency/direction
     event.delta.resize(ntx, vector<double>(nrx, 0.));
     event.beta.resize(ntx, vector<double>(nrx, 0.));
     event.doppler.resize(ntx, vector<double>(nrx, 0.));
     event.angTVP.resize(ntx, vector<double>(nrx, 0.));
     event.angRVP.resize(ntx, vector<double>(nrx, 0.));
     event.angTVR.resize(ntx, vector<double>(nrx, 0.));
     auto d = event.direction;
     //cout<<d.X()<<" "<<d.Y()<<" "<<d.Z()<<endl;
     
     auto p=event.position;
     //cout<<p.X()<<" "<<p.Y()<<" "<<p.Z()<<endl;
     for(int i=0;i<ntx;i++){
       auto p_tx=tx[i].Vect()-p;
       for(int j=0;j<nrx;j++){
	 auto p_rx=rx[j].Vect()-p;
	 auto tx_rx=rx[j].Vect()-tx[i].Vect();
	 auto beta=p_tx.Angle(p_rx);
	 auto alpha=tx_rx.Angle(p_tx);
	 auto gamma=pi-alpha-(beta/2.);
	 auto c=p_tx.Mag();
	 auto a=(c*sin(alpha))/sin(gamma);
	 auto b = (c*sin(beta/2.))/sin(gamma);
	 auto bVec=tx_rx;
	 bVec.SetMag(b);
	 auto B1=tx[i].Vect()+bVec;
	 auto B2=B1-p;
	 auto delta=d.Angle(B2);

	 auto q=p_tx.Cross(p_rx);
	 B2=p_tx;
	 B2.Rotate(beta/2., q);
	 delta=d.Angle(B2);
	 event.delta[i][j]=delta;
	 event.beta[i][j]=beta;
	 event.angTVP[i][j]=d.Angle(p_tx);
	 event.angRVP[i][j]=d.Angle(p_rx);
	 event.angTVR[i][j]=beta;//redundant
       }
     }
   }
   else{
     cout<<"WARNING:::::::::::::::::::::::::::::::::::"<<endl<<endl<<"You have not set the position and direction of the primary. The bistatic angles (and other things) will be undefined."<<endl<<endl<<"WARNING:::::::::::::::::::::::::::::::::::"<<endl;
   }
 }


void RadioScatter::setMakeSummary(double val){
  MAKE_SUMMARY_FILE=(int)val;
}
void RadioScatter::setRecordWindowLength(double nss){
  half_window= nss/2;
}

void RadioScatter::setCalculateUsingAttnLength(double val){
  useAttnLengthFlag=(int)val;
  if((int)val==1){
    std::cout<<"calculating using attenuation length. (Barwick et al)"<<std::endl;
  }
}

int RadioScatter::setNTx(double n){
  ntx=(int)n;
  TLorentzVector vec;
  if(tx.size()==ntx)return 1;
  for(int i=0;i<ntx;i++){
    tx.push_back(vec);
    tx_interface_dist.push_back(0);
  }
  return 1;
}
int RadioScatter::setNRx(double n){
  nrx=(int)n;
  TLorentzVector vec;
  if(rx.size()==nrx)return 1;
  for(int i=0;i<nrx;i++){
    rx.push_back(vec);
    rx_interface_dist.push_back(0);
  }
  return 1;
}
void RadioScatter::setTxPos(double xin, double yin, double zin, int index){
    tx[index].SetX(xin);
    tx[index].SetY(yin);
    tx[index].SetZ(zin);
    RADIOSCATTER_INIT=false;//resize the output histogram
  }
void RadioScatter::setRxPos(double xin, double yin, double zin, int index){
    rx[index].SetX(xin);
    rx[index].SetY(yin);
    rx[index].SetZ(zin);
    RADIOSCATTER_INIT=false;
  }
void RadioScatter::setTxPos(TVector3 in){
  setTxPos(in, TX_ITERATOR);
  TX_ITERATOR+=1;
}
void RadioScatter::setTxPos(TVector3 in, int index){
  tx[index].SetX(in.X());
  tx[index].SetY(in.Y());
  tx[index].SetZ(in.Z());
  RADIOSCATTER_INIT=false;
  }
void RadioScatter::setRxPos(TVector3 in){
  setRxPos(in, RX_ITERATOR);
  RX_ITERATOR+=1;
}
void RadioScatter::setRxPos(TVector3 in, int index){
  rx[index].SetX(in.X());
  rx[index].SetY(in.Y());
  rx[index].SetZ(in.Z());
  RADIOSCATTER_INIT=false;
  }
  void RadioScatter::setTxVals(double f, double power=1., double gain=1.){
    frequency = f*GHz;
    omega = frequency*twoPi;
    period = 1./omega;
    lambda = (c_light/frequency)/m;
    k = (omega/c_light);
    tx_gain = gain;
    tx_voltage = sqrt(power*50.);
  }
 void RadioScatter::setTxFreq(double f){
   TX_FREQ_SET=1;
   frequency = f*GHz;
  omega = frequency*twoPi;
  period = 1./omega;
  lambda = (c_light/f)/m;
  k = (omega/c_light);
  
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
  std::cout<<"tx frequency: "<<f<<std::endl;
  std::cout<<"omega: "<<omega<<std::endl;
  std::cout<<"lambda: "<<lambda<<std::endl;

  std::cout<<"attn length (mm): "<<attnLength<<std::endl;
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

    std::cout<<"receiver gain: "<<gain<<" dB, ("<<gg<<" linear)"<<std::endl;
    //from rice paper
    rxEffectiveHeight=sqrt(lambda*lambda*rx_gain/(480.*pi*pi));
    //from wikipedia antenna factor
    rxEffectiveHeight=lambda*sqrt(rx_gain)/9.73;
    //first principles losses due to spherical spreading (Friis)         
    rxFactor=sqrt(rx_gain*lambda/(4.*pi));
    
    std::cout<<"effective height: "<<rxEffectiveHeight<<std::endl;
  }
}
void RadioScatter::setTxGain(double gain){
  setTransmitterGain(gain);
}
void RadioScatter::setTransmitterGain(double gain){
  if(TX_FREQ_SET==0){
    std::cout<<"WARNING:antenna gain set before transmit frequency, effective height will be incorrect. please set tx freq before setting antenna gain"<<std::endl;
  }
  double gg=pow(10., gain/10.);
  tx_gain=gg;
  if(TX_GAIN_SET!=1){
    TX_GAIN_SET==1;

    std::cout<<"transmitter voltage: "<<tx_voltage<<" transmitter gain: "<<gain<<" dB, ("<<gg<<" linear)"<<std::endl;

    //from RICE but with the right placement of the gain
    txEffectiveHeight=sqrt(lambda*lambda/(tx_gain*480.*pi*pi));
    //from wikipedia antenna factor but with the right placement of gain
    txEffectiveHeight=lambda/(9.73*sqrt(tx_gain));
    //first principles losses due to spherical spreading (Friis)         
    txFactor=sqrt(tx_gain*lambda/(4.*pi));
    std::cout<<"effective height: "<<txEffectiveHeight<<std::endl;
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

  std::cout<<"n primaries: "<<n_primaries<<std::endl;
  NPRIMARIES_SET=1;

}
void RadioScatter::setWeight(double val){
  event.weight=val;
}
void RadioScatter::setCrossSection(double val){
  collisionalCrossSection=val;
}

void RadioScatter::setTargetEnergy(double e){
  event.targetEnergy=e;
  cout<<"target simulation energy: "<<event.targetEnergy/TUtilRadioScatter::GeV<<" GeV"<<endl;
  TARGET_ENERGY_SET=1;
  NPRIMARIES_SET=0;
  ENERGY_SCALING_SET=0;//it is possible that the target energy has changed (from a previous value). re-set the scaling.
}
int RadioScatter::setScaleByEnergy(double val){
  if((int)val==1){
    SCALE_BY_ENERGY=1;
  }
}
  
int RadioScatter::scaleByEnergy(){
  if(NPRIMARIES_SET==1&&PRIMARY_ENERGY_SET==1){
    if(TARGET_ENERGY_SET!=1){
      setTargetEnergy(event.primaryEnergy*event.nPrimaries);
    }
        
    zscale = (3.*log10(event.targetEnergy)+6.)/((log10(event.primaryEnergy)*3.)+6);
    
    tscale = (10.*log10(event.targetEnergy)+22.)/((log10(event.primaryEnergy)*10.)+22);
    std::cout<<"scaling activated. zscale="<<zscale<<" , tscale="<<tscale<<std::endl;

    //cout<<"__________________"<<endl<<" "<<event.primaryEnergy<<" "<<event.targetEnergy<<" "<<event.nPrimaries<<endl<<"______________"<<endl;
    
    ENERGY_SCALING_SET=1;
    return 1;
  }
  // else {
  //   cout<<
  //     }
}
void RadioScatter::setPrimaryEnergy(double e){
  event.primaryEnergy=e;
  cout<<"primary simulated particle energy: "<<event.primaryEnergy/TUtilRadioScatter::GeV<<" GeV"<<endl;
  //cout<<"__________________"<<endl<<e<<" "<<event.nPrimaries<<endl<<"______________"<<endl;
  PRIMARY_ENERGY_SET=1;//the primary energy has been set

}
void RadioScatter::setPrimaryDirection(TVector3 d){
  // if(RSCAT_POSITION_SET==0){
  //   cout<<"position not set yet! please define vertex position before setting the direction."<<endl;
  //   exit(0);
  // }
  event.direction=d;
  // event.delta.resize(ntx, vector<double>(nrx, 0.));
  // event.beta.resize(ntx, vector<double>(nrx, 0.));
  // event.doppler.resize(ntx, vector<double>(nrx, 0.));

  // auto p=event.position;
  // for(int i=0;i<ntx;i++){
  //   auto p_tx=tx[i].Vect()-p;
  //   for(int j=0;j<nrx;j++){
  //     auto p_rx=rx[j].Vect()-p;
  //     auto tx_rx=rx[j].Vect()-tx[i].Vect();
  //     auto beta=p_tx.Angle(p_rx);
  //     auto alpha=tx_rx.Angle(p_tx);
  //     auto gamma=pi-alpha-(beta/2.);
  //     auto c=p_tx.Mag();
  //     auto a=(c*sin(alpha))/sin(gamma);
  //     auto b = (c*sin(beta/2.))/sin(gamma);
  //     auto bVec=tx_rx;
  //     bVec.SetMag(b);
  //     auto B1=tx[i].Vect()+bVec;
  //     auto B2=B1-p;
  //     auto delta=d.Angle(B2);
  //     event.delta[i][j]=delta;
  //     event.beta[i][j]=beta;
  //   }
  // }

  RADIOSCATTER_INIT=false;
  RSCAT_DIRECTION_SET=true;
}

void RadioScatter::setPrimaryPosition(TVector3 p){
  event.position=p;
  
  // for(int i=0;i<ntx;i++){
  //   auto p_tx=tx[i].Vect()-p;
  //   for(int j=0;j<nrx;j++){
  //     auto p_rx=rx[j].Vect()-p;
  //     auto tx_rx=rx[j].Vect()-tx[i].Vect();
  //     auto beta=p_tx.Angle(p_rx);
  //     event.beta[i][j]=beta;
  //   }
  // }
  RADIOSCATTER_INIT=false;
  RSCAT_POSITION_SET=true;
}
void RadioScatter::setPolarization(TVector3 p){//const char * p){
  event.polarization=p;

  std::cout<<"polarization: "<<event.polarization.X()<<" "<<event.polarization.Y()<<" "<<event.polarization.Z()<<std::endl;

}

void RadioScatter::setPlasmaLifetime(double l){
  lifetime=l;
}
void RadioScatter::setRxVals(double s=1., double gain=1.){
    samplerate = s*ns;
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
void RadioScatter::setIndexOfRefraction(double iof){
  n_rel=iof;
  k=omega/(c_light/iof);
  c_light_r=c_light/iof;

}
 void RadioScatter::setRxSampleRate(double rate){
  samplerate=rate*ns;
  samplingperiod = 1./samplerate;
  event.sampleRate=samplerate;
}
void RadioScatter::setTxInterfaceDistX(double dist, int index){
    tx_interface_dist[index] = abs(dist);
    BOUNDARY_FLAG=1;
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
TVector3 RadioScatter::getTxPos(int index){
  return tx[index].Vect();
}
 TVector3 RadioScatter::getRxPos(int index){
  return rx[index].Vect();
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
  double RadioScatter::getTxAmplitude(int index,TLorentzVector point){
    double gain = getTxGain(point.Vect().Theta());
    double amplitude = tx_voltage/txEffectiveHeight;
    return amplitude;
  }


/*
use the calculated refraction vectors (from makeRays()) to sort out the correct amplitude.
 */

 double RadioScatter::getRxAmplitude(int index,TLorentzVector point, TVector3 j1, TVector3 j2, TVector3 l1, TVector3 l2){
   double dist = ((j1.Mag()+j2.Mag())/m)*((l1.Mag()+l2.Mag())/m);
  //refraction things:


  double alpha1 = atan(j1.Z()/j1.X());
  double alpha_prime1 = atan(l1.Z()/l1.X());
  double alpha2 = atan(j2.Z()/j2.X());
  double alpha_prime2 = atan(l2.Z()/l2.X());

  //E for polarization parallel to plane of incidence
  double E1_para = (2.*cos(alpha1))/(cos(alpha1)+(n_rel*cos(alpha_prime1)));
  double E1_perp = (2.*n_rel*cos(alpha1))/(n_rel*n_rel*cos(alpha1)+(n_rel*cos(alpha_prime1)));

  //E for polarization perpendicular to plane of incidence
  double E2_para = (2.*cos(alpha2))/(cos(alpha2)+(n_rel*cos(alpha_prime2)));
  double E2_perp = (2.*n_rel*cos(alpha2))/(n_rel*n_rel*cos(alpha2)+(n_rel*cos(alpha_prime2)));

  //find angle between plane of incidence and polarization vector
  double theta = atan(point.Y()/point.Z());
  double angle_dependence=1.;

  
  TVector3 nhat((point-rx[index]).Vect().Unit());

  //TVector3 vert(0,1.,0), horiz(0,0,1.);
  TVector3 vert(1.,1.,1.), horiz(1.,1.,1.); 
  //  l1plane.setTheta(0);
  //  l1plane.setPhi(0);

  //NB unused  
  if(event.polarization.Z()==1){
    //refraction angle change
    theta = theta+(pi/2.);

    //angle_dependence = vert.Cross(nhat).Mag();
    angle_dependence = nhat.Cross(nhat.Cross(vert)).Mag();

  }
  else{

    //angle_dependence = horiz.Cross(nhat).Mag();
    angle_dependence = nhat.Cross(nhat.Cross(horiz)).Mag();

  }
  double amp1 = sqrt(pow(E1_para*cos(theta), 2)+pow(E1_perp*sin(theta), 2));
  double amp2 = sqrt(pow(E2_para*cos(theta), 2)+pow(E2_perp*sin(theta), 2));


  //double amplitude = ((tx_voltage/txEffectiveHeight)/dist)*amp1*amp2*angle_dependence;
  double amplitude = ((tx_voltage*txFactor)/dist)*amp1*amp2*angle_dependence;

  if(useAttnLengthFlag==1){
    double attn_dist = (j1.Mag()+j2.Mag())+(l1.Mag()+l2.Mag());
    amplitude=amplitude*exp(-attn_dist/attnLength);
  }
  return amplitude;
}

//non-refracted amplitude
double RadioScatter::getRxAmplitude(int txindex,int rxindex, TLorentzVector point){
  double dist = ((tx[txindex].Vect()-point.Vect()).Mag()/m)*((rx[rxindex].Vect()-point.Vect()).Mag()/m);//here we've used the product of the distances as the radiated amplitude E~(E_0/R_1)/R_2. 

  TVector3 one=tx[txindex].Vect()-point.Vect();
  TVector3 two=point.Vect()-rx[rxindex].Vect();

  double angle_dependence=1.;
  TVector3 nhat(two.Unit());
  


  TVector3 pol=-one.Unit().Cross(one.Unit().Cross(event.polarization.Unit()));//for a dipole rad pattern

  angle_dependence = nhat.Cross(nhat.Cross(pol)).Mag();

  //double amplitude = ((tx_voltage/txEffectiveHeight)/dist)*angle_dependence;
  double amplitude = ((tx_voltage*txFactor)/dist)*angle_dependence;
  if(useAttnLengthFlag==1){
    double attn_dist = ((tx[txindex].Vect()-point.Vect()).Mag())+((rx[rxindex].Vect()-point.Vect()).Mag());//here the overall attenuation is just calculated over the full path length. 
    amplitude=amplitude*exp(-attn_dist/attnLength);
  }

  return amplitude;
}

double RadioScatter::getAmplitudeFromAt(double E_0,TLorentzVector from, TLorentzVector at){
  double dist=((tx[0].Vect()-from.Vect()).Mag()/m)*((from.Vect()-at.Vect()).Mag()/m);

  TVector3 one=tx[0].Vect()-from.Vect();
  TVector3 two=from.Vect()-at.Vect();

  double angle_dependence=1.;
  TVector3 nhat(two.Unit());
  TVector3 vert(0,1.,0), horiz(0,0,1.);



  angle_dependence = nhat.Cross(nhat.Cross(event.polarization.Unit())).Mag();
  
  double amplitude = (tx_voltage*m*m/dist)*angle_dependence;
  if(useAttnLengthFlag==1){
    amplitude=amplitude*exp(-dist/attnLength);
  }
  return amplitude;
}

double RadioScatter::getPhaseFromAt(TLorentzVector from, TLorentzVector at){
  double txtime = getTxTime(0,from, 0);//find retarted time
  double txphase = getTxPhase(txtime);//find phase at retarded time
  //time of full flight
  //  double tof = abs(rxtime-txtime);//time of flight
  //time of flight for zero lifetime(phase is fixed at interaction point)
  //  double tof = point.T()-(point.Vect()-tx[index].Vect()).Mag()/c_light;
  double tof=from.T()-txtime;
  TLorentzVector tx_pr=tx[0]-from, pr_rx = from-at;//make vectors
  //wave number addition
  TVector3 kvec1 = k*tx_pr.Vect();
  TVector3 kvec2 = k*pr_rx.Vect();
  TVector3 ktot = kvec1+kvec2;
  double kx = ktot.Mag();
  //calculate compton effects
  //  double inv_omega_c = (1/omega)+(1/omega_e)*(1-cos(tx_pr.Vect().Unit().Angle(pr_rx[index].Vect().Unit())));
  //omega_c = 1/inv_omega_c;
  //    std::cout<<txtime<<" "<<txphase<<" "<<rxtime<<std::endl;
  return ((kx) - omega*tof + txphase);
}


double RadioScatter::getTxPhase(double t_0){
  //    TLorentzVector tx_pr = tx-point;
  double phase = omega*t_0 + phase0;
  return phase;
}

  double RadioScatter::getRxTime(TLorentzVector point, TVector3 l, TVector3 j){
  double dist = j.Mag()+l.Mag();
  double time = point.T()+abs(dist/c_light_r);
    return time;
  }

  double RadioScatter::getRxTime(int index,TLorentzVector point){

    TVector3 sep(rx[index].Vect()-point.Vect());
  double dist = sep.Mag();
  double time = point.T()+abs(dist/c_light_r);

  return time;
  }

  double RadioScatter::getTxTime(int index,TLorentzVector point, int direct=0){
  TVector3 distvec = point.Vect()-tx[index].Vect();
  double dist=0;
  if(direct==0){
    if(BOUNDARY_FLAG==1){
      dist = findPathLengthWithRefraction(tx[index],point, tx_interface_dist[index]);
    }
    else{
      dist= distvec.Mag();
    } 
  }
  else{
    dist = distvec.Mag();
  }
    double time = point.T()-abs(dist/c_light_r);

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

double RadioScatter::getDirectSignalPhase(int txindex, int rxindex, TLorentzVector point){
  double dist = (tx[txindex]-rx[rxindex]).Vect().Mag();
  double t_0 = point.T()-(dist/c_light);
  return t_0*omega;
}
/*
using the refraction vectors (calculated in makeRays), calculate the wave vectors (including with modified k value in the medium) and the 
total time-of-flight and phase. the values n_rel, k_r, and c_light_r are the relative index of refraction, and 
the modified wavenumber and speed of light (for the medium) respectively  
*/
double RadioScatter::getRxPhase(TLorentzVector point, TVector3 j1, TVector3 j2, TVector3 l1, TVector3 l2){
  
  //debug:check that snell's law is satisfied for the found paths.
  // std::cout<<(j1.Z()/j1.Mag())/(l1.Z()/l1.Mag())<<std::endl;
  // std::cout<<(j2.Z()/j2.Mag())/(l2.Z()/l2.Mag())<<std::endl<<std::endl;
  
  //calculate the time of flight using correct values for velocity
  //double  tof = j1.Mag()/c_light + l1.Mag()/c_light_r + j2.Mag()/c_light + l2.Mag()/c_light_r;

   //this is possibly incorrect. for a lifetime of zero, it may be correct to stop at the reflection point
  //and propagate the phase at the point of scattering. so the signal is a delta function with a fixed
  //phase (polarity).
  double  tof = j1.Mag()/c_light + l1.Mag()/c_light_r;

  double txphase = getTxPhase(point.T()-tof);//find phase at retarded time
  //wave vector calculation, with correct phase velocity
  TVector3  kvec1 = k*j1;
  TVector3  kvec2 = k_r*l1;
  TVector3  kvec3 = k_r*l2;
  TVector3  kvec4 = k*j2;
  TVector3  ktot = kvec1+kvec2+kvec3+kvec4;
  double  kx = kvec1.Mag()+kvec2.Mag()+kvec3.Mag()+kvec4.Mag();

  return ((kx) - omega*tof + txphase);
  }

//non-refraction phase finder
double RadioScatter::getRxPhase(int txindex, int rxindex, TLorentzVector point){
  double rxtime = getRxTime(rxindex,point);//find advanced time
  double txtime = getTxTime(txindex,point);//find retarted time
  double txphase = getTxPhase(txtime);//find phase at retarded time
  //time of full flight
  //  double tof = abs(rxtime-txtime);//time of flight
  //time of flight for zero lifetime(phase is fixed at interaction point)
  //  double tof = point.T()-(point.Vect()-tx[index].Vect()).Mag()/c_light;
  double tof=point.T()-txtime;
  TLorentzVector tx_pr=tx[txindex]-point, pr_rx = point-rx[rxindex];//make vectors
  //wave number addition
  TVector3 kvec1 = k*tx_pr.Vect();
  TVector3 kvec2 = k*pr_rx.Vect();
  TVector3 ktot = kvec1+kvec2;
  double kx = ktot.Mag();
  //calculate compton effects
  //  double inv_omega_c = (1/omega)+(1/omega_e)*(1-cos(tx_pr.Vect().Unit().Angle(pr_rx[index].Vect().Unit())));
  //omega_c = 1/inv_omega_c;
  //    std::cout<<txtime<<" "<<txphase<<" "<<rxtime<<std::endl;
  return ((kx) - omega*tof + txphase);
}




/*
simulate the direct signal that would be seen at the receiver for CW
 */
TH1F * RadioScatter::getDirectSignal(int txindex, int rxindex, const TH1F *in){
  TH1F *outhist=(TH1F*)in;
  TLorentzVector point=rx[rxindex];
  double rx_amp, rx_ph, amp;
  TVector3 dist_vec = tx[txindex].Vect()-rx[rxindex].Vect();
  double dist = dist_vec.Mag()/m;//in meters
  rx_amp = tx_voltage*txFactor/dist;
  //  std::cout<<"rx amp: "<<tx_voltage<<" "<<tx_voltage*m*m<<" "<<(tx[index].Vect()-rx[index].Vect()).Mag()<<" "<<rx_amp;
  int size = in->GetNbinsX();
  int start = in->GetXaxis()->GetXmin();
  int end = in->GetXaxis()->GetXmax();

  for(int i=0;i<size;i++){
    point.SetT(in->GetBinCenter(i));

    rx_ph = getDirectSignalPhase(txindex, rxindex, point);
    amp =rx_amp*rxFactor*cos(rx_ph);

    rx[rxindex].SetT(in->GetBinCenter(i));
    if(checkTxOn(getTxTime(txindex, rx[rxindex], 1))==1){

      outhist->Fill(rx[rxindex].T(), (rx_amp*cos(rx_ph)));
    }
  }

  return outhist;
}

/* Function for use to turn the analytic raytracing on or off*/
void RadioScatter::setUseRayTracing(bool flag){
  USE_RAYTRACING=flag;
}

/* This function is used by the analytical raytracer to sort out ray parameters for the two ray solutions out the three possible ones */
double* RadioScatter::getPathAndTimeOfRays(double TxRaySolPar[3][5], double RxRaySolPar[3][5]){
  
  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[4];
  double TxRays[4]={0, 0, 0, 0};
  double RxRays[4]={0, 0, 0, 0};

  /************************Tx Ray Parameters**********************/
  if(TxRaySolPar[0][1]!=0 && TxRaySolPar[1][1]!=0){
    ////****************For Direct and Reflected Ray btw Tx and Shower***********************
    TxRays[0]=TxRaySolPar[0][3];
    TxRays[1]=TxRaySolPar[1][3];
    TxRays[2]=TxRaySolPar[0][2];
    TxRays[3]=TxRaySolPar[1][2];
  }

  if(TxRaySolPar[0][1]!=0 && TxRaySolPar[2][1]!=0){
    ////****************For Direct and Refracted Ray btw Tx and Shower***********************
    TxRays[0]=TxRaySolPar[0][3];
    TxRays[1]=TxRaySolPar[2][3];
    TxRays[2]=TxRaySolPar[0][2];
    TxRays[3]=TxRaySolPar[2][2];
  }

  if(TxRaySolPar[2][1]!=0 && TxRaySolPar[1][1]!=0){
    ////****************For Refracted and Reflected Ray btw Tx and Shower***********************
    TxRays[0]=TxRaySolPar[2][3];
    TxRays[1]=TxRaySolPar[1][3];
    TxRays[2]=TxRaySolPar[2][2];
    TxRays[3]=TxRaySolPar[1][2];
  }

  /************************Rx Ray Parameters**********************/
  if(RxRaySolPar[0][1]!=0 && RxRaySolPar[1][1]!=0){
    ////****************For Direct and Reflected Ray btw Rx and Shower***********************
    RxRays[0]=RxRaySolPar[0][3];
    RxRays[1]=RxRaySolPar[1][3];
    RxRays[2]=RxRaySolPar[0][2];
    RxRays[3]=RxRaySolPar[1][2];
  }

  if(RxRaySolPar[0][1]!=0 && RxRaySolPar[2][1]!=0){
    ////****************For Direct and Refracted Ray btw Rx and Shower***********************
    RxRays[0]=RxRaySolPar[0][3];
    RxRays[1]=RxRaySolPar[2][3];
    RxRays[2]=RxRaySolPar[0][2];
    RxRays[3]=RxRaySolPar[2][2];
  }

  if(RxRaySolPar[2][1]!=0 && RxRaySolPar[1][1]!=0){
    ////****************For Refracted and Reflected Ray btw Rx and Shower***********************
    RxRays[0]=RxRaySolPar[2][3];
    RxRays[1]=RxRaySolPar[1][3];
    RxRays[2]=RxRaySolPar[2][2];
    RxRays[3]=RxRaySolPar[1][2];
  }

  output[0]=TxRays[2];
  output[1]=TxRays[3];
  output[2]=RxRays[2];
  output[3]=RxRays[3];
  
  return output;    
}

/* This function calls the analytical raytracer and calculates propogation times, optical path lengths, launch angles and recieve angles for all the possible ray paths in the given Tx->Shower->Rx configuration */
double *RadioScatter::rayTrace(TLorentzVector Tx, TLorentzVector Rx, TVector3 Shwr){

  //////////////////////////////////////////////////////////////////////
  //divide by meter. m = 1000 *mm
  double Tx_x=Tx.X()/m;
  double Tx_y=Tx.Y()/m;
  double Tx_z=Tx.Z()/m;

  double Rx_x=Rx.X()/m;
  double Rx_y=Rx.Y()/m;
  double Rx_z=Rx.Z()/m;

  double Shwr_x=Shwr.X()/m;
  double Shwr_y=Shwr.Y()/m;
  double Shwr_z=Shwr.Z()/m;
  
  double startpoint=0;////always zero
  double Tx2ShwrDist=sqrt(pow(Tx_x-Shwr_x,2)+ pow(Tx_y-Shwr_y,2));
  double Rx2ShwrDist=sqrt(pow(Rx_x-Shwr_x,2)+ pow(Rx_y-Shwr_y,2));
  
  /* These Two arrays will be filled up with the ray parameters */
  /* for each of the 3 possible rays we have 4 elements*/
  /* element 0 is LaunchAngle */
  /* element 1 is RecieveAngle  */
  /* element 2 is PropogationTime  */
  /* element 3 is PropogationDistance  */
  /* element 4 is the IncidentAngleInIce for Reflected ray. For the other two rays this number goes to zero.  */
  double TxRaySolPar[3][5];
  double RxRaySolPar[3][5];
  
  double* GetTx2ShwrRays=IceRayTracing::IceRayTracing(startpoint,Tx_z,Tx2ShwrDist,Shwr_z);
 
  if(GetTx2ShwrRays[6]!=0){
    //cout<<"we got a direct ray!"<<endl;
    TxRaySolPar[0][0]=GetTx2ShwrRays[0];
    TxRaySolPar[0][1]=GetTx2ShwrRays[6];
    TxRaySolPar[0][2]=GetTx2ShwrRays[3]*s;
    TxRaySolPar[0][3]=GetTx2ShwrRays[3]*IceRayTracing::c_light_ms*m;
    TxRaySolPar[0][4]=0;
  }

  if(GetTx2ShwrRays[7]!=0){
    //cout<<"we got a reflected ray!"<<endl;
    TxRaySolPar[1][0]=GetTx2ShwrRays[1];
    TxRaySolPar[1][1]=GetTx2ShwrRays[7];
    TxRaySolPar[1][2]=GetTx2ShwrRays[4]*s;
    TxRaySolPar[1][3]=GetTx2ShwrRays[4]*IceRayTracing::c_light_ms*m;
    TxRaySolPar[1][4]=GetTx2ShwrRays[11];
  }
  if(GetTx2ShwrRays[8]!=0){
    //cout<<"we got a refracted ray!"<<endl;
    TxRaySolPar[2][0]=GetTx2ShwrRays[2];
    TxRaySolPar[2][1]=GetTx2ShwrRays[8];
    TxRaySolPar[2][2]=GetTx2ShwrRays[5]*s;
    TxRaySolPar[2][3]=GetTx2ShwrRays[5]*IceRayTracing::c_light_ms*m;
    TxRaySolPar[2][4]=0;
  }
  delete []GetTx2ShwrRays;

  if(TxRaySolPar[0][1]!=0 || TxRaySolPar[1][1]!=0 || TxRaySolPar[2][1]!=0){
    double* GetRx2ShwrRays=IceRayTracing::IceRayTracing(startpoint,Rx_z,Rx2ShwrDist,Shwr_z);

    if(GetRx2ShwrRays[6]!=0){
      //cout<<"we got a direct ray!"<<endl;
      RxRaySolPar[0][0]=GetRx2ShwrRays[0];
      RxRaySolPar[0][1]=GetRx2ShwrRays[6];
      RxRaySolPar[0][2]=GetRx2ShwrRays[3]*s;
      RxRaySolPar[0][3]=GetRx2ShwrRays[3]*IceRayTracing::c_light_ms*m;
      RxRaySolPar[0][4]=0;
    }
    if(GetRx2ShwrRays[7]!=0){
      //cout<<"we got a reflected ray!"<<endl;
      RxRaySolPar[1][0]=GetRx2ShwrRays[1];
      RxRaySolPar[1][1]=GetRx2ShwrRays[7];
      RxRaySolPar[1][2]=GetRx2ShwrRays[4]*s;
      RxRaySolPar[1][3]=GetRx2ShwrRays[4]*IceRayTracing::c_light_ms*m;
      RxRaySolPar[1][4]=GetRx2ShwrRays[11];
    }
    if(GetRx2ShwrRays[8]!=0){
      //cout<<"we got a refracted ray!"<<endl;
      RxRaySolPar[2][0]=GetRx2ShwrRays[2];
      RxRaySolPar[2][1]=GetRx2ShwrRays[8];
      RxRaySolPar[2][2]=GetRx2ShwrRays[5]*s;
      RxRaySolPar[2][3]=GetRx2ShwrRays[5]*IceRayTracing::c_light_ms*m;
      RxRaySolPar[2][4]=0;
    }
    delete []GetRx2ShwrRays;
  }

  double *PathAndTimeOfRays=RadioScatter::getPathAndTimeOfRays(TxRaySolPar, RxRaySolPar);    
  double *StraightLineDTimeTx=IceRayTracing::GetDirectRayPar_Cnz(Tx_z,Tx2ShwrDist,Shwr_z, n_rel);
  double *StraightLineDTimeRx=IceRayTracing::GetDirectRayPar_Cnz(Rx_z,Rx2ShwrDist,Shwr_z, n_rel);
  
  delete []PathAndTimeOfRays;
  delete []StraightLineDTimeTx;
  delete []StraightLineDTimeRx;
  
  double *output=new double[4];
  output[0]=PathAndTimeOfRays[0]-StraightLineDTimeTx[2]*s;
  output[1]=PathAndTimeOfRays[1]-StraightLineDTimeTx[2]*s;
  output[2]=PathAndTimeOfRays[2]-StraightLineDTimeRx[2]*s;
  output[3]=PathAndTimeOfRays[3]-StraightLineDTimeRx[2]*s;
  return output;
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

double RadioScatter::makeRays(TLorentzVector point, double e, double l, double e_i){
  if(RADIOSCATTER_INIT==false){
    makeTimeHist();
    RADIOSCATTER_INIT=true;
    RSCAT_HIST_DECLARE=true;

    if(USE_RAYTRACING==true){
      delta_t.resize(ntx,vector<vector<vector<double> > >(nrx,vector<vector<double> >(4,vector<double>(totalShowerPoints))));
      
      TVector3 showerStart=TVector3(event.position);
      TVector3 stretch=TVector3(event.direction);
      double showerAngle=stretch.Theta();
  
      for(int insh=0;insh<totalShowerPoints;insh++){
	showerPointDist2Vertex[insh]=insh*m;
	if(insh==0){  
	  stretch.SetMag(showerPointDist2Vertex[insh] + 0.000001*m);
	}else{
	  stretch.SetMag(showerPointDist2Vertex[insh]);
	}

	TVector3 showerPoint=showerStart+stretch;
	double pointDepth=showerPoint.Z();
	double showerAngle=stretch.Theta();
	double shwrPropTime=(showerPointDist2Vertex[insh]/c_light)*ns;
    
	for(int i=0;i<ntx;i++){
	  for(int j=0;j<nrx;j++){
	    double *rayTraceTimes=rayTrace(tx[i], rx[j], showerPoint);
	    delta_t[i][j][0][insh]=rayTraceTimes[0]+shwrPropTime;
	    delta_t[i][j][1][insh]=rayTraceTimes[1]+shwrPropTime;
	    delta_t[i][j][2][insh]=rayTraceTimes[2]+shwrPropTime;
	    delta_t[i][j][3][insh]=rayTraceTimes[3]+shwrPropTime;
	    delete []rayTraceTimes;
	  }
	}
      }  
      
      spline.resize(ntx,vector<vector<gsl_spline* > >(nrx,vector<gsl_spline* >(4)));
      for(int i=0;i<ntx;i++){
	for(int j=0;j<nrx;j++){    
	  for (int irxtx = 0; irxtx < 4; irxtx++){
	    spline[i][j][irxtx]= gsl_spline_alloc (gsl_interp_cspline, totalShowerPoints);
	    gsl_spline_init (spline[i][j][irxtx], showerPointDist2Vertex, delta_t[i][j][irxtx].data(), totalShowerPoints);
	  }
	}
      }
      
    }///if raytracing is true
  }

  double rx_time, rx_amplitude, rx_phase, point_time, t_step=0.;

  if(NPRIMARIES_SET==0){
      if(TARGET_ENERGY_SET==1){
	setNPrimaries(event.targetEnergy/event.primaryEnergy);
	//cout<<"<<<<<<<<<<<<<<<"<<event.targetEnergy<<" "<<event.nPrimaries<<" "<<event.primaryEnergy<<endl<<">>>>>>>>>>>>>"<<endl;
      }
    }

  if(SCALE_BY_ENERGY==1&&ENERGY_SCALING_SET==0){

    if(PRIMARY_ENERGY_SET==0){
      cout<<"The primary energy has not been set, but scaling has been enabled. Please tell radioscatter what the energy of the primary cascade particle is. If running within GEANT4, do this in RunAction() initialization. If not running in GEANT4, set directly with setPrimaryEnergy(), with the energy provided in MeV."<<endl;
      exit(0);
    }
    
      
    scaleByEnergy();
  }
  
  if(ENERGY_SCALING_SET==1){
    auto pVec=point.Vect()-event.position;
    pVec.SetMag(pVec.Mag()*zscale);
    auto pointNew=pVec+event.position;
    point.SetXYZT(pointNew.X(), pointNew.Y(), pointNew.Z(), point.T()*tscale);
  }

  //double zz=point.Z()*zscale;
  //double tt=point.T()*tscale;

  //std::cout<<point.Z()<<" ";
  //point.SetZ(zz);
  //point.SetVect(point.Vect()*zScaleVec);
  //  std::cout<<zscale<<" "<<point.Z()<<std::endl;
  //point.SetT(tt);
  
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
      nu_col = sqrt(kB*(300)*kelvin/m_e)*collisionalCrossSection*(N_ice);

      event.totNScatterers+=n;//track total number of scatterers

      //the full scattering amplitude pre-factor  
      //double prefactor = -rxEffectiveHeight*n*n_primaries*e_radius*omega/(pow(omega, 2)+pow(nu_col, 2));
      double prefactor = -rxFactor*n*n_primaries*e_radius*omega/(pow(omega, 2)+pow(nu_col, 2));

      //x position of charge w/r/t shower axis
      TVector3 vec=(point.Vect()-event.position);
      double x = vec.Mag()*abs(sin(vec.Unit().Angle(event.direction)));
      
      //only consider charges within 1 m (10 moliere radii, to accelerate sim). uses mm. so does c_light in alpha below, which is unitless
      double x_0=(abs(x)>1000?0:(1000-abs(x)));
      //plasma frequency
      double omega_p=sqrt(plasma_const*n_e)*1e-9;//in ns^-1

      //the screening term. as derived in paper
      double alpha= ((omega_p*omega_p)/(2.*c_light_r))*(nu_col/(omega*omega + nu_col*nu_col))*x_0;

      double attn_factor = exp(-alpha);

      prefactor=prefactor*attn_factor;
      
      TLorentzVector point_temp=point;      
      //are we calculating in a region where there is a boundary? (like in a test-beam setup
      if(BOUNDARY_FLAG==1){
	TVector3  q1 = findRefractionPlane(tx[i], point);//make a plane where the refraction will happen
	TVector3 j1;
	j1.SetZ(findRefractionPointZ(q1.X(), q1.Z(), tx_interface_dist[i]));//find the refraction point in this plane on interface 
	TVector3  q2 = findRefractionPlane(rx[j], point);
	TVector3  j2;
	j2.SetZ(findRefractionPointZ(q2.X(), q2.Z(), rx_interface_dist[j]));
	j1.SetX(tx_interface_dist[i]);
	j1.SetY(0);
	j2.SetX(rx_interface_dist[j]);
	j2.SetY(0);
	TVector3 l1;
	l1.SetXYZ(q1.X()-tx_interface_dist[i], 0., q1.Z()-j1.Z());
	TVector3 l2;
	l2.SetXYZ(q2.X()-rx_interface_dist[j], 0., q2.Z()-j2.Z());

	if(USE_RAYTRACING==false){
	  point_time=point_temp.T();
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
	    point_temp.SetT(point_time);
	  }///while loop
	  
	}else{///add raytracing times since its on
	  for (int irxtx = 0; irxtx < 4; irxtx++){
	    double showerPoint=vec.Mag();
	    double dt_time=gsl_spline_eval(spline[i][j][irxtx], showerPoint, acc);
	  
	    point_time=point_temp.T()+dt_time;
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
	      point_temp.SetT(point_time);
	    }///while loop
	  }///irxtx loop
	}///add raytracing else statement
	
      }
      //assuming transmitter and interaction and receiver are in a medium with same refractive index.    
      else{
	if(USE_RAYTRACING==false){	  
	  point_time=point_temp.T();
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
	    point_temp.SetT(point_time);
	  }///while loop

	}else{///add raytracing times since raytracing is on
	  for (int irxtx = 0; irxtx < 4; irxtx++){
	    double showerPoint=vec.Mag();
	    double dt_time=gsl_spline_eval(spline[i][j][irxtx], vec.Mag(), acc);
	    
	    point_time=point_temp.T()+dt_time;    	
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
	      point_temp.SetT(point_time);

	    }///while loop
	  }///irxtx for loop
	}///add raytracing else statement
	
      }
    }
    return 1;
  }
}



//calculate screeing effects. at the moment, cannot be run in 'realtime' inside of GEANT4.
double RadioScatter::doScreening(TTree * tree, int entry){
  TLorentzVector  point(0,0,0,0);
  TLorentzVector * test=0;
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
    //    std::cout<<i<<std::endl;
    double tof=(test->Vect()-point.Vect()).Mag()/c_light;
    double t_a0=test->T()+tof;
    double t_a1=t_a0+lifetime;
    double t_b0=point.T();
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
//   TLorentzVector  point(0,0,0,0);
//   TLorentzVector * test=0;
//   double n_e;
//   tree->SetBranchAddress("point", &test);
//   tree->SetBranchAddress("n_e", &n_e);
//   tree->GetEntry(entry);
//   //TH3F *hist = new TH3F("asdf", "asdf", 40, 1, -1, 40, 1, -1, 40, 1, -1);
//   std::vector<double> xx, yy, zz;
//   point=*test;//assign "point" to this electron
//   int entries = tree->GetEntries();
//   double E_eff=0;
//     for(int i=0;i<entries;i++){
//     tree->GetEntry(i);
//     //    std::cout<<i<<std::endl;
//     double tof=(test->Vect()-point.Vect()).Mag()/c_light;
//     double t_a0=test->T()+tof;
//     double t_a1=t_a0+lifetime;
//     double t_b0=point.T();
//     double t_b1=t_b0+lifetime;
//     if(t_a0>t_b1||t_a1<t_b0)continue;
//     xx.push_back(test->X());
//     yy.push_back(test->y());
//     zz.push_back(test->z());
//     double amp = getAmplitudeFromAt(tx_voltage, *test, point);
//     double phase = getPhaseFromAt(*test, point);
//     double ee = -(e_radius*amp*cos(phase)*n_e*n_primaries);//negative from polarity flip
//     //	if(ee!=0&&!isnan(ee)&&!isinf(ee))E_eff+=ee;
//     E_eff+=ee;
//     //	std::cout<<i-entry<<std::endl;
//     std::cout<<"amp: "<<amp<<" phase: "<<phase<<" E: "<<E_eff<<" "<<ee<<std::endl;
//     //std::cout<<point.X()<<" "<<test->X()<<" "<<std::endl;

// 	//    std::cout<<E_eff<<std::endl;
//   }
//     std::cout<<"ssdf"<<xx.size()<<std::endl;
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

 TVector3 RadioScatter::findRefractionPlane(TLorentzVector p1, TLorentzVector p2){
  TVector3 k;
  k.SetX(abs(p2.X()-p1.X()));

  k.SetY(abs(p2.Y()-p1.Y()));
  k.SetZ(abs(p2.Z()-p1.Z()));
  // std::cout<<std::endl<<k.X()<<" "<<k.Y()<<" "<<k.Z()<<std::endl;
  // std::cout<<p1.X()<<" "<<p1.Y()<<" "<<p1.Z()<<std::endl;
  // std::cout<<tx[index].X()<<" "<<tx[index].Y()<<" "<<tx[index].Z()<<std::endl;
  // std::cout<<p2.X()<<" "<<p2.Y()<<" "<<p2.Z()<<std::endl<<std::endl;
  double theta = atan(k.Y()/k.Z());
  k.RotateX(theta);
  return k;
}
/*
find the path length including refraction
 */
   double RadioScatter::findPathLengthWithRefraction(TLorentzVector p1, TLorentzVector p2, double interface_dist_x){
  TVector3 k = findRefractionPlane(p1,p2);
  double jz = findRefractionPointZ(k.X(), k.Z(), interface_dist_x);
  TVector3 j(interface_dist_x, 0., jz);//from source point to refraction point
  TVector3 l(k.X()-interface_dist_x, 0., k.Z()-jz);//from refraction point to interaction point
  double mag = j.Mag()+l.Mag();
  //  std::cout<<mag-k.Mag()<<std::endl;
  return mag;
  
}

std::vector<std::vector<TH1F*>> RadioScatter::scaleHist(float num_events=1.){
  std::vector<std::vector<TH1F*>> outhist(ntx, std::vector<TH1F*>(nrx));
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
  std::cout<<std::setprecision(10);
  std::cout<<"total number of scatterers: "<<  event.totNScatterers<<std::endl;
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
  std::vector<double> xvals, yvals;
  
  f->cd();
  
  
  t->Fill();
  if(debug==1){
    std::cout<<"The RadioScatter root file: "<<std::endl<<f->GetName()<<std::endl<<" has been filled. This is run number "<<fRunCounter<<"."<<std::endl;
  }
  
  std::cout<<"Run total N scatterers:"<<event.totNScatterers<<std::endl; 
  
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
    std::cout<<"The RadioScatter root file: "<<std::endl<<f->GetName()<<std::endl<<" has been filled. This is run number "<<fRunCounter<<"."<<std::endl;
    
  }
  
  std::cout<<"Event total N scatterers:"<<event.totNScatterers<<std::endl; 

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
    rss->polarization=rs->polarization;
    rss->nPrimaries=rs->nPrimaries;
    rss->primaryParticleEnergy=rs->primaryParticleEnergy();
    rss->inelasticity=rs->inelasticity;
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
	rss->delta[j][k]=rs->delta[j][k];
	rss->beta[j][k]=rs->beta[j][k];
	rss->angTVP[j][k]=rs->angTVP[j][k];
	rss->angRVP[j][k]=rs->angRVP[j][k];
	rss->angTVR[j][k]=rs->angTVR[j][k];
	rss->doppler[j][k]=2*rs->freq*cos(rs->delta[j][k])*cos(rs->beta[j][k]/2);
      }
    }
    outfile->cd();
    outtree->Fill();
  }
  outfile->Write();
  std::cout<<"the summary file:"<<std::endl<<outfile->GetName()<<std::endl<<"has been written."<<std::endl;
  outfile->Close();
  return 0;
}

    
  void RadioScatter::close(){
    TFile *f=    ((TFile *)(gROOT->GetFile(output_file_name)));
    TString fname = f->GetName();

    f->Write();

    std::cout<<"The RadioScatter root file: "<<std::endl<<fname<<std::endl<<"has been written."<<std::endl;
    f->Close();

    f=TFile::Open(fname);
    if(MAKE_SUMMARY_FILE==1){
      makeSummary(f);
    }
    f->Close();

  }





#include <RadioScatter/RadioScatter.hh>


void doIt(double lifetimens, double frequency, double power){
  TString frequencyStr=TString::Itoa(frequency, 10);
  TString powerStr=TString::Itoa(power, 10);
  TString lifetimeStr=TString::Itoa(lifetimens, 10);
  //the radioscatter object.
  RadioScatter *radio = new RadioScatter();

  //mandatory call before doing anything else. make the output ROOT file. it isn't opened while the simulation is running, so you can do other ROOT things below, but it gets re-opened and filled at the end. 
  radio->makeOutputFile("../doc/raytrace_test_"+frequencyStr+"MHz_"+powerStr+"W_"+lifetimeStr+"ns.root");

  //here we set the radioscatter simulation parameters.  

  //set the transmitter position  
  radio->setTxPos(-10.*m, 10*m, -10*m);
  //set the number of receivers
  radio->setNRx(3);
  
  //can set the receiver positions like this by indexing them
  // radio->setRxPos(10*m,10*m,10*m, 0);//10m x 10m x 10m
  //radio->setRxPos(50*m,10*m,10*m, 1);//50m x 10m x 10m

  //or like this, and the index of receiver increments automatically up to nrx-1
  auto r0=Hep3Vector(10*m, 0*m, -10.*m);
  auto r1=Hep3Vector(20*m, 0*m, -10.*m);
  auto r2=Hep3Vector(30*m, 0*m, -10.*m);

  radio->setRxPos(r0);//sets the 0th receiver as r0
  radio->setRxPos(r1);//sets the 1st receiver as r1
  radio->setRxPos(r2);//sets the 1st receiver as r2

  double nPrimaries=1e9;
  radio->setNPrimaries(nPrimaries);//the number of 'primaries'. used to imitate the density of a charge bunch or the approximate density of a higher-energy primary
  radio->setTxFreq(frequency);//transmitter frequency (set by user)

  radio->setTxPower(power); //transmitter power in Watts (set by user) 
  radio->setTxGain(9);//transmitter gain in dB
  radio->setRxSampleRate(4.);//receiver sample rate
  radio->setRxGain(9);//receiver gain in dB
  radio->setRecordWindowLength(1000);//length of the received window
  radio->setCalculateUsingAttnLength(1);//use ice attenuation length in the calculation?
  radio->setPolarization("horizontal");//antenna polarization. currently vertical = (0,0,1) and horizontal = (0,1,0); 
  radio->setPrimaryEnergy(1e4);//10GeV in MeV. need to set this for the scaling to be correct, if you simulate a higher energy shower than the input file (which was 10GeV)
  radio->setScaleByEnergy(0);//scales the shower longitudinally by a factor to simulate a higher energy shower. not exact, but fast.
  radio->setMakeSummary(1);//make a nice summary file for simple plotting of things like peak power, voltage, etc.
  radio->setPlasmaLifetime(lifetimens);//set the plasma lifetime



  HepLorentzVector pt; //the point of the ionization from which we calculate the individual scatter.
  double ionizationE=.00001;//MeV, the mean ionization energy.
  int num=10;
  double dz=10*m/num;
  for(int i=0;i<num;i++){
    pt.setX(0.);
    pt.setY(0.);
    pt.setZ(-(double)i*dz);
    pt.setT(pt.z()/c_light);
    //    cout<<i<<endl<<endl<<endl<<endl;
    cout<<i<<endl;
    //calculate the scatter from this ionization deposit.
    radio->makeRaysRayTrace(pt, 1, 1, ionizationE);
  }

  //don't forget to write the run!
  radio->writeRun();

  //close the radioscatter root file. not really mandatory.  
  radio->close();
    
}

int main(int argc, char**argv){
  if(argc!=4){
    cout<<"usage ./exampleRayTrace <plasma lifetime [ns]> <frequency [MHz]> <tx power [W]>"<<endl;
    exit (0);
  }
  doIt(stod(argv[1]), stod(argv[2]), stod(argv[3]));
}


#include <RadioScatter/RadioScatter.hh>
/*
This example shows the basic use of radioscatter. the input file is a cascade produced by GEANT4 in which we have saved the 4-vector (x, y, z, t) of each step in each particle's track, as well as the step length and the energy deposited in that step. This particular shower is a 10 GeV shower which we scale to 10PeV.

we then set the simulation parameters. transmitter frequency, power, polarization, antenna gains, the number of primaries (which allows for the simulation of a 'bunch'), whether or not to scale the shower longitudinally by energy (so that a GeV shower can be simulated and scaled instead of a PeV shower which would take forever for GEANT4 to produce) etc. for more information, see the documentation. 

finally, we simply loop through the file and calculate the scatter. for this, compile against the root libraries (`root-config --cflags --glibs --libs`) and the radioscatter library (-lRadioScatter) once you've installed. it should run in about 1 second. try:

./scatterFromCascade 3 500 100

for a 500MHz transmitter at 100W with a 3ns plasma lifetime.

 */

void doIt(double lifetimens, double frequency, double power){
  TFile *ff = TFile::Open("../doc/shower_particleinfo_10gev_single.root");
  TTree * tree =(TTree*)ff->Get("tracks");
  TString frequencyStr=TString::Itoa(frequency, 10);
  TString powerStr=TString::Itoa(power, 10);
  TString lifetimeStr=TString::Itoa(lifetimens, 10);
  double edep=0, steplength=0, x=0, y=0, z=0, t=0;
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("t", &t);
  tree->SetBranchAddress("edep", &edep);
  tree->SetBranchAddress("steplength", &steplength);

  //the radioscatter object.
  RadioScatter *radio = new RadioScatter();

  //mandatory call before doing anything else. make the output ROOT file. it isn't opened while the simulation is running, so you can do other ROOT things below, but it gets re-opened and filled at the end. 
  radio->makeOutputFile("../doc/output_test_"+frequencyStr+"MHz_"+powerStr+"W_"+lifetimeStr+"ns.root");

  //here we set the radioscatter simulation parameters.  

  //set the transmitter position  
  radio->setTxPos(-10.*m, 10*m, 0);
  //set the number of receivers
  radio->setNRx(3);
  
  //can set the receiver positions like this by indexing them
  // radio->setRxPos(10*m,10*m,10*m, 0);//10m x 10m x 10m
  //radio->setRxPos(50*m,10*m,10*m, 1);//50m x 10m x 10m

  //or like this, and the index of receiver increments automatically up to nrx-1
  auto r0=TVector3(10*m, 0*m, 0.);
  auto r1=TVector3(20*m, 0*m, 0.);
  auto r2=TVector3(30*m, 0*m, 0.);

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
  radio->setRecordWindowLength(100);//length of the received window
  radio->setCalculateUsingAttnLength(1);//use ice attenuation length in the calculation?
  radio->setPolarization("horizontal");//antenna polarization. currently vertical = (0,0,1) and horizontal = (0,1,0); 
  radio->setPrimaryEnergy(1e4);//10GeV in MeV. need to set this for the scaling to be correct, if you simulate a higher energy shower than the input file (which was 10GeV)
  radio->setScaleByEnergy(0);//scales the shower longitudinally by a factor to simulate a higher energy shower. not exact, but fast.
  radio->setMakeSummary(1);//make a nice summary file for simple plotting of things like peak power, voltage, etc.
  radio->setPlasmaLifetime(lifetimens);//set the plasma lifetime


  int entries=tree->GetEntries();
  TLorentzVector pt; //the point of the ionization from which we calculate the individual scatter.
  double ionizationE=.000059;//MeV, the mean ionization energy.
  for(int i=0;i<entries;i++){
    tree->GetEntry(i);
    pt.SetX(x);
    pt.SetY(y);
    pt.SetZ(z);
    pt.SetT(t);
    //calculate the scatter from this ionization deposit.
    radio->makeRays(pt, edep, steplength, ionizationE);
  }

  //don't forget to write the run!
  radio->writeRun();

  //close the radioscatter root file. not really mandatory.  
  radio->close();
    
}

int main(int argc, char**argv){
  if(argc!=4){
    cout<<"usage ./scatterFromCascade <plasma lifetime [ns]> <frequency [MHz]> <tx power [W]>"<<endl;
    exit (0);
  }
  doIt(stod(argv[1]), stod(argv[2]), stod(argv[3]));
}


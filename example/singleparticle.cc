/*
this simple example shows scattering from a single particle. 
it's a test of the module, to be sure scaling is correct, as a sanity check. 
the tx and rx are set at 1 m from the vertex. the output voltage is 1 volt (set in mV).
therefore, the received voltage should be equal in magnitude to the classical electron radius, 
roughly 2.8x10^-12 mV.

the output file will have a TTree called 'tree' with a single branch, 
which is a RadioScatterEvent object. 
 */

#include <RadioScatter/RadioScatter.hh>


int doscatt(TString filename, int lifetimens){
  //this makes a new radio scatter object
  RadioScatter *radio = new RadioScatter();
  radio->makeOutputFile(filename);

  //set the transmitter and receiver positions (in mm)
  Hep3Vector tx(1000,0,0);
  Hep3Vector rx(1000,0,0);

  //make a default 4 vector for the charge
  HepLorentzVector point(0,0,0,0);

  //set the parameters of the simulation
  radio->setRxPos(rx);
  radio->setTxPos(tx);
  radio->setTxVoltage(1000);
  radio->setTxFreq(1400);
  radio->setRxSampleRate(2.4);
  radio->setPlasmaLifetime(lifetimens);
  radio->setPolarization("horizontal");
  radio->setRecordWindowLength(100);

  //calculate the scattering from this single particle.
  //if you wanted to do more complicated stuff, make a collection of points in 4 space and loop them all!
  radio->makeRays(point, 1., 1., 1.);

  //write the 'event'
  radio->writeEvent();
  //close the file.
  radio->close();
  
  return 1;
}

int main(int argc, char**argv){
  if(argc!=3){
    std::cout<<"usage: ./singleparticle filename lifetime(ns)"<<std::endl;
    exit(0);
  }
  doscatt(argv[1], stoi(argv[2]));
}

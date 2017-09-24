#include "RSmessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "globals.hh"

RSmessenger::RSmessenger(RadioScatter *rscat)
  :rs(rscat)
{
  RSDir = new G4UIdirectory("/RS/");
  RSDir->SetGuidance("directory for radio scatter commands");

  freqCommand = new G4UIcmdWithADouble("/RS/setFreq", this);
  freqCommand->SetGuidance("set freq!");
  freqCommand->SetParameterName("choice", false);

  sampleRateCommand = new G4UIcmdWithADouble("/RS/setSampleRate", this);
  sampleRateCommand->SetGuidance("set sampleRate!");
  sampleRateCommand->SetParameterName("choice", false);

  //  freqCommand->AvailableForStates(G4State_PreInit,G4State_Idle);

  voltageCommand = new G4UIcmdWithADouble("/RS/setVoltage", this);
  voltageCommand->SetGuidance("set voltage!");
  voltageCommand->SetParameterName("choice", false);

  polarizationCommand = new G4UIcmdWithAString("/RS/setPolarization", this);
  polarizationCommand->SetGuidance("set Polarization!");
  polarizationCommand->SetParameterName("choice", false);

  nPrimariesCommand = new G4UIcmdWithADouble("/RS/setNprimaries", this);
  nPrimariesCommand->SetGuidance("set Nprimaries!");
  nPrimariesCommand->SetParameterName("choice", false);

  showCWCommand = new G4UIcmdWithADouble("/RS/setShowCWFlag", this);
  showCWCommand->SetGuidance("set ShowCW!");
  showCWCommand->SetParameterName("choice", false);
  
  setTxOnCommand = new G4UIcmdWithADouble("/RS/setTxOnTime", this);
  setTxOnCommand->SetGuidance("set tx on time!");
  setTxOnCommand->SetParameterName("choice", false);

  txpositionCmd = new G4UIcmdWith3VectorAndUnit("/RS/setTxPos",this);
  txpositionCmd->SetGuidance("Set position of the tx.");
  txpositionCmd->SetParameterName("x","y","z",true,true);
  txpositionCmd->SetDefaultUnit("m");
  txpositionCmd->SetUnitCandidates("cm m km");

  rxpositionCmd = new G4UIcmdWith3VectorAndUnit("/RS/setRxPos",this);
  rxpositionCmd->SetGuidance("Set position of the rx.");
  rxpositionCmd->SetParameterName("x","y","z",true,true);
  rxpositionCmd->SetDefaultUnit("m");
  rxpositionCmd->SetUnitCandidates("cm m km");
  
  
  setTxOffCommand = new G4UIcmdWithADouble("/RS/setTxOffTime", this);
  setTxOffCommand->SetGuidance("set tx off time!");
  setTxOffCommand->SetParameterName("choice", false);

  indexOfRefractionCommand = new G4UIcmdWithADouble("/RS/setIndexOfRefraction", this);
  indexOfRefractionCommand->SetGuidance("set indexOfRefraction!");
  indexOfRefractionCommand->SetParameterName("choice", false);

}


RSmessenger::~RSmessenger()
{
  delete freqCommand;
  delete voltageCommand;
  delete nPrimariesCommand;
  delete polarizationCommand;
  delete indexOfRefractionCommand;
  delete showCWCommand;
}

void RSmessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command==polarizationCommand)rs->setPolarization((char*)newValue.c_str());
  
  double val = (double)StoD(newValue);
  if(command==freqCommand)rs->setTxFreq(val);
  if(command==voltageCommand)rs->setTxVoltage(val);

  if(command==indexOfRefractionCommand)rs->setRelativeIndexOfRefraction(val);
  if(command==sampleRateCommand)rs->setSampleRate(val);
  if(command==nPrimariesCommand)rs->setNprimaries(val);
  if(command==showCWCommand)rs->setShowCWFlag(val);
  if(command==setTxOnCommand)rs->setTxOnTime(val);
  if(command==setTxOffCommand)rs->setTxOffTime(val);
  if(command==txpositionCmd)rs->setTxPos(txpositionCmd->GetNew3VectorValue(newValue));
  if(command==rxpositionCmd)rs->setRxPos(rxpositionCmd->GetNew3VectorValue(newValue));
  

}
  

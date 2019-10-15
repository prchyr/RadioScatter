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

  freqCommand = new G4UIcmdWithADouble("/RS/setTxFreq", this);
  freqCommand->SetGuidance("set freq!");
  freqCommand->SetParameterName("choice", false);

  sampleRateCommand = new G4UIcmdWithADouble("/RS/setRxSampleRate", this);
  sampleRateCommand->SetGuidance("set sampleRate!");
  sampleRateCommand->SetParameterName("choice", false);

  receiverGainCommand = new G4UIcmdWithADouble("/RS/setReceiverGain", this);
  receiverGainCommand->SetGuidance("set receiverGain!");
  receiverGainCommand->SetParameterName("choice", false);

  transmitterGainCommand = new G4UIcmdWithADouble("/RS/setTransmitterGain", this);
  transmitterGainCommand->SetGuidance("set transmitterGain!");
  transmitterGainCommand->SetParameterName("choice", false);

  windowLengthCommand = new G4UIcmdWithADouble("/RS/setRecordWindowLength", this);
  windowLengthCommand->SetGuidance("set the length of receiver record window");
  windowLengthCommand->SetParameterName("choice", false);

  setNTxCommand = new G4UIcmdWithADouble("/RS/setNTx", this);
  setNTxCommand->SetGuidance("set the number of transmitters");
  setNTxCommand->SetParameterName("choice", false);

  setNRxCommand = new G4UIcmdWithADouble("/RS/setNRx", this);
  setNRxCommand->SetGuidance("set the number of receivers");
  setNRxCommand->SetParameterName("choice", false);
  //  freqCommand->AvailableForStates(G4State_PreInit,G4State_Idle);

  voltageCommand = new G4UIcmdWithADouble("/RS/setTxVoltage", this);
  voltageCommand->SetGuidance("set voltage!");
  voltageCommand->SetParameterName("choice", false);

  powerCommand = new G4UIcmdWithADouble("/RS/setTxPower", this);
  powerCommand->SetGuidance("set power!");
  powerCommand->SetParameterName("choice", false);

  lifetimeCommand = new G4UIcmdWithADouble("/RS/setPlasmaLifetime", this);
  lifetimeCommand->SetGuidance("set lifetime!");
  lifetimeCommand->SetParameterName("choice", false);

  polarizationCommand = new G4UIcmdWithAString("/RS/setPolarization", this);
  polarizationCommand->SetGuidance("set Polarization!");
  polarizationCommand->SetParameterName("choice", false);

    setParticleInfoFilenameCommand = new G4UIcmdWithAString("/RS/setParticleInfoFilename", this);
  setParticleInfoFilenameCommand->SetGuidance("set SetParticleInfoFilename!");
  setParticleInfoFilenameCommand->SetParameterName("choice", false);


  
  nPrimariesCommand = new G4UIcmdWithADouble("/RS/setNPrimaries", this);
  nPrimariesCommand->SetGuidance("set NPrimaries!");
  nPrimariesCommand->SetParameterName("choice", false);

  makeSummaryCommand = new G4UIcmdWithADouble("/RS/setMakeSummary", this);
  makeSummaryCommand->SetGuidance("set makeSummary!");
  makeSummaryCommand->SetParameterName("choice", false);
  
  showCWCommand = new G4UIcmdWithADouble("/RS/setShowCWFlag", this);
  showCWCommand->SetGuidance("set ShowCW!");
  showCWCommand->SetParameterName("choice", false);

  setFillByEventCommand = new G4UIcmdWithADouble("/RS/setFillByEvent", this);
  setFillByEventCommand->SetGuidance("calculate by event not by run!");
  setFillByEventCommand->SetParameterName("choice", false);

  setScaleByEnergyCommand = new G4UIcmdWithADouble("/RS/setScaleByEnergy", this);
  setScaleByEnergyCommand->SetGuidance("set scale by energy!");
  setScaleByEnergyCommand->SetParameterName("choice", false);

  
  setFillParticleInfoCommand = new G4UIcmdWithADouble("/RS/setFillParticleInfo", this);
  setFillParticleInfoCommand->SetGuidance("fill the particle info tuples");
  setFillParticleInfoCommand->SetParameterName("choice", false);

  setCalculateUsingAttnLengthCommand = new G4UIcmdWithADouble("/RS/setCalculateUsingAttnLength", this);
  setCalculateUsingAttnLengthCommand->SetGuidance("set use attn length!");
  setCalculateUsingAttnLengthCommand->SetParameterName("choice", false);

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
  delete setNTxCommand;
  delete setNRxCommand;
  delete lifetimeCommand;
  delete receiverGainCommand;
  delete windowLengthCommand;
  delete setFillByEventCommand;
  delete makeSummaryCommand;
  delete setScaleByEnergyCommand;
  delete setFillParticleInfoCommand;
  delete setParticleInfoFilenameCommand;
  delete voltageCommand;
  delete powerCommand;
  delete nPrimariesCommand;
  delete polarizationCommand;
  delete indexOfRefractionCommand;
  delete showCWCommand;
}

void RSmessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command==polarizationCommand)rs->setPolarization((char*)newValue.c_str());
  if(command==setParticleInfoFilenameCommand)rs->setParticleInfoFilename((char*)newValue.c_str());
  
  double val = (double)StoD(newValue);
  if(command==freqCommand)rs->setTxFreq(val);
  if(command==voltageCommand)rs->setTxVoltage(val);
  if(command==powerCommand)rs->setTxPower(val);
  if(command==setNTxCommand)rs->setNTx(val);
  if(command==setNRxCommand)rs->setNRx(val);
  if(command==receiverGainCommand)rs->setReceiverGain(val);
  if(command==transmitterGainCommand)rs->setTransmitterGain(val);
  if(command==makeSummaryCommand)rs->setMakeSummary(val);
  if(command==setScaleByEnergyCommand)rs->setScaleByEnergy(val);
  if(command==lifetimeCommand)rs->setPlasmaLifetime(val);
  if(command==setFillByEventCommand)rs->setFillByEvent(val);
  if(command==setFillParticleInfoCommand)rs->setFillParticleInfo(val);
  if(command==windowLengthCommand)rs->setRecordWindowLength(val);
  if(command==indexOfRefractionCommand)rs->setRelativeIndexOfRefraction(val);
  if(command==sampleRateCommand)rs->setRxSampleRate(val);
  if(command==setCalculateUsingAttnLengthCommand)rs->setCalculateUsingAttnLength(val);
  if(command==nPrimariesCommand)rs->setNPrimaries(val);
  if(command==showCWCommand)rs->setShowCWFlag(val);
  if(command==setTxOnCommand)rs->setTxOnTime(val);
  if(command==setTxOffCommand)rs->setTxOffTime(val);
  if(command==txpositionCmd)rs->setTxPos(txpositionCmd->GetNew3VectorValue(newValue));
  if(command==rxpositionCmd)rs->setRxPos(rxpositionCmd->GetNew3VectorValue(newValue));
  

}
  

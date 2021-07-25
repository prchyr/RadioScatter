#ifndef RSmessenger_flag
#define RSmessenger_flag 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "RadioScatter/RadioScatter.hh"

class RadioScatter;
class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class RSmessenger: public G4UImessenger
{
public:
  RSmessenger(RadioScatter*);
  ~RSmessenger();
  void SetNewValue(G4UIcommand*, G4String);

private:
  RadioScatter* rs;
  G4UIdirectory*             RSDir;
  G4UIcmdWithADouble*        freqCommand;
  G4UIcmdWithADouble*        setNTxCommand;
  G4UIcmdWithADouble*        setNRxCommand;
  G4UIcmdWithADouble*        voltageCommand;
  G4UIcmdWithADouble*        makeSummaryCommand;
  G4UIcmdWithADouble*        receiverGainCommand;
  G4UIcmdWithADouble*        rxGainCommand;
  G4UIcmdWithADouble*        transmitterGainCommand;
  G4UIcmdWithADouble*        txGainCommand;
  G4UIcmdWithADouble*        powerCommand;
  G4UIcmdWithADoubleAndUnit*        setPrimaryEnergyCommand;
  G4UIcmdWithADoubleAndUnit*        setTargetEnergyCommand;
  G4UIcmdWithADouble*        setScaleByEnergyCommand;
  G4UIcmdWithADouble*        setUseRayTracingCommand;
  G4UIcmdWithADouble*        setCrossSectionCommand;
  G4UIcmdWithADouble*        lifetimeCommand;
  G4UIcmdWithADouble*        windowLengthCommand;
  G4UIcmdWith3Vector*        polarizationCommand;
  G4UIcmdWithADouble*        setFillByEventCommand;
  G4UIcmdWithADouble*        setFillParticleInfoCommand;
  G4UIcmdWithAString*        setParticleInfoFilenameCommand;
  G4UIcmdWithADouble*        nPrimariesCommand;
  G4UIcmdWithADouble*        sampleRateCommand;
  G4UIcmdWithADouble*        indexOfRefractionCommand;
  G4UIcmdWithADouble*        setCalculateUsingAttnLengthCommand;
  G4UIcmdWithADouble*        showCWCommand;
  G4UIcmdWithADouble*        setTxOnCommand;
  G4UIcmdWithADouble*        setTxOffCommand;
  G4UIcmdWith3VectorAndUnit*        txpositionCmd;
  G4UIcmdWith3VectorAndUnit*        rxpositionCmd;
  G4UIcmdWithAString*        setOutputDirectory;
  G4UIcmdWithAString*        setOutputFileName;
  
};
#endif

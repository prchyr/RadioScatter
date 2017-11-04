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
  G4UIcmdWithADouble*        voltageCommand;
  G4UIcmdWithADouble*        windowLengthCommand;
  G4UIcmdWithAString*        polarizationCommand;
  G4UIcmdWithADouble*        nPrimariesCommand;
  G4UIcmdWithADouble*        sampleRateCommand;
  G4UIcmdWithADouble*        indexOfRefractionCommand;
  G4UIcmdWithADouble*        showCWCommand;
  G4UIcmdWithADouble*        setTxOnCommand;
  G4UIcmdWithADouble*        setTxOffCommand;
  G4UIcmdWith3VectorAndUnit*        txpositionCmd;
  G4UIcmdWith3VectorAndUnit*        rxpositionCmd;
  G4UIcmdWithAString*        setOutputDirectory;
  G4UIcmdWithAString*        setOutputFileName;
  
};
#endif

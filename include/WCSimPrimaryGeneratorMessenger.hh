#ifndef WCSimPrimaryGeneratorMessenger_h
#define WCSimPrimaryGeneratorMessenger_h 1

#include "G4UIcmdWithADouble.hh"

class WCSimPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool; //jl145

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4Tokenizer.hh"

class WCSimPrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* mpga);
  ~WCSimPrimaryGeneratorMessenger();
  
 public:
  void     SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
  
 private:
  WCSimPrimaryGeneratorAction* myAction;
  
 private: //commands
  G4UIdirectory*      mydetDirectory;
  G4UIcmdWithAString* genCmd;
  G4UIcmdWithAString* fileNameCmd;
  G4UIcmdWithAString* isotopeCmd;
  G4UIcmdWithAString* radonCmd;
  G4UIcmdWithADouble* radioactive_time_window_Cmd;
  G4UIcmdWithAString* volgenCmd;
  G4UIcmdWithAString* fileNameCmdCosmics;
  
  G4UIcmdWithABool* StorePhotons;

  void IsotopeCommand(G4String newValue);
  void RadonCommand(G4String newValue);
  void VolumeGenCommand(G4String newValue);
  
};

#endif



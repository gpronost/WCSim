#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "G4UIcmdWithABool.hh" //jl145

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{
  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  //T. Akiri: Addition of laser
  //G. Pronost: Addition of radioactive element (from F. Nova code)
  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, radioactive, radon, volume");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  //T. Akiri: Addition of laser
  //G. Pronost: Addition of radioactive element (from F. Nova code)
  genCmd->SetCandidates("muline gun laser gps radioactive volume cosmics radon");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  fileNameCmdCosmics = new G4UIcmdWithAString("/mygen/cosmicsfile",this);
  fileNameCmdCosmics->SetGuidance("Select the file of cosmics.");
  fileNameCmdCosmics->SetGuidance(" Enter the file name of the cosmics file");
  fileNameCmdCosmics->SetParameterName("fileName",true);
  fileNameCmdCosmics->SetDefaultValue("inputvectorfile");
  
  //G. Pronost (from F. Nova code) Radioactivity:

  radioactive_time_window_Cmd = new G4UIcmdWithADouble("/mygen/radioactive_time_window",this);
  radioactive_time_window_Cmd->SetGuidance("Select time window for radioactivity");
  radioactive_time_window_Cmd->SetParameterName("radioactive_time_window",true);
  radioactive_time_window_Cmd->SetDefaultValue(0.);

  isotopeCmd = new G4UIcmdWithAString("/mygen/isotope",this);
  isotopeCmd->SetGuidance("Select properties of radioactive isotope");
  isotopeCmd->SetGuidance("[usage] /mygen/isotope ISOTOPE LOCATION ACTIVITY ACTIVITY_MIN ");
  isotopeCmd->SetGuidance("     ISOTOPE : Tl208, Bi214, K40");
  isotopeCmd->SetGuidance("     LOCATION : water PMT");
  isotopeCmd->SetGuidance("     ACTIVITY : (double) activity of isotope (Bq) in case of water: Activity at the border of the ID in Bq/m3");
  isotopeCmd->SetGuidance("     ACTIVITY_MIN : (double) activity of isotope (Bq) in case of water: Activity at the center of the ID in Bq/m3");
  isotopeCmd->SetGuidance("     FORCE : (int) force single decay (0: normal, 1: force)");
  G4UIparameter* param;
  param = new G4UIparameter("ISOTOPE",'s',true);
  param->SetDefaultValue("Tl208");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("LOCATION",'s',true);
  param->SetDefaultValue("water");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("ACTIVITY",'d',true);
  param->SetDefaultValue("0");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("ACTIVITY_MIN",'d',true);
  param->SetDefaultValue("0");
  isotopeCmd->SetParameter(param);
  param = new G4UIparameter("FORCE",'i',true);
  param->SetDefaultValue("0");
  isotopeCmd->SetParameter(param);
  
  radonScalingCmd = new G4UIcmdWithAString("/mygen/radon_scaling",this);
  radonScalingCmd->SetGuidance("Select scalling scenario");
  radonScalingCmd->SetGuidance("[usage] /mygen/radon SCENARIO ");
  radonScalingCmd->SetGuidance("     SCENARIO : A, B, C");
  param = new G4UIparameter("SCENARIO",'s',true);
  param->SetDefaultValue("C");
  radonScalingCmd->SetParameter(param);
  
  //G. Pronost: Addition of Volume generator
  volgenCmd = new G4UIcmdWithAString("/mygen/volume",this);
  volgenCmd->SetGuidance("Simulate particle randomly in a given volume or on a given surface");
  volgenCmd->SetGuidance("[usage] /mygen/volume TYPE PARTICLE ENERGY LOCATION WALLDIST");
  volgenCmd->SetGuidance("     TYPE     : VOLUME SURFACE");
  volgenCmd->SetGuidance("     PARTICLE : PDG Id of the particle");
  volgenCmd->SetGuidance("     ENERGY   : (double) energy of the particle (MeV)");
  volgenCmd->SetGuidance("     LOCATION : ID OD");
  volgenCmd->SetGuidance("     WALLDIST : (double) minimum (maximum if neg) distance from wall (mm)");
  //G4UIparameter* param;
  param = new G4UIparameter("TYPE",'s',false);
  param->SetDefaultValue("VOLUME");
  volgenCmd->SetParameter(param);
  param = new G4UIparameter("PARTICLE",'s',false);
  param->SetDefaultValue("22");
  volgenCmd->SetParameter(param);
  param = new G4UIparameter("ENERGY",'s',false);
  param->SetDefaultValue("4.0");
  volgenCmd->SetParameter(param);
  param = new G4UIparameter("LOCATION",'s',false);
  param->SetDefaultValue("ID");
  volgenCmd->SetParameter(param);
  param = new G4UIparameter("WALLDIST",'s',true);
  param->SetDefaultValue("0.0");
  volgenCmd->SetParameter(param);

}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
  
  delete radioactive_time_window_Cmd;
  delete isotopeCmd;
  delete radonScalingCmd;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    if (newValue == "muline")
    {
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "radioactive") //G. Pronost: Addition of Radioactivity (from F. Nova code)
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(true);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "radon" ) //G. Pronost: Addition of Radioactivity (from F. Nova code) with SK distribution in the water
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(true);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "volume") //G. Pronost: Addition of Volume generator
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(true);
      myAction->SetCosmicsGenerator(false);
    }
    else if ( newValue == "cosmics")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
      myAction->SetRadioactiveEvtGenerator(false);
      myAction->SetRadonEvtGenerator(false);
      myAction->SetVolumeEvtGenerator(false);
      myAction->SetCosmicsGenerator(true);
    }
  }

  if( command == fileNameCmd || command == fileNameCmdCosmics )
  {
    myAction->OpenVectorFile(newValue);
    G4cout << "Input vector file set to " << newValue << G4endl;
  }
  
  //G. Pronost: Addition of Radioactivity (from F. Nova code)
  if( command==isotopeCmd )
  {
    IsotopeCommand(newValue);
  }
  
  //G. Pronost: Addition of Radon generator
  if( command==radonScalingCmd )
  {
    RadonScalingCommand(newValue);
  }

  if( command==radioactive_time_window_Cmd )
  {
    myAction->SetRadioactiveTimeWindow(StoD(newValue));
  }

  //G. Pronost: Addition of Volume generator
  if( command==volgenCmd )
  {
    VolumeGenCommand(newValue);
  }
}

void WCSimPrimaryGeneratorMessenger::IsotopeCommand(G4String newValue)
{

   G4Tokenizer next( newValue );

   G4String isotope = next();
   G4String location = next();
   G4double activity     = StoD(next());
   G4double activity_min = StoD(next());
   G4int    force        = StoI(next());

   myAction->AddRadioactiveSource(isotope, location, activity, activity_min, force);


}

void WCSimPrimaryGeneratorMessenger::RadonScalingCommand(G4String newValue)
{

   G4Tokenizer next( newValue );

  G4String scenario = next();
  G4int iScenario = 0;
   
  if ( scenario == "A" ) iScenario = 1; // Relative scaling with respect to full ID volume (Pessimistic)
  if ( scenario == "B" ) iScenario = 2; // Relative scaling with respect to fiducial volume
  if ( scenario == "C" ) iScenario = 3; // Absolute scaling with respect to ID border (Optimistic)
   
   myAction->SetRadonScenario(iScenario);


}

void WCSimPrimaryGeneratorMessenger::VolumeGenCommand(G4String newValue)
{

   G4Tokenizer next( newValue );

   G4String vol_type  = next();
   G4int    vol_pdgid = StoI(next());
   G4double vol_e     = StoD(next());
   G4String vol_loc   = next();
   G4double vol_wall  = StoD(next());
  

   myAction->SetVolumeSource(vol_type, vol_pdgid, vol_e, vol_loc, vol_wall);
   
   G4cout << " Type: " << vol_type << " PDGID " << vol_pdgid << " E " << vol_e << " LOC " << vol_loc << " WALL " << vol_wall << G4endl;
}



G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingGunEvtGenerator())
      { cv = "gun"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = "gps"; }
    else if(myAction->IsUsingRadioactiveEvtGenerator())
      { cv = "radioactive"; } //G. Pronost: Addition of Radioactivity (from F. Nova code)
    else if(myAction->IsUsingRadonEvtGenerator())
      { cv = "radon"; } //G. Pronost: Addition of Radioactivity (from F. Nova code)
    else if(myAction->IsUsingCosmicsGenerator())
      { cv = "cosmics"; }
  }
  
  return cv;
}


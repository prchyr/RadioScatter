//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4GeneralParticleSourceData.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "RSmessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//    std::vector<G4double> x, y, z, t;  
RunAction::RunAction(RadioScatter *radio, DetectorConstruction * dc)
 : G4UserRunAction()
{

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     
  fRadio=radio;
  fDetConstruction = dc;
  fRSmessenger = new RSmessenger(fRadio);
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");

  //  if(fRadio->FILL_PARTICLE_INFO==1){
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

    // Book histograms, ntuple
    //
  
    // Creating histograms
    //analysisManager->CreateH1("pulse", "received radio pulse", 250*10, 550., 800., "ns");
    // analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
    // analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
    // analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
    // analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

    // Creating ntuple
    //
    analysisManager->CreateNtuple("event_tots", "Edep and TrackL");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Egap");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleDColumn("Lgap");
    analysisManager->CreateNtupleDColumn("Eice");
    analysisManager->CreateNtupleDColumn("Lice");
    analysisManager->CreateNtupleDColumn("Eps");
    analysisManager->CreateNtupleDColumn("Lps");
    analysisManager->CreateNtupleDColumn("etot");
    analysisManager->CreateNtupleIColumn("neplus");
    analysisManager->CreateNtupleIColumn("neminus");
    analysisManager->CreateNtupleDColumn("inite");
    analysisManager->FinishNtuple();

    //  EventAction *eventAction;
    analysisManager->CreateNtuple("tracks", "particle tracks and deposited energy");
    analysisManager->CreateNtupleIColumn("trackid");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleDColumn("t");//global time
    analysisManager->CreateNtupleSColumn("pname");//p name
    analysisManager->CreateNtupleSColumn("ptype");
    analysisManager->CreateNtupleDColumn("charge");
    analysisManager->CreateNtupleDColumn("edep");
    analysisManager->CreateNtupleDColumn("enonionizing");
    analysisManager->CreateNtupleSColumn("lv");//logical volume of step
    analysisManager->CreateNtupleDColumn("steplength");
    analysisManager->CreateNtupleDColumn("ke");
    analysisManager->CreateNtupleDColumn("tote");
    analysisManager->CreateNtupleDColumn("p");
    analysisManager->CreateNtupleDColumn("inite");
    analysisManager->CreateNtupleDColumn("pt");
    analysisManager->CreateNtupleDColumn("eta");
    analysisManager->FinishNtuple();

    analysisManager->CreateNtuple("stats", "particle stats, numbers etc.");
    analysisManager->CreateNtupleIColumn("trackid");
    analysisManager->CreateNtupleSColumn("ptype");
    analysisManager->CreateNtupleSColumn("pname");
    analysisManager->CreateNtupleDColumn("charge");
    analysisManager->CreateNtupleDColumn("etot");
    analysisManager->CreateNtupleSColumn("lv");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleDColumn("t");
    analysisManager->CreateNtupleDColumn("einit");
    analysisManager->CreateNtupleDColumn("tracklength");
    // analysisManager->CreateNtupleDColumn("evno");
    analysisManager->FinishNtuple();
    //  }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  //delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* r)
{


  
  /*
*************
set the RadioScatter parameters. 
    
these can be set here, or more simply (and efficiently) in the macro file using the RSmessenger. 
that way, they can be changed at run-time, without the need to re-compile. 
set them here if you want some default values between runs, without needing to re-define in the macro, if you want. 
    
there is one exception: the tx/interacion medium/rx distances, which must come from the GEANT4 system,
must be sent to radioscatter if the refraction machinery is to work. 

if these distances are not provided, the simulation assumes that the shower is happening in free space, with no refraction calculation.

set the transmitter frequency, output power, and amplification

*******

some example commands below, all of which can (and should) be set in the macro file.

  */
  //fRadio->setTxVoltage(200.);
  //fRadio->setTxFreq(1400.);
  /*//    fRadio->setRxVals(20., 1.);//samplerate (ns^-1), gain (not db, 100=40db)
    //fRadio->setTxPos(4.5*m, 0.*m, -1.*m);
    //fRadio->setRxPos(-4.3*m, 0.*m, -1.*m);
    //fRadio->makeTimeHist();
    // fRadio->setTxOnTime(0);
    //fRadio->setTxOffTime(30);
    
    */

  //THESE CALLS ARE MANDATORY FOR THE REFRACTION CALCULATIONS
  //if not set, it is assumed everything happens in free space.
  //they tell RadioScatter the distance to the refraction boundary
  //  fRadio->setTxInterfaceDistX(fDetConstruction->GetTgtPV()->GetLogicalVolume()->GetSolid()->DistanceToIn(fRadio->getTxPos()));
  //fRadio->setRxInterfaceDistX(fDetConstruction->GetTgtPV()->GetLogicalVolume()->GetSolid()->DistanceToIn(fRadio->getRxPos()));

    auto gpsDat=G4GeneralParticleSourceData::Instance();
    auto gps = gpsDat->GetCurrentSource();
    fRadio->setPrimaryEnergy(gps->GetParticleEnergy());
    TVector3 dirV(gps->GetParticleMomentumDirection().x(), gps->GetParticleMomentumDirection().y(),gps->GetParticleMomentumDirection().z());
    fRadio->setPrimaryDirection(dirV);
    TVector3 pos(gps->GetParticlePosition().x(), gps->GetParticlePosition().y(),gps->GetParticlePosition().z());
    
    fRadio->setPrimaryPosition(pos);


    //fRadio->setPrimaryDirection(gps->GetParticleMomentumDirection());
    //fRadio->setPrimaryPosition(gps->GetParticlePosition());
  
  if(fRadio->FILL_PARTICLE_INFO==1){
    //redundant root files saved through geant root system
    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    // Open an output file
    G4String runno = G4UIcommand::ConvertToString(r->GetRunID());
    //G4String dir="$HOME/doc/root/";
    //    G4String fileName = "slac_rf_"+runno+"_";
    //    G4String fileName = "slac_rf_photon";
    //    G4String fileName = "shower_particleinfo";
    G4String fileName = (G4String)fRadio->PARTICLE_INFO_FILENAME;
    //  std::cout<<r->GetRunID()<<std::endl;
    analysisManager->OpenFile(fileName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* r)
{
  //finish the run with RadioScatter()  

  //optional second arguments scales the histogram output by the number of events, so that output amplitude is correct
  //setting the filename here gives acces to G4 variables for the filename, as below


  G4String runno = G4UIcommand::ConvertToString(r->GetRunID());
  G4String freq = G4UIcommand::ConvertToString(fRadio->getFreq());

  //   auto gpsDat=G4GeneralParticleSourceData::Instance();
  // auto gps = gpsDat->GetCurrentSource();
  // fRadio->setPrimaryEnergy(gps->GetParticleEnergy());
  // fRadio->setPrimaryDirection(gps->GetParticleMomentumDirection());
  // fRadio->setPrimaryPosition(gps->GetParticlePosition());
  //be sure to send writeRun the number of events in the run so that it scales the output histogram properly!!! assumes you have opened the root file already in the main program
  if(fRadio->FILL_BY_EVENT==0){
    auto gpsDat=G4GeneralParticleSourceData::Instance();
    auto gps = gpsDat->GetCurrentSource();
    fRadio->setPrimaryEnergy(gps->GetParticleEnergy());

    TVector3 dirV(gps->GetParticleMomentumDirection().x(), gps->GetParticleMomentumDirection().y(),gps->GetParticleMomentumDirection().z());
    fRadio->setPrimaryDirection(dirV);

    TVector3 pos(gps->GetParticlePosition().x(), gps->GetParticlePosition().y(),gps->GetParticlePosition().z());
    
    fRadio->setPrimaryPosition(pos);

    //    fRadio->setPrimaryDirection(gps->GetParticleMomentumDirection());
    //fRadio->setPrimaryPosition(gps->GetParticlePosition());
    fRadio->writeRun((float)r->GetNumberOfEvent());
    //    fRadio->writeRun(1);
  }
  //  fRadio->writeEvent("$HOME/Documents/physics/geant/root/slac_rf_rs_"+runno+"_freq_"+freq+".root", (float)r->GetNumberOfEvent());

  //optionally, set the filename in the macro by /RS/setOutputFileName and then here:
  
  
  //  fRadio->printEventStats(); 
  

  //redundant filling of histogram through geant
  if(fRadio->FILL_PARTICLE_INFO==1){
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// $Id: SteppingAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"
#include "RunAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"
//#include "CLHEP/Units/PhysicalConstants.h"
//#include "CLHEP/Vector/LorentzVector.h"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include <vector> 


 

G4double x, y, z, t, v, e, gt,tote, tracklength, trackid, inite, edep, advanced_t, step_length, E_i;
double rx_amplitude;
G4String pname, ptype, lv;
HepLorentzVector p;
G4int charge;
G4Track * track;
//TH1F * hist;


auto analysisManager = G4AnalysisManager::Instance();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int lastTrackID, thisTrackID;

SteppingAction::SteppingAction(const DetectorConstruction* detectorConstruction,EventAction* eventAction, RadioScatter * radio)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    fRadio(radio)
{
  
  // double interface_dist = fDetConstruction->GetTgtPV()->GetLogicalVolume()->GetSolid()->DistanceToIn(fRadio->getTxPos());
  // G4ThreeVector pos = fRadio->getTxPos();
     // cout<<"***************"<<endl<<endl<<interface_dist<<endl<<endl<<"********"<<endl;
 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void SteppingAction::UserSteppingAction(const G4Step* step)
{

// auto gpsdata=G4GeneralParticleSourceData::Instance();  
// G4ThreeVector pos = gpsdata->GetCurrentSource()->GetParticlePosition();
//  cout<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<endl;

  
  track = step->GetTrack();
  
  // energy deposit
  edep = step->GetTotalEnergyDeposit();
  v = track->GetVelocity();
  //position, time
  p.setX(track->GetPosition().x());
  p.setY(track->GetPosition().y());
  p.setZ(track->GetPosition().z());
  p.setT(track->GetGlobalTime());
  t = track->GetGlobalTime();
  step_length = step->GetStepLength();
  E_i = track->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();


  /*
make the rays. 

will fill a histogram with default parameters if makeHist was called in run action. 
otherwise, will output the resultant E_field, and can be put into a G4 histogram. 

or both. 
   */

  rx_amplitude = fRadio->makeRays(p, edep, step_length, E_i);  
  //  fRadio->makeRaysReal(p, edep, step_length, E_i);
  //fRadio->makeRaysImag(p, edep, step_length, E_i);

  if(fRadio->FILL_PARTICLE_INFO==1){
    //auto manager = G4RunManager::GetRunManager();
    //auto run = manager->GetCurrentRun();
    
    //if (run->GetNumberOfEvent()<10){
      fillParticleInfoTuples(step);
      // }
  }
  // advanced_t = fRadio->getRxTime(p);
  //analysisManager->FillH1(0, advanced_t, rx_amplitude);

  //if(  track->GetTrackID()==10){


  //lastTrackID=thisTrackID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::fillParticleInfoTuples(const G4Step *step){
  track = step->GetTrack();
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  edep = step->GetTotalEnergyDeposit();
  thisTrackID=step->GetTrack()->GetTrackID();
  auto pos = step->GetPreStepPoint()->GetPosition();

  //for each step
  if(volume == fDetConstruction->GetTgtPV())
    {

    //once per track, fill with last track's values. 
    if(thisTrackID != lastTrackID && thisTrackID!=1){
      fEventAction->fillStatVector(
				   trackid,
				   ptype,
				   pname,
				   charge,
				   e,
				   lv,
				   x,
				   y,
				   z,
				   gt,
				   inite,
				   tracklength
				   );
      //identify the charge discrepancy for shower particles
      if(track->GetParticleDefinition()->GetParticleName()=="e+"&&track->GetTotalEnergy()>70.)fEventAction->IncrementEplus();
      if(track->GetParticleDefinition()->GetParticleName()=="e-"&&track->GetTotalEnergy()>70.)fEventAction->IncrementEminus();
    

    }
    //for every step, record all the things. 
    fEventAction->fillTrajectoryVector(track->GetTrackID(),
				       track->GetPosition().x(),
				       track->GetPosition().y(),
				       track->GetPosition().z(),
				       // pos.x(),
				       // pos.y(),
				       // pos.z(),
				       track->GetGlobalTime(),
				       track->GetParticleDefinition()->GetParticleName(),
				       track->GetParticleDefinition()->GetParticleType(),
				       track->GetParticleDefinition()->GetPDGCharge(),
				       //step->GetTotalEnergyDeposit(),
				       edep,
				       step->GetNonIonizingEnergyDeposit(),
				       track->GetVolume()->GetLogicalVolume()->GetName(),
				       step->GetStepLength(),
				       track->GetKineticEnergy(),
				       track->GetTotalEnergy(),
				       track->GetMomentum().mag(),
				       track->GetVertexKineticEnergy(),
				       track->GetMomentum().perp(),
				       track->GetMomentum().pseudoRapidity()
				       );
    //accumulate total energy
    fEventAction->doTotE(edep);

    //  analysisManager->FillH1(0, track->GetPosition().z(), edep);
    //set parameters for full track. to be filled at the last step of the track.
    // trackid = track->GetTrackID();
    // ptype = track->GetParticleDefinition()->GetParticleType();
    // pname = track->GetParticleDefinition()->GetParticleName();
    // charge = track->GetParticleDefinition()->GetPDGCharge();
    // lv = track->GetVolume()->GetLogicalVolume()->GetName();
    // e += edep;
    // x=      pos.x();
    // y=      pos.y();
    // z=      pos.z();
    // gt=      track->GetGlobalTime();
    // inite = track->GetVertexKineticEnergy();
    // tote=      track->GetTotalEnergy();
    // tracklength = track->GetTrackLength();
   

    //set initial particle energy.
    if(thisTrackID ==1)fEventAction->setInitE(track->GetVertexKineticEnergy());
    }   
}

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
// $Id: EventAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4VTrajectory.hh"
#include "G4UIsession.hh"
#include "G4VisExecutive.hh"
#include "G4GeneralParticleSourceData.hh"
#include "Randomize.hh"
#include <iomanip>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RadioScatter *radio)
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   fEnergyPs(0.),
   fEnergyIce(0.),
   fTrackLPs(0.),
   fTrackLIce(0.),
   fEnergyTot(0.)
{
  fRadio = radio;
  //G4cout<<track->GetTrackID()<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fEnergyPs=0.;
  fEnergyIce=0.;
  fTrackLPs=0.;
  fTrackLIce=0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  fEnergyTot =0.;
  
  if(fRadio->FILL_BY_EVENT==1){
    //set some parameters of the shower for radioscatter
    auto gpsDat=G4GeneralParticleSourceData::Instance();
    auto gps = gpsDat->GetCurrentSource();
    fRadio->setPrimaryEnergy(gps->GetParticleEnergy());
    TVector3 dir(gps->GetParticleMomentumDirection().x(), gps->GetParticleMomentumDirection().y(),gps->GetParticleMomentumDirection().z());
    fRadio->setPrimaryDirection(dir);
    TVector3 pos(gps->GetParticlePosition().x(), gps->GetParticlePosition().y(),gps->GetParticlePosition().z());
    
    fRadio->setPrimaryPosition(pos);
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
// get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  // analysisManager->FillH1(0, fEnergyAbs);
  // analysisManager->FillH1(1, fEnergyGap);
  // analysisManager->FillH1(2, fTrackLAbs);
  // analysisManager->FillH1(3, fTrackLGap);
  //  analysisManager->FillH1(4, fEnergyIce);

  //fill tup
  // analysisManager->FillNtupleDColumn(0, fEnergyAbs);
  // analysisManager->FillNtupleDColumn(1, fEnergyGap);
  // analysisManager->FillNtupleDColumn(2, fTrackLAbs);
  // analysisManager->FillNtupleDColumn(3, fTrackLGap);
  // analysisManager->FillNtupleDColumn(4, fEnergyIce);
  // analysisManager->FillNtupleDColumn(5, fTrackLIce);
  // analysisManager->FillNtupleDColumn(6, fEnergyPs);
  // analysisManager->FillNtupleDColumn(7, fTrackLPs);
  // analysisManager->FillNtupleDColumn(8, fEnergyTot);
  // analysisManager->FillNtupleIColumn(9, fTotEplus);
  // analysisManager->FillNtupleIColumn(10, fTotEminus);
  // analysisManager->FillNtupleDColumn(11, fInitE);
  // analysisManager->AddNtupleRow();  

  if (fRadio->FILL_PARTICLE_INFO==1){
  // //for each step in each track in this event,
  // //record the coordinates, time, deposited energy, and logical volume

  int size = stepEdep.size();
  for(int i=0;i<size;i++){
    analysisManager->FillNtupleIColumn(1, 0, trackID[i]);
    analysisManager->FillNtupleDColumn(1, 1, trackXVec[i]);
    analysisManager->FillNtupleDColumn(1, 2, trackYVec[i]);
    analysisManager->FillNtupleDColumn(1, 3, trackZVec[i]);
    analysisManager->FillNtupleDColumn(1, 4, trackTVec[i]);
    analysisManager->FillNtupleSColumn(1, 5, trackPNameVec[i]);
    analysisManager->FillNtupleSColumn(1, 6, trackPTypeVec[i]);   
    analysisManager->FillNtupleDColumn(1, 7, trackCVec[i]);
    analysisManager->FillNtupleDColumn(1, 8, stepEdep[i]);
    analysisManager->FillNtupleDColumn(1, 9, stepEion[i]);
    analysisManager->FillNtupleSColumn(1, 10, stepLV[i]);   
    analysisManager->FillNtupleDColumn(1, 11, stepLength[i]);
    analysisManager->FillNtupleDColumn(1, 12, stepKE[i]);
    analysisManager->FillNtupleDColumn(1, 13, stepTotE[i]);
    analysisManager->FillNtupleDColumn(1, 14, trackPVec[i]);
    analysisManager->FillNtupleDColumn(1, 15, trackInitE[i]);
    analysisManager->FillNtupleDColumn(1, 16, trackPTVec[i]);
    analysisManager->FillNtupleDColumn(1, 17, trackEtaVec[i]);
    analysisManager->AddNtupleRow(1);
  }

  int size2 = idVec.size();
  for(int i=0;i<size2;i++){
    analysisManager->FillNtupleIColumn(2, 0, idVec[i]);
    analysisManager->FillNtupleSColumn(2, 1, ptypeVec[i]);
    analysisManager->FillNtupleSColumn(2, 2, pnameVec[i]);
    analysisManager->FillNtupleDColumn(2, 3, chargeVec[i]);
    analysisManager->FillNtupleDColumn(2, 4, etotVec[i]);
    analysisManager->FillNtupleSColumn(2, 5, lvVec[i]);
    analysisManager->FillNtupleDColumn(2, 6, xVec[i]);
    analysisManager->FillNtupleDColumn(2, 7, yVec[i]);
    analysisManager->FillNtupleDColumn(2, 8, zVec[i]);
    analysisManager->FillNtupleDColumn(2, 9, tVec[i]);
    analysisManager->FillNtupleDColumn(2, 10, eVec[i]);
    analysisManager->FillNtupleDColumn(2, 11, trackLengthVec[i]);
    analysisManager->AddNtupleRow(2);
  }
  //
  clearTrajectoryVector();  
  clearStatVector();
  }
  if(fRadio->FILL_BY_EVENT==1){
    //    cout<<"here"<<endl;
    fRadio->writeEvent();
  }
  
 }  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

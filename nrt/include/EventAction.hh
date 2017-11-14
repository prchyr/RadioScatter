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
// $Id: EventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "Analysis.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4VTrajectory.hh"
#include <vector>
/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
  void fillTrajectoryVector(G4int id, G4double x, G4double y, G4double z, G4double t, G4String pname, G4String ptype, G4double charge, G4double etot, G4double eion, G4String lv, G4double steplength, G4double ke, G4double tote, G4double p, G4double inite, G4double pt, G4double eta);
  void fillStatVector(G4int id, G4String ptype, G4String pname, G4double charge, G4double etot, G4String lv, G4double x, G4double y, G4double z, G4double t, G4double e, G4double tracklength);
  void clearTrajectoryVector();
  void clearStatVector();
  void AddAbs(G4double de, G4double dl);
  void AddGap(G4double de, G4double dl);
  void AddPs(G4double de, G4double dl);
  void AddIce(G4double de, G4double dl);
  void doTotE(G4double e);
  void IncrementEplus();
  void IncrementEminus();
  void setInitE(G4double e);
  std::vector<G4double>& x(){return trackXVec;}
  std::vector<G4double>& y(){return trackYVec;}
  std::vector<G4double>& z(){return trackZVec;}
  std::vector<G4double>& t(){return trackTVec;}

private:
  G4double  fEnergyAbs;
  G4double  fEnergyGap;
  G4double  fTrackLAbs; 
  G4double  fTrackLGap;
  G4double fEnergyPs;
  G4double fEnergyIce;
  G4double fTrackLPs;
  G4double fTrackLIce;
  G4double fEnergyTot;
  G4int fTotEplus;
  G4int fTotEminus;
  G4int steptot=0;
  G4double  fInitE=0.;

  std::vector<G4double> trackXVec, trackYVec, trackZVec, trackTVec, trackCVec, stepEdep, stepEion, stepLength, stepKE, stepTotE, trackPVec, trackInitE, trackPTVec, trackEtaVec;
  std::vector<G4String>stepLV, trackPNameVec, trackPTypeVec;
  std::vector<G4int>trackID;

  std::vector<G4double> chargeVec, etotVec, xVec, yVec, zVec, tVec, eVec, trackLengthVec;
  std::vector<G4int> idVec;
  std::vector<G4String> ptypeVec, pnameVec, lvVec;
};

// inline functions
inline void EventAction::fillTrajectoryVector(G4int id, G4double x, G4double y, G4double z, G4double t, G4String pname, G4String ptype, G4double charge, G4double etot, G4double eion, G4String lv, G4double steplength, G4double ke, G4double tote, G4double p, G4double inite, G4double pt, G4double eta){
  trackID.push_back(id);
  trackXVec.push_back(x);
  trackYVec.push_back(y);
  trackZVec.push_back(z);
  trackTVec.push_back(t);
  trackPNameVec.push_back(pname);
  trackPTypeVec.push_back(ptype);
  trackCVec.push_back(charge);
  stepEdep.push_back(etot);
  stepEion.push_back(eion);
  stepLV.push_back(lv);
  stepLength.push_back(steplength);
  stepKE.push_back(ke);
  stepTotE.push_back(tote);
  trackPVec.push_back(p);
  trackInitE.push_back(inite);
  trackPTVec.push_back(pt);
  trackEtaVec.push_back(eta);
  //std::cout<<t<<std::endl;
}

inline void EventAction::fillStatVector(G4int id, G4String ptype, G4String pname, G4double charge, G4double etot, G4String lv, G4double x, G4double y, G4double z, G4double t, G4double e, G4double tracklength){
  idVec.push_back(id);
  ptypeVec.push_back(ptype);
  pnameVec.push_back(pname);
  chargeVec.push_back(charge);
  etotVec.push_back(etot);
  lvVec.push_back(lv);
  xVec.push_back(x);
  yVec.push_back(y);
  zVec.push_back(z);
  tVec.push_back(t);
  eVec.push_back(e);
  trackLengthVec.push_back(tracklength);
}

inline void EventAction::clearTrajectoryVector(){
  trackID.clear();
  trackTVec.clear();
  trackXVec.clear();
  trackYVec.clear();
  trackZVec.clear();
  trackTVec.clear();
  trackPNameVec.clear();
  trackPTypeVec.clear();
  trackCVec.clear();
  stepEdep.clear();
  stepEion.clear();
  stepLV.clear();
  stepLength.clear();
}

inline void EventAction::clearStatVector(){
  idVec.clear();
  ptypeVec.clear();
  pnameVec.clear();
  chargeVec.clear();
  etotVec.clear();
  lvVec.clear();
  xVec.clear();
  yVec.clear();
  zVec.clear();
  tVec.clear();
}

inline void EventAction::setInitE(G4double e){
  fInitE= e;
}
inline void EventAction::AddAbs(G4double de, G4double dl) {
  fEnergyAbs += de; 
  fTrackLAbs += dl;
}

inline void EventAction::AddGap(G4double de, G4double dl) {
  fEnergyGap += de; 
  fTrackLGap += dl;
}
inline void EventAction::AddIce(G4double de, G4double dl) {
  fEnergyIce += de; 
  fTrackLIce += dl;
}
inline void EventAction::AddPs(G4double de, G4double dl) {
  fEnergyPs += de; 
  fTrackLPs += dl;
}
inline void EventAction::doTotE(G4double e){
  fEnergyTot+=e;
  steptot++;
}
inline void EventAction::IncrementEplus(){
  fTotEplus++;
}
inline void EventAction::IncrementEminus(){
  fTotEminus++;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

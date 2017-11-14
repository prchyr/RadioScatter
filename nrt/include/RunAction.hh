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
// $Id: RunAction.hh 74265 2013-10-02 14:41:20Z gcosmo $
// 
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "RadioScatter/RadioScatter.hh"
#include "RSmessenger.hh"
class G4Run;
class DetectorConstruction;

class RunAction : public G4UserRunAction
{
public:
  RunAction(RadioScatter *radio, DetectorConstruction *);
  virtual ~RunAction();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);
  // void setRSFreq(double freq);
  // void setRSVoltage(double v);
  // void setRSIndexOfRefraction(double iof);


private:
  RadioScatter * fRadio;
  RSmessenger *fRSmessenger;
  DetectorConstruction *fDetConstruction;
};

// inline setRSFreq(double freq){
//   fRadio->setTxFreq(freq);
// }
// inline setRSVoltage(double v){
//   fRadio->SetTxVoltage(v);
// }
// inline setRSIndexOfRefraction(double iof){
//   fRadio->setRelativeIndexOfRefraction(iof);
// }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


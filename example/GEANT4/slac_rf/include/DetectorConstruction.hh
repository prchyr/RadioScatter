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
// $Id: DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // get methods
    //
  const G4VPhysicalVolume* GetPsPV() const;
  const G4VPhysicalVolume* GetTankPV() const;
  const G4VPhysicalVolume* GetAbsorberPV() const;
  const G4VPhysicalVolume* GetAirPV() const;
  const G4VPhysicalVolume* GetTankCapFrontPV() const;
  const G4VPhysicalVolume* GetTankCapBackPV() const;
  const G4VPhysicalVolume* GetBossPV() const;
     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger
  G4VPhysicalVolume*   fPsPV;     // the pre-shower physical volume
  G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
  G4VPhysicalVolume*   fTankPV;      //the tank
  G4VPhysicalVolume*   fAirPV;      // the air physical volume
  G4VPhysicalVolume*   fTankCapFrontPV;
  G4VPhysicalVolume*   fTankCapBackPV;
  G4VPhysicalVolume*   fAirCapFrontPV;
  G4VPhysicalVolume*   fAirCapBackPV;
  G4VPhysicalVolume*   fBossPV;
  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const { 
  return fAbsorberPV; 
}

inline const G4VPhysicalVolume* DetectorConstruction::GetTankPV() const  { 
  return fTankPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetPsPV() const  { 
  return fPsPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetAirPV() const  { 
  return fAirPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetTankCapFrontPV() const  {
  return fTankCapFrontPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetTankCapBackPV() const  {
  return fTankCapBackPV; 
}
inline const G4VPhysicalVolume* DetectorConstruction::GetBossPV() const  {
  return fBossPV; 
}
     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


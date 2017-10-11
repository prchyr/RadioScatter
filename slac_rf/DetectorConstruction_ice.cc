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
// $Id: DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"
//#include "CLHEP/Vector/ThreeVector.h"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(0),
   fIcePV(0),
   fAirPV(0),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Hg");
  nistManager->FindOrBuildMaterial("G4_He");
  nistManager->FindOrBuildMaterial("G4_Ar");
  nistManager->FindOrBuildMaterial("G4_Sn");
  nistManager->FindOrBuildMaterial("G4_WATER");
  nistManager->FindOrBuildMaterial("G4_AIR");
  // Liquid argon material
  G4String name, symbol;
  G4int natoms, ncomponents;
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density, temperature, fractionmass;
  //some elements
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen",symbol="O" , z= 8., a);
  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen",symbol="N", z=7., a);
  a = 4.*g/mole;
  G4Element* elHe = new G4Element(name="He", symbol="He", z=2., a);
  //  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  
  density=.92*g/cm3;
  temperature=238*kelvin;
  G4Material *Ice = new  G4Material("Ice", density, ncomponents=2, kStateSolid, temperature);
  Ice->AddElement(elH, natoms=2);
  Ice->AddElement(elO, natoms=1);

 
  //  density=.001225*g/cm3;
  density=.396*g/cm3;
  G4Material *compressedAir = new G4Material("compressedAir", density, ncomponents=2, kStateGas, 273.*kelvin, 306.2*atmosphere);
  compressedAir->AddElement(elN, fractionmass=0.7);
  compressedAir->AddElement(elO, fractionmass=0.3);

  G4Material *denseHe = new G4Material("denseHe", density, ncomponents=1, kStateGas, 273*kelvin);
  denseHe->AddElement(elHe, fractionmass=1.);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int nofLayers = 1;
  G4double psThickness = 5.5*cm;//pre-shower
  G4double absoThickness = 20.*cm;
  G4double gapThickness =  1.*mm;
  G4double iceThickness =  10.3*m;
  G4double targetSizeXY  = 1*m;

    auto layerThickness = psThickness + absoThickness + gapThickness + iceThickness;
  //    auto layerThickness = iceThickness;
  //    auto layerThickness = absoThickness + gapThickness;
  auto targetThickness = nofLayers * layerThickness;
  auto worldSizeXY = 3. * targetThickness;
  auto worldSizeZ  = 3. * targetThickness; 
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("G4_AIR");
  auto iceMaterial = G4Material::GetMaterial("Ice");
  auto airMaterial = G4Material::GetMaterial("G4_AIR");
  auto psMaterial = G4Material::GetMaterial("G4_AIR");
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial || ! iceMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  G4RotationMatrix *rot = new G4RotationMatrix();
  rot->set(G4ThreeVector(1, 0, 0), pi/2);
  
  auto txS = new G4Tubs("tx", 0., 10*cm, 100*cm, 0, twopi);
  auto txLV = new G4LogicalVolume(txS, absorberMaterial, "tx");
  new G4PVPlacement(rot, G4ThreeVector(-5*m,0, 5*m), txLV, "tx", worldLV, false, 0, fCheckOverlaps);
  
  auto rxS = new G4Tubs("rx", 0., 10*cm, 100*cm, 0, twopi);
  auto rxLV = new G4LogicalVolume(rxS, absorberMaterial, "rx");
  new G4PVPlacement(rot, G4ThreeVector(5*m,0, 5*m), rxLV, "rx", worldLV, false, 0, fCheckOverlaps);

  
  //                               
  // Target
  //  
  auto targetS
    = new G4Box("Target",     // its name
                 targetSizeXY/2, targetSizeXY/2, targetThickness/2); // its size
                         
  auto targetLV
    = new G4LogicalVolume(
                 targetS,     // its solid
                 defaultMaterial,  // its material
                 "Target");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,0,targetThickness/2-psThickness-gapThickness),  // at (0,0,0)
                 targetLV,          // its logical volume                         
                 "Target",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 targetSizeXY/2, targetSizeXY/2, layerThickness/2); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 targetLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica

  //pre-shower plate
  auto psS 
    = new G4Box("Ps",             // its name
                 targetSizeXY/2, targetSizeXY/2, psThickness/2); // its size
                         
  auto psLV
    = new G4LogicalVolume(
                 psS,             // its solid
                 psMaterial,      // its material
                 "Ps");           // its name
                                   
  fPsPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-iceThickness/2-absoThickness/2-gapThickness/2), // its position
                 psLV,            // its logical volume
                 "Ps",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  
  //air gap
  auto airS 
    = new G4Box("Air",             // its name
                 targetSizeXY/2, targetSizeXY/2, gapThickness/2); // its size
                         
  auto airLV
    = new G4LogicalVolume(
                 airS,             // its solid
                 airMaterial,      // its material
                 "Air");           // its name
                                   
  fAirPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-iceThickness/2-absoThickness/2+psThickness/2), // its position
                 airLV,            // its logical volume
                 "Air",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

 
  //                               
  // ice
  //
  auto iceS 
    = new G4Box("Ice",             // its name
                 targetSizeXY/2, targetSizeXY/2, iceThickness/2); // its size
                         
  auto iceLV
    = new G4LogicalVolume(
                 iceS,             // its solid
                 iceMaterial,      // its material
                 "Ice");           // its name
                                   
  fIcePV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-absoThickness/2+psThickness/2+gapThickness/2), // its position
                 iceLV,            // its logical volume
                 "Ice",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Absorber
  //
  auto absorberS 
    = new G4Box("Abso",            // its name
                 targetSizeXY/2, targetSizeXY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Abso");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., iceThickness/2+psThickness/2+gapThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The target is " << nofLayers << " layers of: [ "
    << psThickness/mm << "mm of " << psMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName()
    << " + "
    << iceThickness/mm << "mm of " << iceMaterial->GetName()
    << " + "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName()<< " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  //      worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  // auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // simpleBoxVisAtt->SetVisibility(true);
  // absorberLV->SetVisAttributes(simpleBoxVisAtt);
  auto iceVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0, .2));
  iceVisAtt->SetVisibility(true);
  iceVisAtt->SetForceSolid(true);
  iceLV->SetVisAttributes(iceVisAtt);
  
  auto airVisAtt= new G4VisAttributes(G4Colour(0.1,0.1,1.0));
  airVisAtt->SetVisibility(true);
  airVisAtt->SetForceSolid(true);
  airLV->SetVisAttributes(airVisAtt);

  auto psVisAtt= new G4VisAttributes(G4Colour(0.9,0.8,1.0));
  psVisAtt->SetVisibility(true);
  psVisAtt->SetForceSolid(true);
  psLV->SetVisAttributes(psVisAtt);
  
  auto absoVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  absoVisAtt->SetVisibility(true);
  absoVisAtt->SetForceSolid(true);
  absoVisAtt->SetForceAuxEdgeVisible(true);
  absorberLV->SetVisAttributes(absoVisAtt);

  
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

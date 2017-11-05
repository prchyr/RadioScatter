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

#include "DetectorConstructionT510tgt.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Transform3D.hh" 

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
   fPsPV(0),
   fAbsorberPV(0),
   fTgtPV(0),
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
  nistManager->FindOrBuildMaterial("G4_KEVLAR");
  nistManager->FindOrBuildMaterial("G4_CONCRETE");
  nistManager->FindOrBuildMaterial("G4_BRASS");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE");
  nistManager->FindOrBuildMaterial("G4_LUNG_ICRP");
  nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material * cellulose = nistManager->FindOrBuildMaterial("G4_CELLULOSE_CELLOPHANE");
  G4Material * carbon = nistManager->FindOrBuildMaterial("G4_C");

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

  //ice  
  density=.92*g/cm3;
  temperature=238*kelvin;
  G4Material *Ice = new  G4Material("Ice", density, ncomponents=2, kStateSolid, temperature);
  Ice->AddElement(elH, natoms=2);
  Ice->AddElement(elO, natoms=1);

  //compressed air
  //  density = .4397*g/cm3;//for 340atm
    density=.396*g/cm3;//for 306.2atm
    G4Material *compressedAir = new G4Material("compressedAir", density, ncomponents=2, kStateGas, 273.*kelvin, 306.2*atmosphere);
    //G4Material *compressedAir = new G4Material("compressedAir", density, ncomponents=2, kStateGas, 273.*kelvin, 340.*atmosphere);
  compressedAir->AddElement(elN, fractionmass=0.7);
  compressedAir->AddElement(elO, fractionmass=0.3);
  //dense helium
  G4Material *denseHe = new G4Material("denseHe", density, ncomponents=1, kStateGas, 273*kelvin);
  denseHe->AddElement(elHe, fractionmass=1.);

  //wood
  density = .85*g/cm3;
  G4Material *wood = new G4Material("wood", density, ncomponents = 4, kStateSolid, temperature);
  wood->AddMaterial(cellulose, fractionmass=.64);
  wood->AddMaterial(carbon, fractionmass = .23);
  wood->AddElement(elO, fractionmass = .11);
  wood->AddElement(elH, fractionmass = .02);

  
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int nofLayers = 1;
  G4double psThickness = 1.8*cm;
  G4double absoThickness = 1.*m;

  //for the T510 poly target
  // G4double dx1= .5*61*cm;
  // G4double dx2= .5*61*cm;
  // G4double dy2= .5*97*cm;//Because of getting the Half of the trapezoid
  // G4double dy1= .5*25.4*cm;
  // G4double dz=  .5*396*cm;

  //for the poly cylinder
  G4double tgtRadius = 1.5*2.54*cm;//1.5 inches
  G4double tgtLength = 12*12*2.54*cm;//12 feet

  //for a full world volume (no refraction)
  G4double dz=500*m;
  G4double dy=500*m;
  G4double dx=500*m;
  G4double dx1=dx;


  //  G4double gapThickness = 3.*tankRadius;
  G4double bossLength = 8.*cm;
  G4double bossRadius = 3.*cm;

  G4double targetSizeXY = 2.*dx1;
  auto layerThickness = 2.*dz;
  //  auto layerThickness = tgtLength;
  //    auto layerThickness = iceThickness;
  //    auto layerThickness = absoThickness + gapThickness;
  auto targetThickness = (nofLayers * layerThickness);
  auto worldSizeXY = 10. * targetThickness;
  auto worldSizeZ  = 10. * targetThickness; 
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_CONCRETE");
  auto airMaterial = G4Material::GetMaterial("G4_AIR");
  auto targetMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  auto compressedAirMaterial = G4Material::GetMaterial("G4_POLYPROPYLENE");
  //  auto compressedAirMaterial = G4Material::GetMaterial("G4_CONCRETE");
  
  // auto compressedAirMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  //  auto compressedAirMaterial = G4Material::GetMaterial("wood");
  auto iceMaterial = G4Material::GetMaterial("Ice");
  auto tankMaterial = G4Material::GetMaterial("G4_AIR");
  auto psMaterial = G4Material::GetMaterial("G4_Fe");
  auto bossMaterial = G4Material::GetMaterial("G4_Fe");
  if ( ! defaultMaterial || ! absorberMaterial || ! tankMaterial ) {
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
  G4RotationMatrix *rot2 = new G4RotationMatrix();
  rot2->set(G4ThreeVector(-1, 0, 0), pi/2);
  
  auto txS = new G4Tubs("tx", 0., 10*cm, 30*cm, 0, twopi);
  auto txLV = new G4LogicalVolume(txS, absorberMaterial, "tx");
  new G4PVPlacement(rot, G4ThreeVector(-5*m,0, 2*m), txLV, "tx", worldLV, false, 0, fCheckOverlaps);
  
  auto rxS = new G4Tubs("rx", 0., 10*cm, 30*cm, 0, twopi);
  auto rxLV = new G4LogicalVolume(rxS, absorberMaterial, "rx");
  new G4PVPlacement(rot, G4ThreeVector(5*m,0, 2*m), rxLV, "rx", worldLV, false, 0, fCheckOverlaps);

  auto daqS = new G4Box("daq", 15*cm, 10*cm, 20*cm);
  auto daqLV = new G4LogicalVolume(daqS, psMaterial, "daq");
  new G4PVPlacement(0, G4ThreeVector(5*m,0, 2.5*m), daqLV, "daq", worldLV, false, 0, fCheckOverlaps);
                                
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
                 G4ThreeVector(0,0,targetThickness/2),  // at (0,0,0)
                 targetLV,          // its logical volume                         
                 "Target",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

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
                 G4ThreeVector(0., 0.,-psThickness), // its position
                 psLV,            // its logical volume
                 "Ps",            // its name
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
                 airMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 targetLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica


  //------------------------------ 
  // Target
  //------------------------------

  // Define Target Volume

  //uncomment for the full world volume
   auto box = new G4Box("box", dx, dy, dz);
  auto tgtLV = new G4LogicalVolume(box, iceMaterial, "Tgt");
  fTgtPV = new G4PVPlacement(0, G4ThreeVector(0., -.5*dy, 0.), tgtLV, "Tgt", layerLV, false, 0, fCheckOverlaps);


  //this is the T510 target
  // auto trapezoid = new G4Trd("trapezoid", dx1, dx2, dy1, dy2, dz);
  // auto box = new G4Box("box", 2.*dx1, dy2, 2.*dz);
  // auto tgtS = new G4SubtractionSolid("Tgt", trapezoid, box, G4Translate3D(0., -dy1-dy2, 0.));
  // auto tgtLV = new G4LogicalVolume(tgtS, targetMaterial, "Tgt");
  //  fTgtPV = new G4PVPlacement(0, G4ThreeVector(0., -.5*dy2, 0.), tgtLV, "Tgt", layerLV, false, 0, fCheckOverlaps);


  //for a cylinder instead of the T510 target (comment the above stuff)
  // auto tgtS = new G4Tubs("Tgt", 0., tgtRadius, tgtLength/2, 0., twopi);
  // auto tgtLV = new G4LogicalVolume(tgtS, targetMaterial, "Tgt");
  // fTgtPV = new G4PVPlacement(0, G4ThreeVector(0,0,0), tgtLV, "Tgt", layerLV, false, 0, fCheckOverlaps);
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
                 G4ThreeVector(0., 0., targetThickness+absoThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

    //                                        
  // Visualization attributes
  //
  //      worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  
  auto airVisAtt= new G4VisAttributes(G4Colour(0.1,0.1,1.0));
  airVisAtt->SetVisibility(true);
  airVisAtt->SetForceSolid(true);
  //  airLV->SetVisAttributes(airVisAtt);

  auto psVisAtt= new G4VisAttributes(G4Colour(0.9,0.8,1.0));
  psVisAtt->SetVisibility(true);
  psVisAtt->SetForceSolid(true);
  psLV->SetVisAttributes(psVisAtt);
  // airLV->SetVisAttributes(psVisAtt);
  //  daqLV->SetVisAttributes(psVisAtt);
  // tankLV->SetVisAttributes(psVisAtt);
  // tank_cap_front_LV->SetVisAttributes(psVisAtt);
  // tank_cap_back_LV->SetVisAttributes(psVisAtt);
  auto tgtVisAtt = new G4VisAttributes(G4Colour(.9, .7, .4));
  tgtVisAtt->SetVisibility(true);
  tgtVisAtt->SetForceSolid(true);
  tgtLV->SetVisAttributes(tgtVisAtt);
  
  auto absoVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  absoVisAtt->SetVisibility(true);
  absoVisAtt->SetForceSolid(true);
  absoVisAtt->SetForceAuxEdgeVisible(true);
  absorberLV->SetVisAttributes(absoVisAtt);

  auto bossVisAtt = new G4VisAttributes(G4Colour(.5, .5, 0));
  //  bossLV->SetVisAttributes(bossVisAtt);
  //  airLV->SetVisAttributes(bossVisAtt);

  
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

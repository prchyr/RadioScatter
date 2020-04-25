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
// $Id: example.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file example.cc
/// \brief Main program of the  example

//choose the correct detector construction

//#include "DetectorConstruction.hh"
#include "DetectorConstructionT510tgt.hh"
#include "ActionInitialization.hh"
//#include "physicslist.hh"
#include "SimplePhysicsList.hh"
//#include "OpNovicePhysicsList.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "NuBeam.hh"
#include <errno.h> 
#include "Randomize.hh"
#include "RadioScatter/RadioScatter.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " example [-m macro ] [-u UIsession] [-t nThreads] [-f root filename]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
  G4String filename;
  long seed=rand();
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) filename = argv[i+1];
    else if ( G4String(argv[i]) == "-s" ) seed = atoi(argv[i+1]);
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
    
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(seed);  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  auto runManager = new G4MTRunManager;
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  auto runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  //  auto physicsList = new OpNovicePhysicsList;
  auto physicsList = new SimplePhysicsList;
  runManager->SetUserInitialization(physicsList);

  //  auto ftfp = new FTFP_BERT;
  //runManager->SetUserInitialization(ftfp);
  // auto qgsp = new QGSP_BERT;
  // runManager->SetUserInitialization(qgsp);
  // auto nubeam = new NuBeam;
  // runManager->SetUserInitialization(nubeam);
  //  auto qgsp_hp = new QGSP_BERT_HP;
  //runManager->SetUserInitialization(qgsp_hp);

  //declare the RadioScatter module
  RadioScatter *radio = new RadioScatter();

  if(! filename.size() ){
    radio->makeOutputFile("$HOME/Documents/physics/geant/root/G4_RadioScatter_output.root");
  }
  else{
    radio->makeOutputFile(filename);
  }
  auto actionInitialization = new ActionInitialization(detConstruction, radio);
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  //
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
 
  
  if ( macro.size() ) {
    // batch mode
  

    //    isInteractiveSession=true;
   G4String command = "/control/execute ";
   UImanager->ApplyCommand(command+macro);

   if(macro=="test.mac"){
     //do nothing
   }

   if(macro=="run1.mac"){
     Hep3Vector offset(-600,0,2000.);
     Hep3Vector tx, rx;

     double theta=0, phi=0;
     rx.setRThetaPhi(4000,55*degree, pi);
     radio->setRxPos(rx+offset, 0);
     rx.setRThetaPhi(6000, 75*degree, pi);
     radio->setRxPos(rx+offset, 1);
     rx.setRThetaPhi(6000, 65*degree, pi);
     radio->setRxPos(rx+offset, 2);

     tx.setRThetaPhi(4000, 125*degree, pi);
     radio->setTxPos(tx+offset, 0);
     runManager->BeamOn(1);
   }

   if (macro=="t576.mac"){

     Hep3Vector tx, rx;
      double theta=0, phi=0;
      rx.setRThetaPhi(4000,115*degree, 0.);
      radio->setRxPos(rx, 0);
      rx.setRThetaPhi(4000, 125*degree, 0.);
      radio->setRxPos(rx, 1);
      rx.setRThetaPhi(4000, 85*degree, 0.);
      radio->setRxPos(rx, 2);
      
      double rr=2.;
      while(rr<7.){
	tx.setRThetaPhi(rr*1000.,  65.*degree, 0.);

	radio->setTxPos(tx,0);
	runManager->BeamOn(1);

	rr+=.25;
      }

   }

   if(macro=="sweep.mac"){
     Hep3Vector tx, rx;
     double theta=0, phi=0;
     rx.setRThetaPhi(4000,115*degree, 0.);
     radio->setRxPos(rx, 0);
     rx.setRThetaPhi(4000, 125*degree, 0.);
     radio->setRxPos(rx, 1);
     rx.setRThetaPhi(4000, 85*degree, 0.);
     radio->setRxPos(rx, 2);

     for(int i=0;i>-7;i--){
       radio->setTxPower(i*10.);
       runManager->BeamOn(1); 
     }

   }

   if(macro=="desy.mac"){
     //do nothing
   }

   // else{

   if(macro=="surf.mac"){
     Hep3Vector tx, rx;
     double theta=0, phi=0;
     Hep3Vector offset(-600, 0, 2000);
     for(int i=0;i<11;i++){
       rx.setRThetaPhi(5000,(double)(i+4)*10.*degree, pi);
       radio->setRxPos(rx+offset, i);
     }

     tx.setRThetaPhi(6000., 65.*degree, pi);
     radio->setTxPos(tx+offset, 0);
     runManager->BeamOn(10);
 }
   if(macro=="anglet576.mac"){
     Hep3Vector tx, rx;
     double theta=0, phi=0;
     for(int i=0;i<11;i++){
       rx.setRThetaPhi(6000,(double)(i+4)*10.*degree, 0.);
       radio->setRxPos(rx, i);
     }

     tx.setRThetaPhi(6000., 65.*degree, 0.);
     radio->setTxPos(tx, 0);
     runManager->BeamOn(1);
   }
   if(macro.contains("surfacesensor")){

     Hep3Vector rx(1, 1, 1);
     for(int i=0;i<100;i++){
       auto x=(-5.+((double)i*.1))*m;
       auto y=0;
       rx.setX(x);
       rx.setY(y);
       for(int j=0;j<100;j++){
         auto z=(-5.+((double)j*.1))*m;
         rx.setZ(z);
         radio->setRxPos(rx);
       }
     }
     runManager->BeamOn(1);
   }

   if(macro=="run6.mac"||macro=="run62100.mac"||macro=="run61750.mac"){
     Hep3Vector offset(-600,0,2000.);
     Hep3Vector tx, rx;

     double theta=0, phi=0;
     rx.setRThetaPhi(4000,55*degree, pi);
     radio->setRxPos(rx+offset, 0);
     rx.setRThetaPhi(6000, 75*degree, pi);
     radio->setRxPos(rx+offset, 1);
     rx.setRThetaPhi(6000, 65*degree, pi);
     radio->setRxPos(rx+offset, 2);

     tx.setRThetaPhi(4000, 125*degree, pi);
     radio->setTxPos(tx+offset, 0);
     runManager->BeamOn(10);
   }

   if(macro=="surfsensor.mac"){
     Hep3Vector offset(0.,0.,0.);
     Hep3Vector tx, rx;

     int N=60;
     double theta=0, phi=0;
     double minz=-10000.;
     double maxz=10000.;
     double dz=abs(maxz-minz)/((double)N);

     double minx=-10000.;
     double maxx=10000.;
     double dx=abs(maxx-minx)/((double)N);

     cout<<"dx: "<<dx<<endl;
     for(int i=0;i<N;i++){
       for(int j=0;j<N;j++){
	 int num=(i*N)+j;
	 //	 cout<<num<<endl;
	 rx.set(minx+((double)i*dx), 0, minz+((double)j*dz));
	 radio->setRxPos(rx+offset, num);
       }
     }

     // tx.set(-4000, 0, -2000);
     //radio->setTxPos(tx+offset, 0);
     runManager->BeamOn(1);
   }
   
   if(macro=="fanoutt576.mac"){
     Hep3Vector offset(-600,0,2000);//coordinate system is 0,0,0 at front of target. this offset is for the coordinate system used at ESA
     Hep3Vector tx, rx;
     double theta=0, phi=0;
     rx.setRThetaPhi(6000,75.*degree, pi);
     radio->setRxPos(rx+offset, 1);
     rx.setRThetaPhi(6000,145.*degree, pi);
     radio->setRxPos(rx+offset, 2);
     for(int i=0;i<9;i++){
       rx.setRThetaPhi(4000,(double)(i+5)*10.*degree, pi);
       radio->setRxPos(rx+offset, 0);
       tx.setRThetaPhi(4000., (double)(8-i+5)*10.*degree, pi);
       radio->setTxPos(tx+offset, 0);
       runManager->BeamOn(1);
     }
   }
   
  }

  else{
    UImanager->ApplyCommand("/control/execute "+macro);
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  //close the RadioScatter object, to close the root file
  radio->close();
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
	

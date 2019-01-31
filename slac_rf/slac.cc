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
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) filename = argv[i+1];
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



   if (macro=="run1.mac"){

     Hep3Vector tx, rx;
      double theta=0, phi=0;
      int num=10;
      int xmax=40000;
      int xmin=-40000;
      int xstep=(xmax-xmin)/num;

      int ymax=40000;
      int ymin=-40000;
      int ystep=(ymax-ymin)/num;

      int zmax=40000;
      int zmin=-40000;
      int zstep=(zmax-zmin)/num;

      for(int j=0;j<num;j++){
	for(int i=0;i<num;i++){
	  for(int k=0;k<num;k++){

	  // tx=radio->getTxPos();
	  // rx=radio->getRxPos();
	  // cout<<endl<<endl<<"tx position"<<tx.x()<<" "<<tx.y()<<" "<<tx.z()<<endl;
	  // cout<<"rx position"<<rx.x()<<" "<<rx.y()<<" "<<rx.z()<<endl<<endl<<endl;
	  // //	tx.setRThetaPhi(j*m, angle, 0);
	  rx.setX(xmin+(j*xstep));
	  rx.setY(ymin+(i*ystep));
	  rx.setZ(zmin+(k*zstep));
	  // phi=pi;
	  // rx.setRThetaPhi(j*m, theta, phi);
	  // rx.setZ(rx.z()+3.5*m);
	  //	radio->setTxPos(tx);
	  radio->setRxPos(rx);
	  runManager->BeamOn(10);
	  //	theta+=2*pi/num;
	  }
	}
      }
    }


   if (macro=="t576.mac"){

     Hep3Vector tx, rx;
      double theta=0, phi=0;
      rx.setRThetaPhi(4000,0,115*degree);
      radio->setRxPos(rx, 0);
      rx.setRThetaPhi(4000, 0,125*degree);
      radio->setRxPos(rx, 1);
      rx.setRThetaPhi(4000, 0,85*degree);
      radio->setRxPos(rx, 2);
      
      double rr=3.;
      while(rr<7.){
	tx.setRThetaPhi(rr*1000., 0, 1.13446);

	radio->setTxPos(tx,0);
	runManager->BeamOn(10);

	rr+=.25;
      }

   }
    // else{
   // if(macro=="anglemac.mac"){

   // }
    // }
  }
  else{
    UImanager->ApplyCommand("/control/execute init_vis.mac");
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
	

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

#include "DetectorConstruction.hh"
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
#include "TRandom3.h"
#include "RadioScatter/RadioScatter.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4GeneralParticleSource.hh"
#include "G4GeneralParticleSourceData.hh"
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
  long seed=rand();
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

   if(macro.contains("anglemacro.mac")){
     TRandom *ran = new TRandom();
     Hep3Vector rx(1., 1., 1.), tx(1., 1., 1.);
    double mag=100000;
    rx.setMag(mag);
    rx.setTheta(pi/2);
    rx.setPhi(pi);
    radio->setNTx(1);

    //    double div=25.;

    
    
    
    //    double halfstep=step/2
    //    for(double i=0;i<twopi;i+=radianstep){
    double i=0.0001;
    tx.setPhi(pi);
    tx.setTheta(pi/2.);
    int n=0;

    //make a sphere
    double div=40;
    double radianstep=twopi/div;
    radio->setNRx(div*(div/2));
    while(i<twopi){

      double j=0.0001;
      while(j<pi){
    	//for(double j=radianstep;j<pi;j+=radianstep){

    	rx.setTheta(j);
    	rx.setPhi(i);

    	rx.setMag(mag);
    	radio->setRxPos(rx, n);
    	n++;
    	cout<<rx.x()<<" "<<rx.y()<<" "<<rx.z()<<endl;
    	j+=radianstep;
    	//	cout<<rx.r()<<" "<<rx.theta()<<" "<<rx.phi()<<endl;
      }
      i+=radianstep;
    }
    
    //make a filled cube
    // double div=6.;
    //double step = mag/div;
    // radio->setNRx(div*div*div);
    // for(int i=0;i<div;i++){
    //   for(int j=0;j<div;j++){
    // 	for(int k=0;k<div;k++){
    // 	  rx.setX(i*step);
    // 	  rx.setY(j*step);
    // 	  rx.setZ(k*step);
    // 	  radio->setRxPos(rx, n);
    // 	  n++;
    // 	}
    //   }
    // }
    cout<<endl<<n<<endl<<endl;

    // double stepp=pi/100.;
    // auto gpsDat=G4GeneralParticleSourceData::Instance();
    // auto gps = gpsDat->GetCurrentSource();
    // for(int i=1;i<100;i++){
    //   //  tx.setPhi(ran->Rndm()*twopi);
    //   //tx.setTheta(ran->Rndm()*(pi/2));
    //   //      tx.setMag(10000.*i);
    //   //tx.setMag(ran->Rndm()*1000000.);
    //   tx.setPhi(0);
    //   tx.setTheta((i*stepp)-(pi/2.));
    //   tx.setMag(mag);
    //   //radio->setTxPos(tx, 0);
    //   gps->GetAngDist()->SetParticleMomentumDirection(tx);  
    //   runManager->BeamOn(1);
    //   //      cout<<tx.x()<<" "<<tx.y()<<" "<<tx.z()<<endl;
    // }
    //    runManager->BeamOn(1);
    runManager->BeamOn(100);
  }


   if (macro.contains("run1.mac")){

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

   if(macro.contains("effectivevol.mac")){
     //nothing here
     
   }

   if(macro.contains("effvol_surface")){

     }

   if(macro.contains("surface_array")){
     TRandom3 *rann=new TRandom3();
     int num=50;
     radio->setNRx(num*num);
     int xmax=200000;
     int xmin=-xmax;
     int xstep=(xmax-xmin)/num;

     int ymax=200000;
     int ymin=-xmax;
     int ystep=(ymax-ymin)/num;

     int z=0;
     Hep3Vector tx=Hep3Vector(0,0,0);
     Hep3Vector rx;
     rx.setZ(0.);
     for(int j=0;j<num;j++){
       double xx=xmin+(j*xstep);
       rx.setX(xx);
       for(int i=0;i<num;i++){
	 double yy=ymin+(i*ystep);
	 rx.setY(yy);
	 radio->setRxPos(rx);
       }
     }

     int nThrow=10;
     for(int i=0;i<nThrow;i++){
       auto logE= rann->Uniform(5., 9.);
       radio->setNPrimaries(pow(10, logE));
       runManager->BeamOn(1);
     }
   }

   
   if (macro.contains("threshmacro.mac")){

     Hep3Vector tx, rx;
     double theta=0, phi=0;
     int num=5;
     int xmax=10000000;
     int xmin=-10000000;
     int xstep=(xmax-xmin)/num;

     int ymax=10000000;
     int ymin=-10000000;
     int ystep=(ymax-ymin)/num;

     int zmax=10000000;
     int zmin=-10000000;
     int zstep=(zmax-zmin)/num;

     for(int j=0;j<num;j++){
       for(int i=0;i<num;i++){
	 for(int k=0;k<num;k++){
	   for(int l=0;l<num+4;l++){
	     // tx=radio->getTxPos();
	     // rx=radio->getRxPos();
	     // cout<<endl<<endl<<"tx position"<<tx.x()<<" "<<tx.y()<<" "<<tx.z()<<endl;
	     // cout<<"rx position"<<rx.x()<<" "<<rx.y()<<" "<<rx.z()<<endl<<endl<<endl;
	     // //	tx.setRThetaPhi(j*m, angle, 0);
	     rx.setX(1000*pow(10, j));
	     rx.setY(1000*pow(10, i));
	     rx.setZ(1000*pow(10, k));
	     // rx.setX(xmin+(j*xstep));
	     // rx.setY(ymin+(i*ystep));
	     // rx.setZ(zmin+(k*zstep));
	     // phi=pi;
	     // rx.setRThetaPhi(j*m, theta, phi);
	     // rx.setZ(rx.z()+3.5*m);
	     //	radio->setTxPos(tx);
	     double n=100*pow(10, l);
	     radio->setRxPos(rx);
	     radio->setNPrimaries(n);
	     //this way we have a different run for each event, and so we can make a threshold curve as a function of position and primary energy
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     runManager->BeamOn(1);
	     //	theta+=2*pi/num;
	   }
	 }
       }
     }
   }
    // else{
    
    // }
  
    // else{
    
    // }
  

  
     if(macro.contains("collisionmacro.mac")){
    auto gpsDat=G4GeneralParticleSourceData::Instance();
    auto gps = gpsDat->GetCurrentSource();
    radio->setPrimaryEnergy(gps->GetParticleEnergy());
 
    //    for(int i=15;i>1;i--){
      for(int i=2;i<16;i++){
      radio->setNPrimaries(pow(10, i));
      //      radio->setScaleByEnergy(1);
      //	cout<<rx.phi()<<" "<<radio->rx.phi()<<endl;
    	runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    	// runManager->BeamOn(1);
    }
  }
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
	

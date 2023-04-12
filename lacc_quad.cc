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
// 	Author: yskim
//

// G4Runmanager and mandatory classes


#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "ActionInitialization.hh"
#include "LacaDesignFactor.hh"

// Randomize class to set seed number
#include "Randomize.hh"

// UI and visualization classes
#include "G4UImanager.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "G4Threading.hh"
#define CCMODE 2	// 0: scatter only, 1: absorber same as scatter, 2: absorber similar to scatter
#define SCINTON 0	// 0: no light, 1: simulate lights
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
namespace{
void PrintUsage(){
	G4cerr << " Usage: " << G4endl
			<< " ProjectName [-option1 value1] [-option2 value2] ..." << G4endl;
	G4cerr << "\t--- Option lists ---"
			<< "\n\t[-m] <Set macrofile> default: vis.mac, inputtype: String"
			<< "\n\t[-r] <Set seednum> default: 'current time', inputtype: int"
			<< "\n\t[-n] <Set ProgressN> default: 1000, inputtype: int"
			<< "\n\t[-p] <Set physics> default: 'code', inputtype: String"
			<< "\n\t[-t] <Set nThreads> default: 1, inputtype: int"
			<< G4endl;
	G4cerr << "\tNote: -t option is available only for multi-threaded mode. Max: "
			<< G4Threading::G4GetNumberOfCores() << G4endl;
}
}

G4int nThreads;
G4int ProgressN;
G4String Seed;

int main(int argc, char** argv){
	// Evaluate arguments
	if (argc>21){
		PrintUsage();
		return 1;
	}
	G4String macro;
	G4String physName = "code";
	G4int SeedNumber = time(NULL);
	ProgressN = 1000000;
	Seed = "X";

	LacaDesignFactor testdesign;
	G4double m2d=-1;
	G4double mask_t;
    G4double mask_w;
    G4int rank;
    G4int mode=0;


#ifdef G4MULTITHREADED
	nThreads = 1;
#endif
	for(G4int i=1; i<argc; i=i+2){
		if(G4String(argv[i])=="-m") macro = argv[i+1];
		else if(G4String(argv[i])=="-p") physName = argv[i+1];
		else if (G4String(argv[i])=="-r"){
			Seed = argv[i+1];
			SeedNumber = atoi(Seed);
		}
		else if(G4String(argv[i])=="-n") ProgressN = atoi(argv[i+1]);
		else if(G4String(argv[i])=="-m2d") {
					m2d = atof(argv[i+1]);
					testdesign.setM2D_Dist(m2d);
				}
				else if(G4String(argv[i])=="-mt") {
					mask_t = atof(argv[i+1]);
					testdesign.setMask_T(mask_t);
				}
				else if(G4String(argv[i])=="-npix") {
					rank = atoi(argv[i+1]);
					testdesign.setMask_nPix(rank);
				}
				else if(G4String(argv[i])=="-mw") {
					mask_w = atoi(argv[i+1]);
					testdesign.setMask_W(mask_w);
				}
				else if(G4String(argv[i])=="-mode"){
					mode = atoi(argv[i+1]);
					testdesign.setmode(mode);
				}
#ifdef G4MULTITHREADED
		else if(G4String(argv[i])=="-t"){
			nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
		}
#endif
		else{
			PrintUsage();
			return 1;
		}
	}

	// Seed number setting
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(SeedNumber);

	// Construct runmanager
#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	if (nThreads>0) runManager->SetNumberOfThreads(nThreads);
	else runManager->SetNumberOfThreads(1);
#else
	G4RunManager* runManager = new G4RunManager;
#endif

	DetectorConstruction* _DC = new DetectorConstruction(CCMODE,testdesign);
	runManager->SetUserInitialization(_DC);
	G4VModularPhysicsList* phys;
	if(physName=="code"){
		phys = new PhysicsList(SCINTON);
	}else{
		G4PhysListFactory factory;
		phys = factory.GetReferencePhysList(physName);
	}
	runManager->SetUserInitialization(phys);
	runManager->SetUserInitialization(new ActionInitialization(CCMODE, SCINTON,_DC));
	runManager->Initialize();

	// Construct UI and visualization manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4VisManager* visManager = new G4VisExecutive();
	visManager->Initialize();

	if(macro.size()){
		// Batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}else{
		// interactive mode: define UI session
		G4UIExecutive* ui = new G4UIExecutive(argc, argv, "");
		UImanager->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
		delete ui;
	}

	// Free the store
	delete visManager;
	delete runManager;

	return 0;
}

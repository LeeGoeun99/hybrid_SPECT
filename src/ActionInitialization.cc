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

#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction_PG.hh"
#include "PrimaryGeneratorAction_GPS.hh"

#include "CCRunAction.hh"
#include "CCEventAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include <fstream>
#include <ostream>
using namespace std;


/*ActionInitialization::ActionInitialization(int mode, bool scint)
:G4VUserActionInitialization(),
 fmode(mode),
 fscint(scint){

}*/
ActionInitialization::ActionInitialization(int mode, bool scint, DetectorConstruction* _DC)
:G4VUserActionInitialization(),
 fmode(mode),
  fscint(scint),
  fDetectorConstruction(_DC){
	PGmode=_DC->GetPGmode();
//	mode=1;
//	num = 0;
//	ifstream ifp;
//	ifp.open("PhaseSpace.csv");

//	while(!ifp.eof())
//	{
//		int evt;
//		char delim;
//		double _vx, _vy, _vz;
//		double _energy;
//		ifp >> evt >> delim >> _vx >> delim >> _vy >> delim >> _vz >> delim >> _energy;
//		vx.push_back(_vx);
//		vy.push_back(_vy);
//		vz.push_back(_vz);
//		energy.push_back(_energy);
//		num++;
		//		G4cout << energy[i] << ", " << vx[i] << G4endl;
//	}

//	ifp.close();
}

ActionInitialization::~ActionInitialization(){

}

void ActionInitialization::BuildForMaster() const{
	if(fmode)
		SetUserAction(new CCRunAction());
	else
		SetUserAction(new RunAction());
}

void ActionInitialization::Build() const{
/*	if (PGmode==1)	SetUserAction(new PrimaryGeneratorAction_PG(fDetectorConstruction));
	else			SetUserAction(new PrimaryGeneratorAction_GPS());*/
	//SetUserAction(new PrimaryGeneratorAction_GPS());
	SetUserAction(new PrimaryGeneratorAction_PG(fDetectorConstruction));

	if(fmode){
		SetUserAction(new CCRunAction());
		SetUserAction(new CCEventAction(fscint));
	}else{
		SetUserAction(new RunAction());
		SetUserAction(new EventAction(fscint));
	}

//	SetUserAction(new TrackingAction());
//	SetUserAction(new SteppingAction());
}

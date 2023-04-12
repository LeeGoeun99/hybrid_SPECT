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

#include "PrimaryGeneratorAction_PG.hh"
#include "G4Event.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction_PG::PrimaryGeneratorAction_PG(DetectorConstruction* _fDetectorConstruction)
:G4VUserPrimaryGeneratorAction(),
 fDetectorConstruction(_fDetectorConstruction)
 {
 	fPrimary = new G4ParticleGun();
 	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	fPrimary->SetParticleDefinition(particleTable->FindParticle("gamma"));

	//G4double Energy[5] = {59.5,356,662,1170,1330};
 }

PrimaryGeneratorAction_PG::~PrimaryGeneratorAction_PG(){
	delete fPrimary;
}

void PrimaryGeneratorAction_PG::GeneratePrimaries(G4Event* anEvent)
{
	G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
	// point source PG
/*	fPrimary->SetParticleEnergy(662*keV);
	auto direction = G4RandomDirection(cos(8*pi/180));
	fPrimary->SetParticleMomentumDirection(direction); 
	fPrimary->SetParticlePosition(G4ThreeVector(0, 0, 5*m));
	fPrimary->GeneratePrimaryVertex(anEvent);
*/
// bkg PG
/*
	G4int n = fDetectorConstruction->GetPSN();
	G4int randN = (G4int)floor(n*G4UniformRand());
	G4double gammaE = fDetectorConstruction->GetPSE(randN);
	G4double vecX = fDetectorConstruction->GetPSX(randN);
	G4double vecY = fDetectorConstruction->GetPSY(randN);
	G4double vecZ = fDetectorConstruction->GetPSZ(randN);
	G4ThreeVector vec(vecX, vecY, vecZ);
	G4double x, y, z;
	G4double randN2 = G4UniformRand();

	if(randN2 < (1.0/6.0)){
		x = (3*G4UniformRand() - 1.5)*m;
		y = -1.5*m;
		z = (3*G4UniformRand() - 1.5)*m;
	}
	else if(randN2 < (2.0/6.0)){
		x = (3*G4UniformRand() - 1.5)*m;
		y = 1.5*m;
		z = (3*G4UniformRand() - 1.5)*m;
		vec.rotateX(180*deg);
	}
	else if(randN2 < (3.0/6.0)){
		x = -1.5*m;
		y = (3*G4UniformRand() - 1.5)*m;
		z = (3*G4UniformRand() - 1.5)*m;
		vec.rotateZ(-90*deg);
	}
	else if(randN2 < (4.0/6.0)){
		x = 1.5*m;
		y = (3*G4UniformRand() - 1.5)*m;
		z = (3*G4UniformRand() - 1.5)*m;
		vec.rotateZ(90*deg);
	}
	else if(randN2 < (5.0/6.0)){
		x = (3*G4UniformRand() - 1.5)*m;
		y = (3*G4UniformRand() - 1.5)*m;
		z = -1.5*m;
		vec.rotateX(90*deg);
	}
	else{
		x = (3*G4UniformRand() - 1.5)*m;
		y = (3*G4UniformRand() - 1.5)*m;
		z = 1.5*m;
		vec.rotateX(-90*deg);
	}
	fPrimary->SetParticlePosition(G4ThreeVector(x, y, z));
	fPrimary->SetParticleMomentumDirection(vec);
	fPrimary->SetParticleEnergy(gammaE);
	double eventN = anEvent->GetEventID();
	double emitT_interval = 1250* ns; 
	fPrimary->GeneratePrimaryVertex(anEvent);
*/
	G4double innersource_xpos = 0*mm;
	G4double innersource_ypos = 0*mm;
	G4double innersource_zpos = 0*mm;
	G4double x = 406.5;
	//source direction is set to be half sphere
	fPrimary->SetParticleEnergy(662.*keV);
	auto direction = G4RandomDirection(cos(90*pi/180));
	fPrimary->SetParticleMomentumDirection(direction); 
	innersouce_length = sqrt(innersource_xpos*innersource_xpos+innersource_zpos*innersource_zpos);
	ceta = atan2(innersource_xpos,innersource_zpos);
	ceta_degree = ceta * 180 / PI;
	degree_pos = ceta_degree + runID * 1;
	xpos = x*mm + innersouce_length*sin(degree_pos * PI/180);
	zpos = innersouce_length*cos(degree_pos*PI/180);
	fPrimary->SetParticlePosition(G4ThreeVector(xpos, innersource_ypos, zpos));
   	fPrimary->GeneratePrimaryVertex(anEvent);
}


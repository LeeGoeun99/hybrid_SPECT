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

#include "DEPosSD.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#define SIGPROCT_ns 15000
#define ABSOL(x) ((x<0)? -x:x)

DEPosSD::DEPosSD(G4String name):
	G4VSensitiveDetector(name),
	fHCID(-1){
	collectionName.insert("DEPosColl");
}

DEPosSD::~DEPosSD(){
}

void DEPosSD::Initialize(G4HCofThisEvent* hce){
	fHitsCollection = new DEPosHitsCollection(SensitiveDetectorName, collectionName[0]);
	if(fHCID<0){
		fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
	}
	hce->AddHitsCollection(fHCID, fHitsCollection);

//	DEPosHit* hit = new DEPosHit();
//	fHitsCollection->insert(hit);
}

G4bool DEPosSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
	G4String ptcName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
	if(ptcName=="opticalphoton") return true;

	G4double currDE = aStep->GetTotalEnergyDeposit();
	if(currDE==0.) return true;

	G4int DetNum = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0)-400;
	/*G4cout << DetNum <<G4endl;*/
	G4double interT = aStep->GetTrack()->GetGlobalTime();
	G4int n;
	DEPosHit* hit;
	for(n=0; n<fHitsCollection->GetSize(); n++){
		if(DetNum==(*fHitsCollection)[n]->GetDet()
				&& ABSOL(interT - (*fHitsCollection)[n]->GetT())<=SIGPROCT_ns*ns){
			hit = (*fHitsCollection)[n];
			break;
		}
	}
	if(n==fHitsCollection->GetSize())
	{
		hit = new DEPosHit(DetNum);
		hit->SetT(interT);
		fHitsCollection->insert(hit);
	}

	G4ThreeVector postPos = aStep->GetPostStepPoint()->GetPosition();
	hit->UpdateDEPos(postPos, currDE);

	return true;
}

void DEPosSD::EndOfEvent(G4HCofThisEvent*){
}

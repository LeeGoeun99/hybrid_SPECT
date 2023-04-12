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

#include "PhotCntSD.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

PhotCntSD::PhotCntSD(G4String name):
	G4VSensitiveDetector(name),
	fHCID(-1){
	collectionName.insert("PhotCntColl");
}

PhotCntSD::~PhotCntSD(){
}

void PhotCntSD::Initialize(G4HCofThisEvent* hce){
	fHitsCollection = new PhotCntHitsCollection(SensitiveDetectorName, collectionName[0]);
	if(fHCID<0){
		fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
	}
	hce->AddHitsCollection(fHCID, fHitsCollection);

//	PhotCntHit* hit = new PhotCntHit();
//	fHitsCollection->insert(hit);
}

G4bool PhotCntSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){
	G4String ptcName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
	if(ptcName!="opticalphoton") return true;

	G4double currDE = aStep->GetTotalEnergyDeposit();
	if(currDE==0.) return true;

	G4int DetNum = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(6);
	G4int n;
	PhotCntHit* hit;
	for(n=0; n<fHitsCollection->GetSize(); n++){
		if(DetNum==(*fHitsCollection)[n]->GetDet()){
			hit = (*fHitsCollection)[n];
			break;
		}
	}
	if(n==fHitsCollection->GetSize()){
		hit = new PhotCntHit(DetNum);
		fHitsCollection->insert(hit);
	}

	G4int PMTid = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(3)-100;
	hit->addCnt(PMTid);

	return true;
}

void PhotCntSD::EndOfEvent(G4HCofThisEvent*){
}

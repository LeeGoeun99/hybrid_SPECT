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

#include "EventAction.hh"

#include "DEPosHit.hh"
#include "PhotCntHit.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4CsvAnalysisManager.hh"

EventAction::EventAction(bool scint)
:G4UserEventAction(),
 fEPHCID(-1),
 fPCHCID(-1),
 fscint(scint){
}

EventAction::~EventAction(){
}

void EventAction::BeginOfEventAction(const G4Event*){
}

void EventAction::EndOfEventAction(const G4Event* anEvent){
	if(fEPHCID==-1){
		fEPHCID = G4SDManager::GetSDMpointer()->GetCollectionID("DEPos/DEPosColl");
	}
	if(fPCHCID==-1){
		fPCHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PhotCnt/PhotCntColl");
	}

	G4HCofThisEvent* hce = anEvent->GetHCofThisEvent();
	if(!hce){
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found.\n";
		G4Exception("EndOfEventAction()", "Code001", JustWarning, msg);
		return;
	}

	DEPosHitsCollection* epHC
	= static_cast<DEPosHitsCollection*>(hce->GetHC(fEPHCID));
	PhotCntHitsCollection* pcHC
	= static_cast<PhotCntHitsCollection*>(hce->GetHC(fPCHCID));

	if (!epHC||!pcHC){
		G4ExceptionDescription msg;
		msg << "Some of hits collections of this event not found.\n";
		G4Exception("EndOfEventAction()", "Code001", JustWarning, msg);
		return;
	}

	G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();

	G4int ep_n_hit = epHC->entries();
	G4int pc_n_hit = pcHC->entries();
//	if(ep_n_hit!=pc_n_hit) return;

	for(G4int i=0; i<ep_n_hit; i++){
		DEPosHit* ephit = (*epHC)[i];
		G4double GammaE = G4RunManager::GetRunManager()->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
		if(ephit->GetDE()!=GammaE) return;	// for LUT generation (chk for phot evt)
		if(ephit->GetDE()<=0) return;

		analysisManager->FillNtupleIColumn(0, 0, anEvent->GetEventID());
		analysisManager->FillNtupleDColumn(0, 1, ephit->GetPos().x()/mm);
		analysisManager->FillNtupleDColumn(0, 2, ephit->GetPos().y()/mm);
		analysisManager->FillNtupleDColumn(0, 3, ephit->GetPos().z()/mm);
		analysisManager->FillNtupleDColumn(0, 4, ephit->GetDE()/MeV);
		analysisManager->AddNtupleRow(0);
	}

	if(fscint){
		G4int nPMTs[2] = {36, 22};
		for(G4int i=0; i<pc_n_hit; i++){
			PhotCntHit* pchit = (*pcHC)[i];
			analysisManager->FillNtupleIColumn(1, 0, anEvent->GetEventID());
			for(G4int n=0; n<nPMTs[0]; n++){
				analysisManager->FillNtupleIColumn(1, n+1, pchit->GetCnt(n+1));
			}
			analysisManager->AddNtupleRow(1);
		}
	}
}


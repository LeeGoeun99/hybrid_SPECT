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

#include "CCEventAction.hh"
#include "DEPosHit.hh"
#include "PhotCntHit.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4CsvAnalysisManager.hh"
using std::vector;
using namespace std;
#define COINWINDOW_ns 400
#define ABSOL(x) ((x<0)? -x:x)

CCEventAction::CCEventAction(bool scint)
:G4UserEventAction(),
 fEPHCID(-1),
 fPCHCID(-1),
 fscint(scint){
}
CCEventAction::~CCEventAction(){
}
vector<G4int>duplication(vector<G4int>& a, vector <G4int>& b)
	{   vector<G4int>::iterator iter;
	    vector<G4int>::iterator iter_b;
	    vector<G4int> c = a; //a의 값 복사
	    for (iter_b = b.begin(); iter_b != b.end(); iter_b++)
	    {
	        for (iter = c.begin(); iter != c.end();)
	        {
	            if (*iter == *iter_b)
	                iter = c.erase(iter); //중복 제거
	            else
	                iter++;
	        }
	    }
	    return c; //결과 반환
	}

void CCEventAction::BeginOfEventAction(const G4Event*){
}

void CCEventAction::EndOfEventAction(const G4Event* anEvent)
{
	if(fEPHCID==-1)
	{
		fEPHCID = G4SDManager::GetSDMpointer()->GetCollectionID("DEPos/DEPosColl");
	}
	if(fPCHCID==-1)
	{
		fPCHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PhotCnt/PhotCntColl");
	}

	G4HCofThisEvent* hce = anEvent->GetHCofThisEvent();
	if(!hce)
	{
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found.\n";
		G4Exception("EndOfEventAction()", "Code001", JustWarning, msg);
		return;
	}

	DEPosHitsCollection* epHC
	= static_cast<DEPosHitsCollection*>(hce->GetHC(fEPHCID));
	PhotCntHitsCollection* pcHC
	= static_cast<PhotCntHitsCollection*>(hce->GetHC(fPCHCID));

	if (!epHC||!pcHC)
	{
		G4ExceptionDescription msg;
		msg << "Some of hits collections of this event not found.\n";
		G4Exception("EndOfEventAction()", "Code001", JustWarning, msg);
		return;
	}

	G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
	G4int ep_n_hit = epHC->entries();
	G4int pc_n_hit = pcHC->entries();
	vector<DEPosHit*> epHCtmp;
	vector<G4int>epDetIDM;
	vector<G4int>AdetIDM={1,2, 3, 4, 5, 6, 7, 8};
	vector<G4int>ScatIDM;
	vector<G4int>AbsoIDM;

	for(G4int i=0; i<ep_n_hit; i++)	{
		DEPosHit* ephit = (*epHC)[i];
		if(ephit->GetDE()<=0) continue;
		G4int epDetID = ephit->GetDet();
		epHCtmp.push_back(ephit);
		epDetIDM.push_back(epDetID);

		if(epDetID<=4){
		ScatIDM.push_back(epDetID);
		analysisManager->FillNtupleIColumn(2, 0, anEvent->GetEventID());
		analysisManager->FillNtupleDColumn(2, 1, ephit->GetPos().x()/mm);
		analysisManager->FillNtupleDColumn(2, 2, ephit->GetPos().y()/mm);
		analysisManager->FillNtupleDColumn(2, 3, ephit->GetPos().z()/mm);
		analysisManager->FillNtupleDColumn(2, 4, ephit->GetDE()/MeV);
		analysisManager->FillNtupleDColumn(2, 5, ephit->GetT()/ns);
		analysisManager->AddNtupleRow(2);
		}
		else{
		AbsoIDM.push_back(epDetID);
		analysisManager->FillNtupleIColumn(3, 0, anEvent->GetEventID());
		analysisManager->FillNtupleDColumn(3, 1, ephit->GetPos().x()/mm);
		analysisManager->FillNtupleDColumn(3, 2, ephit->GetPos().y()/mm);
		analysisManager->FillNtupleDColumn(3, 3, ephit->GetPos().z()/mm);
		analysisManager->FillNtupleDColumn(3, 4, ephit->GetDE()/MeV);
		analysisManager->FillNtupleDColumn(3, 5, ephit->GetT()/ns);
		analysisManager->AddNtupleRow(3);}

		G4int Scattersize=ScatIDM.size();
		G4int Absorsize=AbsoIDM.size();

		if(i==ep_n_hit-1)
		{
			if(Scattersize>=1 && Absorsize>=1)
			{analysisManager->FillNtupleIColumn(0, 0, anEvent->GetEventID());
			for(G4int x=0;x<8;x++){
			 		 analysisManager->FillNtupleDColumn(0,1+(AdetIDM[x]-1)*5, 0/mm);
			 		 analysisManager->FillNtupleDColumn(0,2+(AdetIDM[x]-1)*5, 0/mm);
			 		 analysisManager->FillNtupleDColumn(0,3+(AdetIDM[x]-1)*5, 0/mm);
			 		 analysisManager->FillNtupleDColumn(0,4+(AdetIDM[x]-1)*5, 0/MeV);
			 		 analysisManager->FillNtupleDColumn(0,5+(AdetIDM[x]-1)*5, 0/ns);
			 		}
			 for(G4int l=0; l<epHCtmp.size(); l++)	
			 	{
				  	if(epDetID!=(epHCtmp)[l]->GetDet()	&& ABSOL(ephit->GetT() - (epHCtmp)[l]->GetT())<=COINWINDOW_ns*ns)
			 	 	{	analysisManager->FillNtupleDColumn(0, 1+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().x()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 2+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().y()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 3+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().z()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 4+(epDetIDM[l]-1)*5, epHCtmp[l]->GetDE()/MeV);
			 	 		analysisManager->FillNtupleDColumn(0, 5+(epDetIDM[l]-1)*5, epHCtmp[l]->GetT()/ns);  
					}
			 	 	if(epDetID==(epHCtmp)[l]->GetDet())
			 	 	{
			 	 		analysisManager->FillNtupleDColumn(0, 1+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().x()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 2+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().y()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 3+(epDetIDM[l]-1)*5, epHCtmp[l]->GetPos().z()/mm);
			 	 		analysisManager->FillNtupleDColumn(0, 4+(epDetIDM[l]-1)*5, epHCtmp[l]->GetDE()/MeV);
			 	 		analysisManager->FillNtupleDColumn(0, 5+(epDetIDM[l]-1)*5, epHCtmp[l]->GetT()/ns);
			 	 		analysisManager->AddNtupleRow(0); 
					}
			 	}
			}
	if(fscint)	{
		G4int nPMTs[2] = {36, 36};
		G4double sumCnt = 0;
		vector<G4int>PMTM;

			for(G4int i=0; i<pc_n_hit; i++)
			{
				PhotCntHit* pchit = (*pcHC)[i];
				G4int pcDetID = pchit->GetDet();
				PMTM.push_back(pcDetID);
				if(i!=pc_n_hit-1)	{
					for(G4int r=0; r<PMTM.size(); r++)	{
						for(G4int n=0; n<nPMTs[PMTM[r]-1]; n++){
							sumCnt += pchit->GetCnt(n+1);}
						if(sumCnt>10){
						analysisManager->FillNtupleIColumn(1, 0, anEvent->GetEventID());
						for(G4int n=0; n<nPMTs[PMTM[r]-1]; n++){
						analysisManager->FillNtupleIColumn(1, n+1+(PMTM[r]-1)*nPMTs[0], pchit->GetCnt(n+1));}
								}
					}
				}
				if(i==pc_n_hit-1)	{
					for(G4int n=0; n<nPMTs[pcDetID-1]; n++)	{
						sumCnt += pchit->GetCnt(n+1);}
					if(sumCnt>10){
						analysisManager->FillNtupleIColumn(1, 0, anEvent->GetEventID());
						for(G4int n=0; n<nPMTs[pcDetID-1]; n++)	{
						analysisManager->FillNtupleIColumn(1, n+1+(pcDetID-1)*nPMTs[0], pchit->GetCnt(n+1));}
						analysisManager->AddNtupleRow(1);
								}
				}
		}
		}
	}}}

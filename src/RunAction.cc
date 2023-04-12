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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4Threading.hh"

extern G4int nThreads;
extern G4int ProgressN;
extern G4String Seed;

RunAction::RunAction()
:G4UserRunAction(){
	G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
	analysisManager->SetVerboseLevel(0);
	G4cout << "Using " << analysisManager->GetType() << G4endl;

	analysisManager->CreateNtuple("DEPos", "DEPos");
	analysisManager->CreateNtupleIColumn(0, "EvtN");
	analysisManager->CreateNtupleDColumn(0, "PosX");
	analysisManager->CreateNtupleDColumn(0, "PosY");
	analysisManager->CreateNtupleDColumn(0, "PosZ");
	analysisManager->CreateNtupleDColumn(0, "Edep");
	analysisManager->FinishNtuple(0);

	analysisManager->CreateNtuple("PhotCnt", "PhotCnt");
	analysisManager->CreateNtupleIColumn(1, "EvtN");
	for(G4int i=0; i<36; i++){
		std::stringstream colName;
		colName << "PMT_" << i+1;
		analysisManager->CreateNtupleIColumn(1, colName.str());
	}
	analysisManager->FinishNtuple(1);
}

RunAction::~RunAction(){
	delete G4CsvAnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*){
	G4RunManager::GetRunManager()->SetPrintProgress(ProgressN);

	G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
	std::stringstream fNameInit;
//	fNameInit << "output/seed" << Seed;
	fNameInit << "output/seed" << Seed << "_pos" << std::setfill('0') << std::setw(3) << runID;

	std::vector < G4String > NtupleName;
	NtupleName.push_back("DEPos");
	NtupleName.push_back("PhotCnt");

	fFileName.clear();
	for(std::vector < G4String >::iterator iter=NtupleName.begin(); iter!=NtupleName.end(); iter++){
		std::stringstream ss;
		ss << fNameInit.str() << "_nt_" << *iter;
		fFileName.push_back(ss.str());
	}

	G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
	analysisManager->OpenFile(fNameInit.str());
}

void RunAction::EndOfRunAction(const G4Run* aRun){
	G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
	analysisManager->CloseFile();

	G4int nEvents = aRun->GetNumberOfEvent();
	if(IsMaster()){
		G4cout << "\n----- End of global run -----"
				<< "\n The run was " << nEvents << " events." << G4endl;

		for(std::vector < G4String >::iterator iter=fFileName.begin(); iter!=fFileName.end(); iter++){
			G4cout << ">> File merging... (" << *iter << ")" << G4endl;
			std::ofstream ofsMerge(*iter + ".csv");
			for(G4int i=0; i<nThreads; i++){
				std::stringstream nameEachFile;
				nameEachFile << *iter << "_t" << i << ".csv";
				std::ifstream ifsEachFile(nameEachFile.str());
				if(!ifsEachFile.is_open()){
					G4cout << "!!!!!!!!! noFileError: " << nameEachFile.str() << G4endl;
					continue;
				}
				G4cout << "> File: " << nameEachFile.str() << G4endl;
				while(!ifsEachFile.eof()){
					std::string str;
					std::getline(ifsEachFile, str);
					str.erase(std::remove(str.begin(), str.end(), '\x00'), str.end());
					if((str.find(",")!=std::string::npos))
						ofsMerge << str << G4endl;
				}
				ifsEachFile.close();
			}
			ofsMerge.close();
		}
		G4cout << ">> File merging process completed." << G4endl;
	}
}


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
//    Author: yskim
//

#include "CCRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4Threading.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

extern G4int nThreads;
extern G4int ProgressN;
extern G4String Seed;

CCRunAction::CCRunAction()
:G4UserRunAction(){
   G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
   analysisManager->SetVerboseLevel(0);
   G4cout << "Using " << analysisManager->GetType() << G4endl;
// modify need
   analysisManager->CreateNtuple("DEPos", "DEPos");
   analysisManager->CreateNtupleIColumn(0, "EvtN");
   analysisManager->CreateNtupleDColumn(0, "PosX1");
   analysisManager->CreateNtupleDColumn(0, "PosY1");
   analysisManager->CreateNtupleDColumn(0, "PosZ1");
   analysisManager->CreateNtupleDColumn(0, "Edep1");
   analysisManager->CreateNtupleDColumn(0, "Time1");
   analysisManager->CreateNtupleDColumn(0, "PosX2");
   analysisManager->CreateNtupleDColumn(0, "PosY2");
   analysisManager->CreateNtupleDColumn(0, "PosZ2");
   analysisManager->CreateNtupleDColumn(0, "Edep2");
   analysisManager->CreateNtupleDColumn(0, "Time2");

   analysisManager->CreateNtupleDColumn(0, "PosX3");
   analysisManager->CreateNtupleDColumn(0, "PosY3");
   analysisManager->CreateNtupleDColumn(0, "PosZ3");
   analysisManager->CreateNtupleDColumn(0, "Edep3");
   analysisManager->CreateNtupleDColumn(0, "Time3");
   analysisManager->CreateNtupleDColumn(0, "PosX4");
   analysisManager->CreateNtupleDColumn(0, "PosY4");
   analysisManager->CreateNtupleDColumn(0, "PosZ4");
   analysisManager->CreateNtupleDColumn(0, "Edep4");
   analysisManager->CreateNtupleDColumn(0, "Time4");

   analysisManager->CreateNtupleDColumn(0, "PosX5");
   analysisManager->CreateNtupleDColumn(0, "PosY5");
   analysisManager->CreateNtupleDColumn(0, "PosZ5");
   analysisManager->CreateNtupleDColumn(0, "Edep5");
   analysisManager->CreateNtupleDColumn(0, "Time5");
   analysisManager->CreateNtupleDColumn(0, "PosX6");
   analysisManager->CreateNtupleDColumn(0, "PosY6");
   analysisManager->CreateNtupleDColumn(0, "PosZ6");
   analysisManager->CreateNtupleDColumn(0, "Edep6");
   analysisManager->CreateNtupleDColumn(0, "Time6");

   analysisManager->CreateNtupleDColumn(0, "PosX7");
   analysisManager->CreateNtupleDColumn(0, "PosY7");
   analysisManager->CreateNtupleDColumn(0, "PosZ7");
   analysisManager->CreateNtupleDColumn(0, "Edep7");
   analysisManager->CreateNtupleDColumn(0, "Time7");
   analysisManager->CreateNtupleDColumn(0, "PosX8");
   analysisManager->CreateNtupleDColumn(0, "PosY8");
   analysisManager->CreateNtupleDColumn(0, "PosZ8");
   analysisManager->CreateNtupleDColumn(0, "Edep8");
   analysisManager->CreateNtupleDColumn(0, "Time8");

   analysisManager->FinishNtuple(0);

   analysisManager->CreateNtuple("PhotCnt", "PhotCnt");
   analysisManager->CreateNtupleIColumn(1, "EvtN");
   for(G4int i=0; i<72; i++){
      std::stringstream colName;
      colName << "D1PMT_" << i+1;
      analysisManager->CreateNtupleIColumn(1, colName.str());
   }
   for(G4int i=0; i<72; i++){
      std::stringstream colName;
      colName << "D2PMT_" << i+1;
      analysisManager->CreateNtupleIColumn(1, colName.str());
   }
   analysisManager->FinishNtuple(1);

   analysisManager->CreateNtuple("D1Sing", "D1Sing");
   analysisManager->CreateNtupleIColumn(2, "EvtN");
   analysisManager->CreateNtupleDColumn(2, "PosX1");
   analysisManager->CreateNtupleDColumn(2, "PosY1");
   analysisManager->CreateNtupleDColumn(2, "PosZ1");
   analysisManager->CreateNtupleDColumn(2, "Edep1");
   analysisManager->CreateNtupleDColumn(2, "Time1");
   analysisManager->FinishNtuple(2);

   analysisManager->CreateNtuple("D2Sing", "D2Sing"); 
   analysisManager->CreateNtupleIColumn(3, "EvtN");
   analysisManager->CreateNtupleDColumn(3, "PosX2");
   analysisManager->CreateNtupleDColumn(3, "PosY2");
   analysisManager->CreateNtupleDColumn(3, "PosZ2");
   analysisManager->CreateNtupleDColumn(3, "Edep2");
   analysisManager->CreateNtupleDColumn(3, "Time2");
   analysisManager->FinishNtuple(3);

}

CCRunAction::~CCRunAction(){
   delete G4CsvAnalysisManager::Instance();
}

void CCRunAction::BeginOfRunAction(const G4Run*){
   G4RunManager::GetRunManager()->SetPrintProgress(ProgressN);

   G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
   std::stringstream fNameInit;

   fNameInit << "output/seed" << Seed << "_" << runID;

   std::vector < G4String > NtupleName;
   NtupleName.push_back("DEPos");
   NtupleName.push_back("PhotCnt");
   NtupleName.push_back("D1Sing");
   NtupleName.push_back("D2Sing");


   fFileName.clear();
   for(std::vector < G4String >::iterator iter=NtupleName.begin(); iter!=NtupleName.end(); iter++){
      std::stringstream ss;
      ss << fNameInit.str() << "_nt_" << *iter;
      fFileName.push_back(ss.str());
   }

   G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();
   analysisManager->OpenFile(fNameInit.str());
}

void CCRunAction::EndOfRunAction(const G4Run* aRun)
{
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

         for(G4int j=0;j<nThreads;j++){
            G4cout << "Delete thread files..." << *iter << "_t" << j << G4endl;
            std::stringstream nameEachFile2;
            nameEachFile2 << *iter << "_t" << j << ".csv";
            G4String fileName = nameEachFile2.str();
            std::remove(fileName);
         }
      }
      G4cout << ">> File merging process completed." << G4endl;

      auto pv_Cyl = G4PhysicalVolumeStore::GetInstance()->GetVolume("rot_World");
      auto rotmat = pv_Cyl->GetRotation();
      G4int rot_degree = -1;

      rotmat->rotateZ(rot_degree* deg); //rot_degree=5.*deg;
      pv_Cyl->SetRotation(rotmat);

   }
}
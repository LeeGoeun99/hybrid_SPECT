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

#include "MaskParam.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

MaskParam::MaskParam(G4double _sizeX, G4double _sizeY, G4int _nX, G4int _nY, G4String fileName)
:sizeX(_sizeX),
 sizeY(_sizeY),
 nX(_nX),
 nY(_nY){
    mat_Air = G4Material::GetMaterial("G4_AIR");
    mat_W = G4Material::GetMaterial("G4_W");

    graySol = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	graySol->SetVisibility(true);
	graySol->SetForceSolid(true);
    grayWired = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	grayWired->SetVisibility(true);
	grayWired->SetForceWireframe(true);

    std::ifstream ifp;
    ifp.open(fileName.c_str());
    if(!ifp.is_open()){
        G4cout << "-- MaskParam Error:" << G4endl
            << "---- Unable to open the file." << G4endl;
    }
    for(G4int j=0; j<nY; j++){
        std::vector < G4bool > tmpVec;
        for(G4int i=0; i<nX; i++){
            G4bool tmp;
            ifp >> tmp;
            tmpVec.push_back(tmp);
        }
        mask_Pattern.push_back(tmpVec);
    }
    ifp.close();

    G4cout << "-- Mask pattern:" << G4endl;
    for(G4int j=0; j<nY; j++){
        for(G4int i=0; i<nX; i++){
            G4cout << mask_Pattern[j][i] << " ";
        }
        G4cout << G4endl;
    }
}

MaskParam::~MaskParam(){
}

void MaskParam::ComputeTransformation(const G4int CpNo, G4VPhysicalVolume* physVol) const{
    G4int iX = CpNo%nX;
    G4int iY = (G4int)(CpNo/nX);
    G4double x = -.5*sizeX + .5*sizeX/(G4double)nX + iX*sizeX/(G4double)nX;
    G4double y = -.5*sizeY + .5*sizeY/(G4double)nY + iY*sizeY/(G4double)nY;
    physVol->SetTranslation(G4ThreeVector(x, y, 0.));
}

G4Material* MaskParam::ComputeMaterial(G4VPhysicalVolume* physVol, const G4int CpNo, const G4VTouchable* parentTouch){
    G4int iX = CpNo%nX;
    G4int iY = (G4int)(CpNo/nX);

    if(mask_Pattern[iY][iX]){
        physVol->GetLogicalVolume()->SetVisAttributes(grayWired);
        return mat_Air;
    }else{
        physVol->GetLogicalVolume()->SetVisAttributes(graySol);
        return mat_W;
    }
}

G4int MaskParam::GetNumberOfMaterials() const{
    return 2;
}

G4Material* MaskParam::GetMaterial(G4int CpNo) const{
    G4int iX = CpNo%nX;
    G4int iY = (G4int)(CpNo/nX);
    return mask_Pattern[iY][iX] ? mat_W : mat_Air;
}

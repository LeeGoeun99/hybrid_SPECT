//
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
// Author: yskim
//

#include "DEPosHit.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// -- Mandatory -- //
G4ThreadLocal G4Allocator<DEPosHit>* DEPosHitAllocator;

// -- User defined & overriding -- //
DEPosHit::DEPosHit(G4int iDet)
:G4VHit(),
 fDet(iDet),
 fDE(0.),
 fPos(0),
 fT(0){
}

DEPosHit::DEPosHit(G4int iDet, G4double iDE, G4ThreeVector iPos, G4double iT)
:G4VHit(),
 fDet(iDet),
 fDE(iDE),
 fPos(iPos),
 fT(iT){
}

DEPosHit::DEPosHit(const DEPosHit &right)
:G4VHit(){
	*this = right;
}

DEPosHit::~DEPosHit(){
}

const DEPosHit& DEPosHit::operator=(const DEPosHit &right){
	fDet = right.fDet;
	fDE = right.fDE;
	fPos = right.fPos;
	fT = right.fT;
	return *this;
}

int DEPosHit::operator==(const DEPosHit &right) const{
    return (fDet==right.fDet && fDE==right.fDE && fPos==right.fPos && fT==right.fT);
}


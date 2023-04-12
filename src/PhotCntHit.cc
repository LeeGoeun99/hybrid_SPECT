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

#include "PhotCntHit.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// -- Mandatory -- //
G4ThreadLocal G4Allocator<PhotCntHit>* PhotCntHitAllocator;

// -- User defined & overriding -- //
PhotCntHit::PhotCntHit(G4int iDet)
:G4VHit(),
 fDet(iDet),
 fCnts(){
}

PhotCntHit::PhotCntHit(G4int iDet, std::map < G4int,G4int > iCnts)
:G4VHit(),
 fDet(iDet){
	fCnts = iCnts;
}

PhotCntHit::PhotCntHit(const PhotCntHit &right)
:G4VHit(){
	*this = right;
}

PhotCntHit::~PhotCntHit(){
}

const PhotCntHit& PhotCntHit::operator=(const PhotCntHit &right){
	fDet = right.fDet;
	fCnts = right.fCnts;
	return *this;
}

int PhotCntHit::operator==(const PhotCntHit &right) const{
	return (fDet==right.fDet && fCnts==right.fCnts);
}

void PhotCntHit::addCnt(G4int iPMTid, G4int nPhots){
	std::pair < std::map < G4int,G4int >::iterator, G4bool > chk = fCnts.insert(std::pair < G4int,G4int > (iPMTid, 1));
	if(!chk.second){
		(fCnts.find(iPMTid)->second)+=nPhots;
	}
}

void PhotCntHit::addCnts(std::map < G4int,G4int > iCnts){
	for(std::map < G4int,G4int >::iterator iter=iCnts.begin(); iter!=iCnts.end(); iter++){
		addCnt(iter->first, iter->second);
	}
}

G4int PhotCntHit::GetCnt(G4int iPMTid) const{
	std::map < G4int,G4int >::const_iterator iter = fCnts.find(iPMTid);
	if(iter==fCnts.end())
		return 0;
	else
		return iter->second;
}

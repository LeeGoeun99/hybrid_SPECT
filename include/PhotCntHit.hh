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

#ifndef PhotCntHit_hh_
#define PhotCntHit_hh_

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "G4ThreeVector.hh"

class PhotCntHit: public G4VHit
{
public:
	// -- Mandatory -- //
	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	// -- User defined & overriding -- //
	PhotCntHit(G4int iDet);
	PhotCntHit(G4int iDet, std::map < G4int,G4int > iCnts);
	PhotCntHit(const PhotCntHit &right);

	virtual ~PhotCntHit();

	const PhotCntHit& operator=(const PhotCntHit &right);
	int operator==(const PhotCntHit &right) const;

	G4int GetDet(){return fDet;}

	void addCnt(G4int iPMTid, G4int nPhots=1);
	void addCnts(std::map < G4int,G4int > iCnts);
	G4int GetCnt(G4int iPMTid) const;
	std::map < G4int,G4int > GetCnts() const{return fCnts;}
private:
	G4int fDet;
	std::map < G4int,G4int > fCnts;
};

// -- Mandatory -- //
typedef G4THitsCollection<PhotCntHit> PhotCntHitsCollection;
extern G4ThreadLocal G4Allocator<PhotCntHit>* PhotCntHitAllocator;
inline void* PhotCntHit::operator new(size_t){
	if (!PhotCntHitAllocator)
		PhotCntHitAllocator = new G4Allocator<PhotCntHit>;
	return (void*)PhotCntHitAllocator->MallocSingle();
}
inline void PhotCntHit::operator delete(void* aHit){
	PhotCntHitAllocator->FreeSingle((PhotCntHit*) aHit);
}

#endif

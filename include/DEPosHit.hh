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

#ifndef DEPosHit_hh_
#define DEPosHit_hh_

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "G4ThreeVector.hh"
////
#include "math.h"

class DEPosHit: public G4VHit
{
public:
	// -- Mandatory -- //
	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	// -- User defined & overriding -- //
	DEPosHit(G4int iDet);
	DEPosHit(G4int iDet, G4double iDE, G4ThreeVector iPos, G4double iT);
	DEPosHit(const DEPosHit &right);

	virtual ~DEPosHit();

	const DEPosHit& operator=(const DEPosHit &right);
	int operator==(const DEPosHit &right) const;

	G4int GetDet() const{return fDet;}

	void SetDE(G4double iDE){fDE = iDE;}
	G4double GetDE() const{return fDE;}

	void SetPos(G4ThreeVector iPos){fPos = iPos;}
	G4ThreeVector GetPos() const{return fPos;}

	void UpdateDEPos(G4ThreeVector iPos, G4double iDE){
		fPos = (fPos*(fDE) + iPos*iDE)/(fDE + iDE);
		fDE += iDE;
	}

	void SetT(G4double iT){fT = iT;}
	G4double GetT() const{return fT;}

	////////////////
	void SetDelT(){
		G4double idT;
		idT = -log(rand_N) * 1/(emiss_Y * Actvt);
		fdT = idT;
	}

	G4double GetDelT()
	const{	return fdT;}
	//Act need to be modified by BEAMON
private:
	G4int fDet;
	G4double fDE;
	G4ThreeVector fPos;
	G4double fT;

	G4double fdT;
	G4double rand_N = 0.5;
	const G4double emiss_Y = 0.5;
	G4double Actvt; //varies by time



};

// -- Mandatory -- //
typedef G4THitsCollection<DEPosHit> DEPosHitsCollection;
extern G4ThreadLocal G4Allocator<DEPosHit>* DEPosHitAllocator;
inline void* DEPosHit::operator new(size_t){
	if (!DEPosHitAllocator)
		DEPosHitAllocator = new G4Allocator<DEPosHit>;
	return (void*)DEPosHitAllocator->MallocSingle();
}
inline void DEPosHit::operator delete(void* aHit){
	DEPosHitAllocator->FreeSingle((DEPosHit*) aHit);
}

#endif

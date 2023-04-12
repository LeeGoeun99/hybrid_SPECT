/*
 * LacaDesignFactor.hh
 *
 *  Created on: Nov 11, 2018
 *      Author: user
 */

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"


#ifndef INCLUDE_LACADESIGNFACTOR_HH_
#define INCLUDE_LACADESIGNFACTOR_HH_

class LacaDesignFactor{
public:
	LacaDesignFactor();
	virtual ~LacaDesignFactor();
	G4double getM2D_Dist(){return M2D_Dist;};
	G4double getMask_W(){return Mask_W;};
	G4double getMask_T(){return Mask_T;};
	G4int getMask_nPix(){return Mask_nPix;};
	G4int getmode(){return mode;};
	void setM2D_Dist(G4double var){ M2D_Dist=var;};
	void setMask_W(G4double var){ Mask_W=var;};
	void setMask_T(G4double var){ Mask_T=var;};
	void setMask_nPix(G4int var){ Mask_nPix=var;};
	void setmode(G4int var){ mode=var;};

private:

	G4double M2D_Dist = 55;
	G4double Mask_W = 370;
	G4double Mask_T = 6;
	G4int Mask_nPix = 37;
	G4int mode=0;

};
#endif /* INCLUDE_LACADESIGNFACTOR_HH_ */

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

#ifndef DetectorConstruction_hh_
#define DetectorConstruction_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "LacaDesignFactor.hh"

class G4VisAttributes;

class DetectorConstruction: public G4VUserDetectorConstruction{
public:
	DetectorConstruction(G4int mode);
	DetectorConstruction(G4int mode, LacaDesignFactor DesignFactor);

	virtual ~DetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

	void DefineDimensions();
	void DefineMaterials();

	void VisColor();
	G4double GetPSE(G4int n){return energy[n];}
	G4double GetPSX(G4int n){return vx[n];}
	G4double GetPSY(G4int n){return vy[n];}
	G4double GetPSZ(G4int n){return vz[n];}
	G4int GetPSN(){return num;}
	G4int GetPGmode(){return PGmode;}
private:
	std::vector < G4double > energy;
		std::vector < G4double > vx;
		std::vector < G4double > vy;
		std::vector < G4double > vz;
		G4int num;
	G4int fmode;
	G4VPhysicalVolume* World_Phys;
	G4VPhysicalVolume* pv_rotWorld;
	G4VPhysicalVolume* phantom_pv;
	
	G4LogicalVolume *lv_Drum, *lv_innerPhantom, *lv_innerPhantom_air;

	G4double World_Size;
	// For source
	G4double SrcCase_D, SrcCase_T;
	G4double Src2DetDist;
	// For Collimator
	G4double Collimator_W;
	G4double CollimatorFace_T, CollimatorBox_T;
	G4double CollimatorHole_D;
	// For CC
	G4double Sc2Ab_Dist;
	// For detector type 1 (36 ch.)
	G4double PMT1_Gap, PMTbox;
	G4double CrossBar_W, CrossBar_T, CrossBar_H;
	G4double WholeDet1_W, WholeDet1_D, WholeDet1_H;
	G4double Housing1_W, Housing_T, Housing1_V;
	G4double Crystal_W, Crystal_Wf, Crystal_T;
	G4double OptGlue_T;
	G4double OptWin_W, OptWin1_T;
	G4double Paint_T;
	G4double OptGrease_T;
	G4double RearCase_W, RearCase1_D, RearCase1_H, RearCase_T;
	G4double WholePMT_W, WholePMT_T;
	G4double PMTFrontGlass_T;
	G4double PMTPhotoCathode_W, PMTPhotoCathode_T;
	G4double PMTFrontBase_W, PMTFrontBase_T;
	G4double PMTRearBase_D, PMTRearBase_H;
	G4double PMTTail_D, PMTTail_H,PMTTail2_H;
	G4double PMTPCB_D, PMTPCB_T;
	G4double Pillar1_D, Pillar1_H;
	G4double Holder1_W, Holder1_H, Holder1_T;
	G4double rear_Al_T, rear_Al_W,rear_Al_H;
	// For detector type 1' (36 ch.)
	G4double WholeDet11_W, WholeDet11_D, WholeDet11_H;
	G4double Housing11_W, Housing11_T, Housing11_V;
	G4double Crystal11_W, Crystal11_Wf, Crystal11_T;
	G4double OptGlue11_T;
	G4double OptWin11_W, OptWin11_T;
	G4double Paint11_T;
	G4double OptGrease11_T;
	G4double RearCase11_W, RearCase11_H, RearCase11_D, RearCase11_T;
	G4double WholePMT11_W, WholePMT11_H;
	G4double PMTFrontGlass11_T;
	G4double PMTPhotoCathode11_W, PMTPhotoCathode11_T;
	G4double PMTFrontBase11_W, PMTFrontBase11_H, PMTbox11;
	G4double PMTRearBase11_D, PMTRearBase11_H;
	G4double PMTTail11_D, PMTTail11_H, PMTTail12_H ;
	G4double PMTPCB11_D, PMTPCB11_T;
	G4double Pillar11_D, Pillar11_H;
	G4double Holder11_W, Holder11_H, Holder11_T;
	// For detector type 2 (22 ch.)
	G4double WholeDet2_W, WholeDet2_H;
	G4double Housing2_W, Housing2_H;
	G4double Crystal2_W, Crystal2_Wf, Crystal2_T;
	G4double OptGlue2_T;
	G4double OptWin2_W, OptWin2_Wf, OptWin2_T;
	G4double Paint2_T;
	G4double OptGrease2_T;
	G4double WholePMT2_D, WholePMT2_H;
	G4double PMTFrontGlass2_T, PMTAlBase2_T;
	G4double PMTPhotoCathode2_D, PMTPhotoCathode2_T;
	G4double PMTFrontBase2_D, PMTFrontBase2_H;
	G4double PMTRearBase2_D, PMTRearBase2_H;

	G4VisAttributes* blue;
	G4VisAttributes* gray;
	G4VisAttributes* white;
	G4VisAttributes* red;
	G4VisAttributes* yellow;
	G4VisAttributes* green;
	G4VisAttributes* darkGreen;
	G4VisAttributes* darkOrange3;
	G4VisAttributes* skyBlue;

	// For coded mask
	G4double Mask_W, Mask_T;
	G4int Mask_nPix;
	G4double M2D_Dist;
	G4String MaskPatternFile;
	G4int PGmode;

	//Phantom
	G4double Phantom_outerD1;
	G4double Phantom_outerD2;
	G4double Phantom_H;
	G4double phantomholeD;
	G4float phantom_S;
	G4float phantom_H;

};

#endif

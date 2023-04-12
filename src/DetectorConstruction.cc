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
#include "DetectorConstruction.hh"
#include "DEPosSD.hh"
#include "PhotCntSD.hh"
#include "MaskParam.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
#include <iostream>
#include "G4UImanager.hh"
#include <fstream>
#include <ostream>

using namespace std;


DetectorConstruction::DetectorConstruction(G4int mode, LacaDesignFactor DesignFactor)
:G4VUserDetectorConstruction(),
 World_Phys(0), fmode(mode){
    num = 0;
	ifstream ifp;
	// For coded mask
		Mask_W = DesignFactor.getMask_W()*mm;
		Mask_T = DesignFactor.getMask_T()*mm;
		Mask_nPix = DesignFactor.getMask_nPix();
		M2D_Dist = DesignFactor.getM2D_Dist()*mm;
		PGmode = DesignFactor.getmode();
		std::stringstream filename;
		filename<<"MURA/maskpattern_MURA"<<Mask_nPix<<".txt";
		MaskPatternFile= filename.str();
		std::cout<<MaskPatternFile<<std::endl;
	/*
	//Back ground
	ifp.open("PhaseSpace_1.csv");
    while(!ifp.eof())
    {
            int evt;
            char delim;
            double _vx, _vy, _vz;
            double _x, _y, _z;
            double _energy;
            ifp >> evt >> delim >> _vx >> delim >> _vy >> delim >> _vz >> delim >> _energy;

            vx.push_back(_vx);
            vy.push_back(_vy);
            vz.push_back(_vz);
            energy.push_back(_energy);
            num++;
     }
    ifp.close(); */
}

DetectorConstruction::~DetectorConstruction(){
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
	// -- Initialization -- //
	// Dimensions
	DefineDimensions();

	// Materials
	DefineMaterials();
	G4Material* mat_Air = G4Material::GetMaterial("G4_AIR");
	G4Material* mat_PE = G4Material::GetMaterial("G4_POLYETHYLENE"); // Source case
	G4Material* mat_Al = G4Material::GetMaterial("G4_Al");
	G4Material* mat_W = G4Material::GetMaterial("G4_W");
	//G4Material* mat_Concrete = G4Material::GetMaterial("Concrete");
	//G4Material* mat_Pb = G4Material::GetMaterial("G4_Pb");
	G4Material* mat_MgO = G4Material::GetMaterial("G4_MAGNESIUM_OXIDE"); // Paint
	G4Material* mat_NaITl = G4Material::GetMaterial("NaITl"); // Crystal
	G4Material* mat_SiO2 = G4Material::GetMaterial("FusedSilica"); // Optical window
	G4Material* mat_BC630 = G4Material::GetMaterial("BC630"); // Optical grease
	G4Material* mat_BSiO2 = G4Material::GetMaterial("G4_Pyrex_Glass"); // PMT front glass
	G4Material* mat_BialkCat = G4Material::GetMaterial("BialkaliCathode"); // PMT cathode
	G4Material* mat_PP = G4Material::GetMaterial("G4_POLYPROPYLENE");
	G4Material* mat_Fe = G4Material::GetMaterial("G4_Fe");
	G4Material* mat_Stainless = G4Material::GetMaterial("G4_STAINLESS-STEEL"); // Frame
	G4Material* mat_DAW = G4Material::GetMaterial("DAW");

	// White painted
	G4double paint_Energies[2] = {.1*eV, 10.*eV};
	G4double paint_specularlobe[2] = {1.0, 1.0};		// for polished paint
	G4double paint_specularlobe_Lamb[2] = {0., 0.};	// for ground paint
	G4double paint_specularspike[2] = {0., 0.};
	G4double paint_backscatter[2] = {0., 0.};
	G4double Whitepainted_Reflectivity[2] = {0.93, 0.93};
	G4double Whitepainted_poor_Reflectivity[2] = {0.93, 0.93};

	//polished paint
	G4OpticalSurface* polishedWhitepainted = new G4OpticalSurface("polishedWhitepainted", unified);
	polishedWhitepainted->SetType(dielectric_dielectric);
	polishedWhitepainted->SetModel(unified);
	polishedWhitepainted->SetFinish(polishedfrontpainted);
	polishedWhitepainted->SetSigmaAlpha(1.3*degree);
	G4MaterialPropertiesTable* polishedWhitepainted_Prop = new G4MaterialPropertiesTable();
	polishedWhitepainted_Prop->AddProperty("SPECULARLOBECONSTANT", paint_Energies, paint_specularlobe, 2);
	polishedWhitepainted_Prop->AddProperty("SPECULARSPIKECONSTANT", paint_Energies, paint_specularspike, 2);
	polishedWhitepainted_Prop->AddProperty("BACKSCATTERCONSTANT", paint_Energies, paint_backscatter, 2);
	polishedWhitepainted_Prop->AddProperty("REFLECTIVITY", paint_Energies, Whitepainted_Reflectivity, 2);
	polishedWhitepainted->SetMaterialPropertiesTable(polishedWhitepainted_Prop);

	//ground paint
	G4OpticalSurface* groundWhitepainted = new G4OpticalSurface("groundWhitepainted", unified);
	groundWhitepainted->SetType(dielectric_dielectric);
	groundWhitepainted->SetModel(unified);
	groundWhitepainted->SetFinish(groundfrontpainted);
	G4MaterialPropertiesTable* groundWhitepainted_Prop = new G4MaterialPropertiesTable();
	groundWhitepainted_Prop->AddProperty("SPECULARLOBECONSTANT", paint_Energies, paint_specularlobe_Lamb, 2);
	groundWhitepainted_Prop->AddProperty("SPECULARSPIKECONSTANT", paint_Energies, paint_specularspike, 2);
	groundWhitepainted_Prop->AddProperty("BACKSCATTERCONSTANT", paint_Energies, paint_backscatter, 2);
	groundWhitepainted_Prop->AddProperty("REFLECTIVITY", paint_Energies, Whitepainted_Reflectivity, 2);
	groundWhitepainted->SetMaterialPropertiesTable(groundWhitepainted_Prop);

	//ground paint poor reflectivity
	G4OpticalSurface* groundWhitepainted_poor = new G4OpticalSurface("groundWhitepainted_poor", unified);
	groundWhitepainted_poor->SetType(dielectric_dielectric);
	groundWhitepainted_poor->SetModel(unified);
	groundWhitepainted_poor->SetFinish(groundfrontpainted);
	G4MaterialPropertiesTable* groundWhitepainted_poor_Prop = new G4MaterialPropertiesTable();
	groundWhitepainted_poor_Prop->AddProperty("SPECULARLOBECONSTANT", paint_Energies, paint_specularlobe_Lamb, 2);
	groundWhitepainted_poor_Prop->AddProperty("SPECULARSPIKECONSTANT", paint_Energies, paint_specularspike, 2);
	groundWhitepainted_poor_Prop->AddProperty("BACKSCATTERCONSTANT", paint_Energies, paint_backscatter, 2);
	groundWhitepainted_poor_Prop->AddProperty("REFLECTIVITY", paint_Energies, Whitepainted_poor_Reflectivity, 2);
	groundWhitepainted_poor->SetMaterialPropertiesTable(groundWhitepainted_poor_Prop);

	G4OpticalSurface* CathodeSurface = new G4OpticalSurface("CathodeSurface", glisur, ground, dielectric_metal, 1.);
	G4MaterialPropertiesTable* CathodeSurface_Prop = new G4MaterialPropertiesTable();
	G4double CathodeSurface_Energies[48] = {
			1.75*eV,
			1.75694888353765*eV,
			1.75694888353765*eV,
			1.75234552401747*eV,
			1.75694888353765*eV,
			1.79949383408072*eV,
			1.79466513864043*eV,
			1.81413709312839*eV,
			1.81907128286491*eV,
			1.83907939963336*eV,
			1.87517348130841*eV,
			1.87517348130841*eV,
			1.88044575913777*eV,
			1.89644198960302*eV,
			1.89644198960302*eV,
			1.92371584372004*eV,
			1.95178562743191*eV,
			1.97483821358268*eV,
			2.01045653807615*eV,
			2.06636006694130*eV,
			2.07920790155440*eV,
			2.12546146716102*eV,
			2.12546146716102*eV,
			2.12546146716102*eV,
			2.18090828804348*eV,
			2.18090828804348*eV,
			2.20973086453744*eV,
			2.23185275305895*eV,
			2.29306928571429*eV,
			2.40868622448980*eV,
			2.44389235688185*eV,
			2.50804453125000*eV,
			2.55596894904459*eV,
			2.66813248005319*eV,
			2.75609289148352*eV,
			2.97690745548961*eV,
			3.23618649193548*eV,
			3.36650272651007*eV,
			3.48945326086956*eV,
			3.76441955909944*eV,
			3.96528779644269*eV,
			4.01287125000000*eV,
			4.18880088726513*eV,
			4.32421470905172*eV,
			4.35235493492408*eV,
			4.40974862637363*eV,
			4.56008096590909*eV,
			4.57*eV
	};
	G4double CathodeSurface_Reflectivity[48] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
			0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	G4double CathodeSurface_Efficiency[48] = {
			0.0,
			0.000131721940389569,
			0.000142510267030300,
			0.000187717288986668,
			0.000274629027335151,
			0.000366524123707963,
			0.000529232764955913,
			0.000794849490680038,
			0.00143448217718175,
			0.00221171592363128,
			0.00280087319502880,
			0.00415177317192709,
			0.00546879618263594,
			0.00683530400293710,
			0.00821353126306006,
			0.0112533558260077,
			0.0135224109785807,
			0.0162489839920498,
			0.0214034770079185,
			0.0293248950992785,
			0.0366524123707963,
			0.0423429027220019,
			0.0482792687743996,
			0.0557748930826887,
			0.0661474064123016,
			0.0805346741782019,
			0.0980512176783237,
			0.108902296226373,
			0.129154966501489,
			0.167907936756026,
			0.174648653695795,
			0.191448197616996,
			0.215443469003188,
			0.248892262746700,
			0.262303097295234,
			0.291331513051222,
			0.295179008812718,
			0.280087319502881,
			0.258884120928555,
			0.209863698317641,
			0.184059087644166,
			0.161427383658924,
			0.153174046370208,
			0.113274213165710,
			0.0826758973518709,
			0.0644342546470982,
			0.0434686986040111,
			0};
	CathodeSurface_Prop->AddProperty("REFLECTIVITY", CathodeSurface_Energies, CathodeSurface_Reflectivity, 48);
	CathodeSurface_Prop->AddProperty("EFFICIENCY", CathodeSurface_Energies, CathodeSurface_Efficiency, 48);
	CathodeSurface->SetMaterialPropertiesTable(CathodeSurface_Prop);

	// -- Geometries -- //
	// World

	G4VSolid* sol_World = new G4Box("World", 0.5*World_Size, 0.5*World_Size, 0.5*World_Size);
	G4LogicalVolume* lv_World = new G4LogicalVolume(sol_World, mat_Air, "World");
	G4VPhysicalVolume* pv_World = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0), lv_World, "World", 0, false, 0);

	// Phantom
	G4RotationMatrix* PT_RotMat = new G4RotationMatrix();
	PT_RotMat->rotateX(90.*deg);

	G4VSolid* sol_rotWorld = new G4Tubs("rot_World", 0, 29.32*cm, 90*cm,0, 360*deg);
	G4LogicalVolume* lv_rotWorld = new G4LogicalVolume(sol_rotWorld, mat_Air, "rot_World");
	pv_rotWorld = new G4PVPlacement(PT_RotMat, G4ThreeVector(406.5, 0,0), lv_rotWorld, "rot_World", lv_World, false, 100000);
	G4VSolid* sol_Phantom = new G4Tubs("Phantom", 0, Phantom_outerD1/2, Phantom_H/2, 0., 360.*deg);
	G4LogicalVolume* lv_Phantom = new G4LogicalVolume(sol_Phantom, mat_DAW, "Phantom");
	//G4LogicalVolume* lv_Phantom = new G4LogicalVolume(sol_Phantom, mat_Air, "Phantom");
	phantom_pv = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0), lv_Phantom, "Phantom", lv_rotWorld, false, 98);
	G4VSolid* sol_innerPhantom = new G4Tubs("inPhantom",0,Phantom_outerD2/2,Phantom_H/2,0,360.*deg);
	lv_innerPhantom_air = new G4LogicalVolume(sol_innerPhantom,mat_Air,"inPhantom");
	lv_innerPhantom = new G4LogicalVolume(sol_innerPhantom,mat_Air,"inPhantom");
	new G4PVPlacement(0, G4ThreeVector(0,0,0),lv_innerPhantom_air,"inPhantom",lv_Phantom,false,99);
	new G4PVPlacement(0, G4ThreeVector(phantomholeD,0,0),lv_innerPhantom, "inPhantom", lv_Phantom,false,99);


	// Scatterlv_World
	G4RotationMatrix* Detector_RotMat = new G4RotationMatrix();
	Detector_RotMat->rotateY(270.*deg);
	Detector_RotMat->rotateZ(90.*deg);

	G4VSolid* sol_Det1 = new G4Box("Det1", WholeDet1_W/2, WholeDet1_D/2, WholeDet1_H/2);
	G4LogicalVolume* lv_Det1 = new G4LogicalVolume(sol_Det1, mat_Air, "Det1");
	//new G4PVPlacement(Detector_RotMat, G4ThreeVector(0, 0,-WholeDet1_H/2), lv_Det1, "Det1", lv_World, false, 1);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector(-102.5*mm, 0, 200*mm), lv_Det1, "Det1", lv_World, false, 1);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector(-102.5*mm, 0, -200*mm), lv_Det1, "Det1", lv_World, false, 1);

    //Housing
	G4VSolid* sol_Housing1 = new G4Box("Housing1", Housing1_W/2, Housing1_V/2, Housing_T/2);
	G4LogicalVolume* lv_Housing1 = new G4LogicalVolume(sol_Housing1, mat_Al, "Housing1");
	new G4PVPlacement(0, G4ThreeVector(0, 0, WholeDet1_H/2-Housing_T/2), lv_Housing1, "Housing1", lv_Det1, false, 8);


	//paint
	G4RotationMatrix* fillet_rot = new G4RotationMatrix();
	fillet_rot->rotateZ(45.*deg);
	G4VSolid* sol_Paint_box_main = new G4Box("Paintbox_m", Crystal_W/2+Paint_T, Crystal_W/2+Paint_T, (Paint_T+Crystal_T)/2);
	G4VSolid* sol_Paint1_fillet = new G4Box("Paint1_f", Crystal_Wf/2+Paint_T, Crystal_Wf/2+Paint_T, (Paint_T+Crystal_T)/2);
	G4VSolid* sol_Paint_box = new G4IntersectionSolid("Paint_box", sol_Paint_box_main, sol_Paint1_fillet, fillet_rot, G4ThreeVector());
	G4LogicalVolume* lv_Paint1_1 = new G4LogicalVolume(sol_Paint_box, mat_MgO, "Paint_box1");
	G4LogicalVolume* lv_Paint1_2 = new G4LogicalVolume(sol_Paint_box, mat_MgO, "Paint_box1");
	G4LogicalVolume* lv_Paint1_3 = new G4LogicalVolume(sol_Paint_box, mat_MgO, "Paint_box1");
	G4LogicalVolume* lv_Paint1_4 = new G4LogicalVolume(sol_Paint_box, mat_MgO, "Paint_box1");

	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2),(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), Housing_T/2-1.5*mm-(Paint_T+Crystal_T)/2), lv_Paint1_1, "Paint_box1", lv_Housing1, false, 81);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2),-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), Housing_T/2-1.5*mm-(Paint_T+Crystal_T)/2),lv_Paint1_2, "Paint_box2", lv_Housing1, false, 81);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2),(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), Housing_T/2-1.5*mm-(Paint_T+Crystal_T)/2), lv_Paint1_3, "Paint_box3", lv_Housing1, false, 81);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2),-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), Housing_T/2-1.5*mm-(Paint_T+Crystal_T)/2), lv_Paint1_4, "Paint_box4", lv_Housing1, false, 81);

	//Crystal1
	G4VSolid* sol_Crystal1_main = new G4Box("Crystal1_m", Crystal_W/2, Crystal_W/2, Crystal_T/2);
	G4VSolid* sol_Crystal1_fillet = new G4Box("Crystal1_f", Crystal_Wf/2, Crystal_Wf/2, Crystal_T/2);
	G4VSolid* sol_Crystal1 = new G4IntersectionSolid("Crystal1", sol_Crystal1_main, sol_Crystal1_fillet, fillet_rot, G4ThreeVector());
	G4LogicalVolume* lv_Crystal1 = new G4LogicalVolume(sol_Crystal1, mat_NaITl, "Crystal1");

	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint_T/2), lv_Crystal1, "Crystal1", lv_Paint1_1, false, 401);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint_T/2), lv_Crystal1, "Crystal1", lv_Paint1_2, false, 402);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint_T/2), lv_Crystal1, "Crystal1", lv_Paint1_3, false, 403);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint_T/2), lv_Crystal1, "Crystal1", lv_Paint1_4, false, 404);

	//OpticalGlue1
	G4VSolid* sol_OptGlue1_main = new G4Box("OptGlue1_m", Crystal_W/2+Paint_T, Crystal_W/2+Paint_T, OptGlue_T/2);
	G4VSolid* sol_OptGlue1_fillet = new G4Box("OptGlue1_f", Crystal_Wf/2+Paint_T, Crystal_Wf/2+Paint_T, OptGlue_T/2);
	G4VSolid* sol_OptGlue1 = new G4IntersectionSolid("OptGlue1", sol_OptGlue1_main, sol_OptGlue1_fillet, fillet_rot, G4ThreeVector());
	G4LogicalVolume* lv_OptGlue1 = new G4LogicalVolume(sol_OptGlue1, mat_BC630, "OptGlue1_1");

	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), (CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -Housing_T/2+OptGrease_T+OptWin1_T + OptGlue_T/2), lv_OptGlue1, "OptGlue1_1", lv_Housing1, false, 82);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), (CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -Housing_T/2+OptGrease_T+OptWin1_T+ OptGlue_T/2), lv_OptGlue1, "OptGlue1_2", lv_Housing1, false, 82);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -Housing_T/2+OptGrease_T+OptWin1_T+ OptGlue_T/2), lv_OptGlue1, "OptGlue1_3", lv_Housing1, false, 82);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -(CrossBar_W/2+3.5*mm+Paint_T+Crystal_W/2), -Housing_T/2+OptGrease_T+OptWin1_T+ OptGlue_T/2), lv_OptGlue1, "OptGlue1_4", lv_Housing1, false, 82);

	//OptWin1
	G4VSolid* sol_OptWin1 = new G4Box("OptWin1", OptWin_W/2, OptWin_W/2, OptWin1_T/2);
	G4LogicalVolume* lv_OptWin1 = new G4LogicalVolume(sol_OptWin1, mat_SiO2, "OptWin1_1");
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+OptWin_W)/2, (CrossBar_W+OptWin_W)/2, -Housing_T/2+OptGrease_T+OptWin1_T/2), lv_OptWin1, "OptWin1_1", lv_Housing1, false, 83);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin_W)/2, (CrossBar_W+OptWin_W)/2, -Housing_T/2+OptGrease_T+OptWin1_T/2), lv_OptWin1, "OptWin1_2", lv_Housing1, false, 83);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin_W)/2, -(CrossBar_W+OptWin_W)/2, -Housing_T/2+OptGrease_T+OptWin1_T/2), lv_OptWin1, "OptWin1_3", lv_Housing1, false, 83);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+OptWin_W)/2, -(CrossBar_W+OptWin_W)/2, -Housing_T/2+OptGrease_T+OptWin1_T/2), lv_OptWin1, "OptWin1_4", lv_Housing1, false, 83);

	//OptGrease1
	G4VSolid* sol_OptGrease1 = new G4Box("OptGrease1_1", OptWin_W/2, OptWin_W/2, OptGrease_T/2);
	G4LogicalVolume* lv_OptGrease1 = new G4LogicalVolume(sol_OptGrease1, mat_BC630, "OptGrease1_1");
	new G4PVPlacement(0, G4ThreeVector( (CrossBar_W+OptWin_W)/2, (CrossBar_W+OptWin_W)/2, -(Housing_T/2)+ OptGrease_T/2), lv_OptGrease1, "OptGrease1_1", lv_Housing1, false, 84);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin_W)/2, (CrossBar_W+OptWin_W)/2, -Housing_T/2+ OptGrease_T/2), lv_OptGrease1, "OptGrease1_2", lv_Housing1, false, 84);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin_W)/2,-(CrossBar_W+OptWin_W)/2, -Housing_T/2+ OptGrease_T/2), lv_OptGrease1, "OptGrease1_3", lv_Housing1, false, 84);
	new G4PVPlacement(0, G4ThreeVector( (CrossBar_W+OptWin_W)/2,-(CrossBar_W+OptWin_W)/2, -Housing_T/2+ OptGrease_T/2), lv_OptGrease1, "OptGrease1_4", lv_Housing1, false, 84);

	//RearCase1_Shell
	G4VSolid* sol_RearCase1_Shell = new G4Box("RearCase1_Shell", RearCase_W/2, RearCase1_D/2, RearCase1_H/2);
	G4LogicalVolume* lv_RearCase1_Shell = new G4LogicalVolume(sol_RearCase1_Shell, mat_Al, "RearCase1_Shell");
	new G4PVPlacement(0, G4ThreeVector(0, 0, -WholeDet1_H/2+RearCase1_H/2+rear_Al_T), lv_RearCase1_Shell, "RearCase1_Shell", lv_Det1, false, 9);

	G4VSolid* sol_RearCase1 = new G4Box("RearCase1", RearCase_W/2-RearCase_T, RearCase1_D/2-RearCase_T, (RearCase1_H-RearCase_T)/2);
	G4LogicalVolume* lv_RearCase1 = new G4LogicalVolume(sol_RearCase1, mat_Air, "RearCase1");
	new G4PVPlacement(0, G4ThreeVector(0, 0, RearCase_T/2), lv_RearCase1, "RearCase1", lv_RearCase1_Shell, false, 91);

	//CrossBar
	G4VSolid* sol_CrossBar1 = new G4Box("CrossBar1", CrossBar_W/2, CrossBar_H/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar1 = new G4LogicalVolume(sol_CrossBar1, mat_BC630, "CrossBar1");
	new G4PVPlacement(0, G4ThreeVector(0, 0, -CrossBar_T/2 + RearCase1_H/2), lv_CrossBar1, "CrossBar1", lv_RearCase1_Shell, false, 1000);

	G4VSolid* sol_CrossBar1_H1 = new G4Box("CrossBar1_H1", (CrossBar_H-CrossBar_W)/4, CrossBar_W/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar1_H1 = new G4LogicalVolume(sol_CrossBar1_H1, mat_BC630, "CrossBar1_H1");
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+CrossBar_H)/4, 0, -CrossBar_T/2 + RearCase1_H/2), lv_CrossBar1_H1, "CrossBar1_H1", lv_RearCase1_Shell, false, 1000);

	G4VSolid* sol_CrossBar1_H2 = new G4Box("CrossBar1_H2", (CrossBar_H-CrossBar_W)/4, CrossBar_W/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar1_H2 = new G4LogicalVolume(sol_CrossBar1_H2, mat_BC630, "CrossBar1_H2");
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+CrossBar_H)/4, 0, -CrossBar_T/2 + RearCase1_H/2), lv_CrossBar1_H2, "CrossBar1_H2", lv_RearCase1_Shell, false, 1000);

	//PMT_case
	G4VSolid* sol_PMT9 = new G4Box("PMT_group", PMTbox/2,PMTbox/2, WholePMT_T/2);
	G4LogicalVolume* lv_PMT9_1 = new G4LogicalVolume(sol_PMT9, mat_Air, "PMT1_group");
	G4LogicalVolume* lv_PMT9_2 = new G4LogicalVolume(sol_PMT9, mat_Air, "PMT1_group");
	G4LogicalVolume* lv_PMT9_3 = new G4LogicalVolume(sol_PMT9, mat_Air, "PMT1_group");
	G4LogicalVolume* lv_PMT9_4 = new G4LogicalVolume(sol_PMT9, mat_Air, "PMT1_group");
	new G4PVPlacement(0, G4ThreeVector(((CrossBar_W/2)+PMTbox/2), ((CrossBar_W/2)+PMTbox/2), (RearCase1_H-RearCase_T)/2-WholePMT_T/2), lv_PMT9_1,"PMT_group_1",lv_RearCase1,false,800);
	new G4PVPlacement(0, G4ThreeVector(-((CrossBar_W/2)+PMTbox/2), ((CrossBar_W/2)+PMTbox/2), (RearCase1_H-RearCase_T)/2-WholePMT_T/2), lv_PMT9_2,"PMT_group_2",lv_RearCase1,false,801);
	new G4PVPlacement(0, G4ThreeVector(-((CrossBar_W/2)+PMTbox/2), -((CrossBar_W/2)+PMTbox/2), (RearCase1_H-RearCase_T)/2-WholePMT_T/2), lv_PMT9_3,"PMT_group_3",lv_RearCase1,false,802);
	new G4PVPlacement(0, G4ThreeVector(((CrossBar_W/2)+PMTbox/2), -((CrossBar_W/2)+PMTbox/2), (RearCase1_H-RearCase_T)/2-WholePMT_T/2), lv_PMT9_4,"PMT_group_4",lv_RearCase1,false,803);

	//nine PMT arrangement
	G4VSolid* sol_PMT1 = new G4Box("PMT1", WholePMT_W/2, WholePMT_W/2, WholePMT_T/2);
	G4LogicalVolume* lv_PMT1 = new G4LogicalVolume(sol_PMT1, mat_Air, "PMT1");

	std::vector < G4double > xPos1;
	std::vector < G4double > yPos1;
	for(G4int i=0;i<9;++i){
		if(i<3)	{
			yPos1.push_back(PMTFrontBase_W);
		}
		else if (i > 2 && i < 6){
			yPos1.push_back(0);
		}
		else{
			yPos1.push_back(-PMTFrontBase_W);
		}
	}

	for(G4int i=0;i<9;i++){
		if(i%3==0)	{
			xPos1.push_back(-PMTFrontBase_W);}
		else if ( i % 3 == 1) {
			xPos1.push_back(0);	}
		else{
			xPos1.push_back(PMTFrontBase_W);}
	}

	for(G4int i=0; i<9; i++){
			new G4PVPlacement(0, G4ThreeVector(xPos1[i], yPos1[i], 0), lv_PMT1,"PMT1",lv_PMT9_1,false,101+i);
			new G4PVPlacement(0, G4ThreeVector(xPos1[i], yPos1[i], 0), lv_PMT1,"PMT1",lv_PMT9_2,false,101+9+i);
			new G4PVPlacement(0, G4ThreeVector(xPos1[i], yPos1[i], 0), lv_PMT1,"PMT1",lv_PMT9_3,false,101+9*2+i);
			new G4PVPlacement(0, G4ThreeVector(xPos1[i], yPos1[i], 0), lv_PMT1,"PMT1",lv_PMT9_4,false,101+9*3+i);
	}

	//PMT_inner_components

		G4VSolid* sol_PMTFrontBase1_Shell = new G4Box("PMTFrontBase1_Shell", PMTFrontBase_W/2, PMTFrontBase_W/2, PMTFrontBase_T/2);
		G4LogicalVolume* lv_PMTFrontBase1_Shell = new G4LogicalVolume(sol_PMTFrontBase1_Shell, mat_BSiO2, "PMTFrontBase1_Shell");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT_T/2-PMTFrontBase_T/2), lv_PMTFrontBase1_Shell, "PMTFrontBase1_Shell", lv_PMT1, false, 92);

		G4VSolid* sol_PMTFrontBase1 = new G4Box("PMTFrontBase1", PMTFrontBase_W/2-PMTFrontGlass_T, PMTFrontBase_W/2-PMTFrontGlass_T, 	(PMTFrontBase_T-PMTFrontGlass_T)/2);
		G4LogicalVolume* lv_PMTFrontBase1 = new G4LogicalVolume(sol_PMTFrontBase1, mat_Air, "PMTFrontBase1");
		new G4PVPlacement(0, G4ThreeVector(0, 0, -PMTFrontGlass_T/2), lv_PMTFrontBase1, "PMTFrontBase1", lv_PMTFrontBase1_Shell, false, 93);

		G4VSolid* sol_PMTPhotoCathode1 = new G4Box("PMTPhotoCathode1", PMTPhotoCathode_W/2, PMTPhotoCathode_W/2, PMTPhotoCathode_T/2);
		G4LogicalVolume* lv_PMTPhotoCathode1 = new G4LogicalVolume(sol_PMTPhotoCathode1, mat_BialkCat, "PMTPhotoCathode1");
		new G4PVPlacement(0, G4ThreeVector(0, 0, (PMTFrontBase_T-PMTFrontGlass_T)/2-PMTPhotoCathode_T/2), lv_PMTPhotoCathode1, "PMTPhotoCathode1", 	lv_PMTFrontBase1, false, 100);

		G4VSolid* sol_PMTRearBase1_Shell = new G4Tubs("PMTRearBase1_Shell", 0, PMTRearBase_D/2, PMTRearBase_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTRearBase1_Shell = new G4LogicalVolume(sol_PMTRearBase1_Shell, mat_BSiO2, "PMTRearBase1_Shell");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT_T/2-PMTFrontBase_T-PMTRearBase_H/2), lv_PMTRearBase1_Shell, "PMTRearBase1_Shell", 	lv_PMT1, false, 94);

		G4VSolid* sol_PMTTail1 = new G4Tubs("PMTTail1", 0, PMTTail_D/2, PMTTail_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTTail1 = new G4LogicalVolume(sol_PMTTail1, mat_PE, "PMTTail1");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT_T/2-PMTFrontBase_T-PMTRearBase_H-PMTTail_H/2), lv_PMTTail1, "PMTTail1", lv_PMT1, false, 96);

		G4VSolid* sol_PMTTail2 = new G4Tubs("PMTTail2", 0, PMTTail_D/2, PMTTail2_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTTail2 = new G4LogicalVolume(sol_PMTTail2, mat_PE, "PMTTail2");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT_T/2-PMTFrontBase_T-PMTRearBase_H-PMTTail_H-PMTPCB_T-PMTTail2_H/2), lv_PMTTail2, "PMTTail2", lv_PMT1, false, 96);

		G4VSolid* sol_PMTPCB1 = new G4Tubs("PMTPCB1", 0, PMTPCB_D/2, PMTPCB_T/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTPCB1 = new G4LogicalVolume(sol_PMTPCB1, mat_PE, "PMTPCB1");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT_T/2-PMTFrontBase_T-PMTRearBase_H-PMTTail_H-PMTPCB_T/2), lv_PMTPCB1, "PMTPCB1", lv_PMT1, false, 97);


	//Pillar
	G4VSolid* sol_Pillar1 = new G4Tubs("Pillar1", 0, Pillar1_D/2, Pillar1_H/2, 0, 360.*deg);
	G4LogicalVolume* lv_Pillar1 = new G4LogicalVolume(sol_Pillar1, mat_Al, "Pillar1");

	std::vector < G4double > yPos_pillar;

	for(G4int i=1; i<7;++i){
		if(i<4){
			yPos_pillar.push_back(CrossBar_W/2 + PMTbox - 0.5 * PMTFrontBase_W *(2*i-1) );
		}
		else {
			yPos_pillar.push_back(-(CrossBar_W/2 + 0.5 * PMTFrontBase_W *(2*i-7)) );
		}
	}

	for(G4int i=0; i<6;++i){
		new G4PVPlacement(0,G4ThreeVector((CrossBar_W/2+PMTbox+21.37*mm+Pillar1_D/2) ,yPos_pillar[i], (RearCase1_H-RearCase_T)/2-Pillar1_H/2 ), lv_Pillar1, "Pillar1", lv_RearCase1, false, 98);
		new G4PVPlacement(0,G4ThreeVector(-(CrossBar_W/2+PMTbox+21.37*mm+Pillar1_D/2),yPos_pillar[i], (RearCase1_H-RearCase_T)/2-Pillar1_H/2 ), lv_Pillar1, "Pillar1", lv_RearCase1, false, 98);
	}

	// Holder
	G4VSolid* sol_Holder1 = new G4Box("Holder1", Holder1_W/2, Holder1_H/2, Holder1_T/2);
	G4LogicalVolume* lv_Holder1 = new G4LogicalVolume(sol_Holder1, mat_Al, "Holder1");
	for(G4int i=1; i<7; i++){
			if(i<4){
	new G4PVPlacement(0, G4ThreeVector( 0, (CrossBar_W+PMTFrontBase_W*(2*i-1))/2, (RearCase1_H-RearCase_T)/2-WholePMT_T-Holder1_T/2), lv_Holder1, "Holder1", lv_RearCase1, false, 98);
			}
			else {
	new G4PVPlacement(0, G4ThreeVector(0,  -(CrossBar_W+PMTFrontBase_W*(2*i-7) )/2, (RearCase1_H-RearCase_T)/2-WholePMT_T-Holder1_T/2), lv_Holder1, "Holder1", lv_RearCase1, false, 98);}
		}

	G4VSolid* sol_rear_Al = new G4Box("rear_Al", rear_Al_W/2, rear_Al_H/2, rear_Al_T/2);
	G4LogicalVolume* lv_rear_Al = new G4LogicalVolume(sol_rear_Al, mat_Al, "rear_Al");
	new G4PVPlacement(0, G4ThreeVector(0,0,-WholeDet1_H/2+rear_Al_T/2), lv_rear_Al, "rear_Al", lv_Det1, false, 700);

	new G4LogicalSkinSurface("alSurface1",lv_Housing1, groundWhitepainted);
	new G4LogicalSkinSurface("pmtSurface1",lv_PMTPhotoCathode1, CathodeSurface);

	// Absorber - type2: same as scatter but different size
	G4VSolid* sol_Det11 = new G4Box("Det11", WholeDet11_W/2, WholeDet11_D/2, WholeDet11_H/2);
	G4LogicalVolume* lv_Det11 = new G4LogicalVolume(sol_Det11, mat_Air, "Det11");
	//new G4PVPlacement(Detector_RotMat, G4ThreeVector(0, 0, -WholeDet1_H-Sc2Ab_Dist-WholeDet11_H/2), lv_Det11, "Det11", lv_World, false, 2);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector(-357.5*mm, 0,200*mm), lv_Det11, "Det11", lv_World, false, 2);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector(-357.5*mm, 0,-200*mm), lv_Det11, "Det11", lv_World, false, 2);
	
	//Housing
	G4VSolid* sol_Housing11 = new G4Box("Housing11", Housing11_W/2, Housing11_V/2, Housing11_T/2);
	G4LogicalVolume* lv_Housing11 = new G4LogicalVolume(sol_Housing11, mat_Al, "Housing11");
	new G4PVPlacement(0, G4ThreeVector(0, 0, WholeDet11_H/2-Housing11_T/2), lv_Housing11, "Housing11", lv_Det11, false, 8);

	//Paint
	G4VSolid* sol_Paint_box11_main = new G4Box("Paintbox11_m", Crystal11_W/2+Paint11_T, Crystal11_W/2+Paint11_T, (Paint11_T+Crystal11_T)/2);
	G4VSolid* sol_Paint11_fillet = new G4Box("Paint11_f", Crystal11_Wf/2+Paint11_T, Crystal11_Wf/2+Paint11_T, (Paint11_T+Crystal11_T)/2);
	G4VSolid* sol_Paint_box11 = new G4IntersectionSolid("Paint_box11", sol_Paint_box11_main, sol_Paint11_fillet, fillet_rot, G4ThreeVector());

	G4LogicalVolume* lv_Paint11_1 = new G4LogicalVolume(sol_Paint_box11, mat_MgO, "Paint_box11_1");
	G4LogicalVolume* lv_Paint11_2 = new G4LogicalVolume(sol_Paint_box11, mat_MgO, "Paint_box11_1");
	G4LogicalVolume* lv_Paint11_3 = new G4LogicalVolume(sol_Paint_box11, mat_MgO, "Paint_box11_1");
	G4LogicalVolume* lv_Paint11_4 = new G4LogicalVolume(sol_Paint_box11, mat_MgO, "Paint_box11_1");

	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2),(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), Housing11_T/2-1.5*mm-(Paint11_T+Crystal11_T)/2), lv_Paint11_3, "Paint_box11_3", lv_Housing11, false, 281);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2),(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), Housing11_T/2-1.5*mm-(Paint11_T+Crystal11_T)/2), lv_Paint11_1, "Paint_box11_1", lv_Housing11, false, 281);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2),-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), Housing11_T/2-1.5*mm-(Paint11_T+Crystal11_T)/2), lv_Paint11_2, "Paint_box11_2", lv_Housing11, false, 281);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2),-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), Housing11_T/2-1.5*mm-(Paint11_T+Crystal11_T)/2), lv_Paint11_4, "Paint_box11_4", lv_Housing11, false, 281);


	//Crystal11,
	G4VSolid* sol_Crystal11_main = new G4Box("Crystal11_m", Crystal11_W/2, Crystal11_W/2, Crystal11_T/2);
	G4VSolid* sol_Crystal11_fillet = new G4Box("Crystal11_f", Crystal11_Wf/2, Crystal11_Wf/2, Crystal11_T/2);
	G4VSolid* sol_Crystal11 = new G4IntersectionSolid("Crystal11_1", sol_Crystal11_main, sol_Crystal11_fillet, fillet_rot, G4ThreeVector());
	G4LogicalVolume* lv_Crystal11 = new G4LogicalVolume(sol_Crystal11, mat_NaITl, "Crystal11");

	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint11_T/2), lv_Crystal11, "Crystal11", lv_Paint11_1, false, 405);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint11_T/2), lv_Crystal11, "Crystal11", lv_Paint11_2, false, 406);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint11_T/2), lv_Crystal11, "Crystal11", lv_Paint11_3, false, 407);
	new G4PVPlacement(0, G4ThreeVector(0,0, -Paint11_T/2), lv_Crystal11, "Crystal11", lv_Paint11_4, false, 408);

	//OpticalGlue11,

	G4VSolid* sol_OptGlue11_main = new G4Box("OptGlue11_m", Crystal11_W/2+Paint11_T, Crystal11_W/2+Paint11_T, OptGlue11_T/2);
	G4VSolid* sol_OptGlue11_fillet = new G4Box("OptGlue11_f", Crystal11_Wf/2+Paint11_T, Crystal11_Wf/2+Paint11_T, OptGlue11_T/2);
	G4VSolid* sol_OptGlue11 = new G4IntersectionSolid("OptGlue11", sol_OptGlue11_main, sol_OptGlue11_fillet, fillet_rot, G4ThreeVector());
	G4LogicalVolume* lv_OptGlue11 = new G4LogicalVolume(sol_OptGlue11, mat_BC630, "OptGlue11_1");

	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), (CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -Housing11_T/2+OptGrease11_T+OptWin11_T +OptGlue11_T/2), lv_OptGlue11, "OptGlue11_1", lv_Housing11, false, 82);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), (CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -Housing11_T/2+OptGrease11_T+OptWin11_T	+OptGlue11_T/2), lv_OptGlue11, "OptGlue11_2", lv_Housing11, false, 82);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -Housing11_T/2+OptGrease11_T+OptWin11_T	+OptGlue11_T/2), lv_OptGlue11, "OptGlue11_3", lv_Housing11, false, 82);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -(CrossBar_W/2+3.5*mm+Paint11_T+Crystal11_W/2), -Housing11_T/2+OptGrease11_T+OptWin11_T	+OptGlue11_T/2), lv_OptGlue11, "OptGlue11_4", lv_Housing11, false, 82);

	//OptWin11

	G4VSolid* sol_OptWin11 = new G4Box("OptWin11", OptWin11_W/2, OptWin11_W/2, OptWin11_T/2);
	G4LogicalVolume* lv_OptWin11 = new G4LogicalVolume(sol_OptWin11, mat_SiO2, "OptWin11_1");
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+OptWin11_W)/2, (CrossBar_W+OptWin11_W)/2, -Housing11_T/2+OptGrease11_T+OptWin11_T/2), lv_OptWin11, "OptWin11_1", lv_Housing11, false, 83);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin11_W)/2, (CrossBar_W+OptWin11_W)/2, -Housing11_T/2+OptGrease11_T+OptWin11_T/2), lv_OptWin11, "OptWin11_2", lv_Housing11, false, 83);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin11_W)/2, -(CrossBar_W+OptWin11_W)/2, -Housing11_T/2+OptGrease11_T+OptWin11_T/2), lv_OptWin11, "OptWin11_3", lv_Housing11, false, 83);
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+OptWin11_W)/2, -(CrossBar_W+OptWin11_W)/2, -Housing11_T/2+OptGrease11_T+OptWin11_T/2), lv_OptWin11, "OptWin11_4", lv_Housing11, false, 83);


	//OptGrease11
	G4VSolid* sol_OptGrease11 = new G4Box("OptGrease11", OptWin11_W/2, OptWin11_W/2, OptGrease11_T/2);
	G4LogicalVolume* lv_OptGrease11 = new G4LogicalVolume(sol_OptGrease11, mat_BC630, "OptGrease11_1");
	new G4PVPlacement(0, G4ThreeVector( (CrossBar_W+OptWin11_W)/2, (CrossBar_W+OptWin11_W)/2, -Housing11_T/2+ OptGrease11_T/2), lv_OptGrease11, "OptGrease11_1", lv_Housing11, false, 84);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin11_W)/2, (CrossBar_W+OptWin11_W)/2, -Housing11_T/2+ OptGrease11_T/2), lv_OptGrease11, "OptGrease11_2", lv_Housing11, false, 84);
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+OptWin11_W)/2,-(CrossBar_W+OptWin11_W)/2, -Housing11_T/2+ OptGrease11_T/2), lv_OptGrease11, "OptGrease11_3", lv_Housing11, false, 84);
	new G4PVPlacement(0, G4ThreeVector( (CrossBar_W+OptWin11_W)/2,-(CrossBar_W+OptWin11_W)/2, -Housing11_T/2+ OptGrease11_T/2), lv_OptGrease11, "OptGrease11_4", lv_Housing11, false, 84);

	//RearCase11_Shell
	G4VSolid* sol_RearCase11_Shell = new G4Box("RearCase11_Shell", RearCase11_W/2, RearCase11_D/2, RearCase11_H/2);
	G4LogicalVolume* lv_RearCase11_Shell = new G4LogicalVolume(sol_RearCase11_Shell, mat_Al, "RearCase11_Shell");
	new G4PVPlacement(0, G4ThreeVector(0, 0, -WholeDet11_H/2+RearCase11_H/2+rear_Al_T), lv_RearCase11_Shell, "RearCase11_Shell", lv_Det11, false, 9);

	G4VSolid* sol_RearCase11 = new G4Box("RearCase11", RearCase11_W/2-RearCase11_T, RearCase11_D/2-RearCase11_T, (RearCase11_H-RearCase11_T)/2);
	G4LogicalVolume* lv_RearCase11 = new G4LogicalVolume(sol_RearCase11, mat_Air, "RearCase11");
	new G4PVPlacement(0, G4ThreeVector(0, 0, RearCase11_T/2), lv_RearCase11, "RearCase11", lv_RearCase11_Shell, false, 91);

	//CrossBar
	G4VSolid* sol_CrossBar11 = new G4Box("CrossBar11", CrossBar_W/2, CrossBar_H/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar11 = new G4LogicalVolume(sol_CrossBar11, mat_BC630, "CrossBar11");
	new G4PVPlacement(0, G4ThreeVector(0, 0, -CrossBar_T/2 + RearCase11_H/2), lv_CrossBar11, "CrossBar11", lv_RearCase11_Shell, false, 1000);

	G4VSolid* sol_CrossBar11_H1 = new G4Box("CrossBar11_H1", (CrossBar_H-CrossBar_W)/4, CrossBar_W/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar11_H1 = new G4LogicalVolume(sol_CrossBar11_H1, mat_BC630, "CrossBar11_H1");
	new G4PVPlacement(0, G4ThreeVector(-(CrossBar_W+CrossBar_H)/4, 0, -CrossBar_T/2 + RearCase11_H/2), lv_CrossBar11_H1, "CrossBar11_H1", lv_RearCase11_Shell, false, 1000);

	G4VSolid* sol_CrossBar11_H2 = new G4Box("CrossBar11_H2", (CrossBar_H-CrossBar_W)/4, CrossBar_W/2, CrossBar_T/2);
	G4LogicalVolume* lv_CrossBar11_H2 = new G4LogicalVolume(sol_CrossBar11_H2, mat_BC630, "CrossBar11_H2");
	new G4PVPlacement(0, G4ThreeVector((CrossBar_W+CrossBar_H)/4, 0, -CrossBar_T/2 + RearCase11_H/2), lv_CrossBar11_H2, "CrossBar11_H2", lv_RearCase11_Shell, false, 1000);

	//PMT_case
	G4VSolid* sol_PMT9_11 = new G4Box("PMT_group", PMTbox11/2,PMTbox11/2, WholePMT11_H/2);
	G4LogicalVolume* lv_PMT9_11_1 = new G4LogicalVolume(sol_PMT9_11, mat_Air, "PMT11_group");
	G4LogicalVolume* lv_PMT9_11_2 = new G4LogicalVolume(sol_PMT9_11, mat_Air, "PMT11_group");
	G4LogicalVolume* lv_PMT9_11_3 = new G4LogicalVolume(sol_PMT9_11, mat_Air, "PMT11_group");
	G4LogicalVolume* lv_PMT9_11_4 = new G4LogicalVolume(sol_PMT9_11, mat_Air, "PMT11_group");

	new G4PVPlacement(0, G4ThreeVector(((CrossBar_W/2)+PMTbox11/2), ((CrossBar_W/2)+PMTbox11/2), (RearCase11_H-RearCase11_T)/2-WholePMT11_H/2),  lv_PMT9_11_1,"PMT_group_1",lv_RearCase11,false,800);
	new G4PVPlacement(0, G4ThreeVector(-((CrossBar_W/2)+PMTbox11/2), ((CrossBar_W/2)+PMTbox11/2), (RearCase11_H-RearCase11_T)/2-WholePMT11_H/2), lv_PMT9_11_2,"PMT_group_2",lv_RearCase11,false,801);
	new G4PVPlacement(0, G4ThreeVector(-((CrossBar_W/2)+PMTbox11/2), -((CrossBar_W/2)+PMTbox11/2), (RearCase11_H-RearCase11_T)/2-WholePMT11_H/2), lv_PMT9_11_3,"PMT_group_3",lv_RearCase11,false,802);
	new G4PVPlacement(0, G4ThreeVector(((CrossBar_W/2)+PMTbox11/2), -((CrossBar_W/2)+PMTbox11/2), (RearCase11_H-RearCase11_T)/2-WholePMT11_H/2), lv_PMT9_11_4,"PMT_group_4",lv_RearCase11,false,803);

	//nine PMT arrangement

	G4VSolid* sol_PMT11 = new G4Box("PMT11", WholePMT11_W/2, WholePMT11_W/2, WholePMT11_H/2);
	G4LogicalVolume* lv_PMT11 = new G4LogicalVolume(sol_PMT11, mat_Air, "PMT11");

	std::vector < G4double > xPos11;
	std::vector < G4double > yPos11;


	for(G4int i=0;i<9;++i){
		if(i<3)	{
			yPos11.push_back(PMTFrontBase11_W);
		}
		else if (i > 2 && i < 6){
			yPos11.push_back(0);
		}
		else{
			yPos11.push_back(-PMTFrontBase11_W);
		}
	}
	for(G4int i=0;i<9;i++){
		if(i%3==0)	{
			xPos11.push_back(-PMTFrontBase11_W);}
		else if ( i % 3 == 1) {
			xPos11.push_back(0);	}
		else{
			xPos11.push_back(PMTFrontBase11_W);}
	}

	for(G4int i=0; i<9; i++){
		new G4PVPlacement(0, G4ThreeVector(xPos11[i], yPos11[i], 0), lv_PMT11,"PMT11",lv_PMT9_11_1,false,101+i);
		new G4PVPlacement(0, G4ThreeVector(xPos11[i], yPos11[i], 0), lv_PMT11,"PMT11",lv_PMT9_11_2,false,101+9+i);
		new G4PVPlacement(0, G4ThreeVector(xPos11[i], yPos11[i], 0), lv_PMT11,"PMT11",lv_PMT9_11_3,false,101+9*2+i);
		new G4PVPlacement(0, G4ThreeVector(xPos11[i], yPos11[i], 0), lv_PMT11,"PMT11",lv_PMT9_11_4,false,101+9*3+i);
	}

	//PMT11_inner_components

		G4VSolid* sol_PMTFrontBase11_Shell = new G4Box("PMTFrontBase11_Shell", PMTFrontBase11_W/2, PMTFrontBase11_W/2, PMTFrontBase11_H/2);
		G4LogicalVolume* lv_PMTFrontBase11_Shell = new G4LogicalVolume(sol_PMTFrontBase11_Shell, mat_BSiO2, "PMTFrontBase11_Shell");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT11_H/2-PMTFrontBase11_H/2), lv_PMTFrontBase11_Shell, "PMTFrontBase11_Shell", lv_PMT11, false, 92);

		G4VSolid* sol_PMTFrontBase11 = new G4Box("PMTFrontBase11", PMTFrontBase11_W/2-PMTFrontGlass11_T, PMTFrontBase11_W/2-PMTFrontGlass11_T, (PMTFrontBase11_H-PMTFrontGlass11_T)/2);
		G4LogicalVolume* lv_PMTFrontBase11 = new G4LogicalVolume(sol_PMTFrontBase11, mat_Air, "PMTFrontBase11");
		new G4PVPlacement(0, G4ThreeVector(0, 0, -PMTFrontGlass11_T/2), lv_PMTFrontBase11, "PMTFrontBase11", lv_PMTFrontBase11_Shell, false, 93);

		G4VSolid* sol_PMTPhotoCathode11 = new G4Box("PMTPhotoCathode11", PMTPhotoCathode11_W/2, PMTPhotoCathode11_W/2, PMTPhotoCathode11_T/2);
		G4LogicalVolume* lv_PMTPhotoCathode11 = new G4LogicalVolume(sol_PMTPhotoCathode11, mat_BialkCat, "PMTPhotoCathode11");
		new G4PVPlacement(0, G4ThreeVector(0, 0, (PMTFrontBase11_H-PMTFrontGlass11_T)/2-PMTPhotoCathode11_T/2), lv_PMTPhotoCathode11, "PMTPhotoCathode11", lv_PMTFrontBase11, false, 100);

		G4VSolid* sol_PMTRearBase11_Shell = new G4Tubs("PMTRearBase11_Shell", 0, PMTRearBase11_D/2, PMTRearBase11_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTRearBase11_Shell = new G4LogicalVolume(sol_PMTRearBase11_Shell, mat_BSiO2, "PMTRearBase11_Shell");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT11_H/2-PMTFrontBase11_H-PMTRearBase11_H/2), lv_PMTRearBase11_Shell, "PMTRearBase11_Shell", lv_PMT11, false, 94);

		G4VSolid* sol_PMTTail11 = new G4Tubs("PMTTail11", 0, PMTTail11_D/2, PMTTail11_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTTail11 = new G4LogicalVolume(sol_PMTTail11, mat_PE, "PMTTail11");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT11_H/2-PMTFrontBase11_H-PMTRearBase11_H-PMTTail11_H/2), lv_PMTTail11, "PMTTail11", lv_PMT11, false, 96);

		G4VSolid* sol_PMTTail11_2 = new G4Tubs("PMTTail11_2", 0, PMTTail11_D/2, PMTTail12_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTTail11_2 = new G4LogicalVolume(sol_PMTTail11_2, mat_PE, "PMTTail11_2");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT11_H/2-PMTFrontBase11_H-PMTRearBase11_H-PMTTail11_H-PMTPCB11_T-PMTTail12_H/2), lv_PMTTail11_2, "PMTTail11_2", lv_PMT11, false, 96);

		G4VSolid* sol_PMTPCB11 = new G4Tubs("PMTPCB11", 0, PMTPCB11_D/2, PMTPCB11_T/2, 0, 360.*deg);
		G4LogicalVolume* lv_PMTPCB11 = new G4LogicalVolume(sol_PMTPCB11, mat_PE, "PMTPCB11");
		new G4PVPlacement(0, G4ThreeVector(0, 0, WholePMT11_H/2-PMTFrontBase11_H-PMTRearBase11_H-PMTTail11_H-PMTPCB11_T/2), lv_PMTPCB11, "PMTPCB11", lv_PMT11, false, 97);

	//Pillar
		G4VSolid* sol_Pillar11 = new G4Tubs("Pillar11", 0, Pillar11_D/2, Pillar11_H/2, 0, 360.*deg);
		G4LogicalVolume* lv_Pillar11 = new G4LogicalVolume(sol_Pillar11, mat_Al, "Pillar11");

		std::vector < G4double > yPos11_pillar;
		for(G4int i=1; i<7;++i){
			if(i<4){
				yPos11_pillar.push_back(CrossBar_W/2 + PMTbox11 - 0.5 * PMTFrontBase11_W *(2*i-1) );
			}
			else {
				yPos11_pillar.push_back(-(CrossBar_W/2 + 0.5 * PMTFrontBase11_W *(2*i-7)) );
			}
		}
		for(G4int i=0; i<6;++i){
			new G4PVPlacement(0,G4ThreeVector((CrossBar_W/2+PMTbox11+21.37*mm+Pillar11_D/2) ,yPos11_pillar[i], (RearCase11_H-RearCase11_T)/2-Pillar11_H/2 ), lv_Pillar11, "Pillar11", lv_RearCase11, false, 98);
			new G4PVPlacement(0,G4ThreeVector(-(CrossBar_W/2+PMTbox11+21.37*mm+Pillar11_D/2),yPos11_pillar[i], (RearCase11_H-RearCase11_T)/2-Pillar11_H/2 ), lv_Pillar11, "Pillar11", lv_RearCase11, false, 98);
		}

	// Holder
		G4VSolid* sol_Holder11 = new G4Box("Holder11", Holder11_W/2, Holder11_H/2, Holder11_T/2);
		G4LogicalVolume* lv_Holder11 = new G4LogicalVolume(sol_Holder11, mat_Al, "Holder11");
		for(G4int i=1; i<7; i++){
				if(i<4){
					new G4PVPlacement(0, G4ThreeVector( 0, (CrossBar_W+PMTFrontBase11_W*(2*i-1))/2, (RearCase11_H-RearCase11_T)/2-WholePMT11_H-Holder11_T/2), lv_Holder11, "Holder11", lv_RearCase11, false, 98);
				}
				else {
					new G4PVPlacement(0, G4ThreeVector(0,  -(CrossBar_W+PMTFrontBase11_W*(2*i-7) )/2, (RearCase11_H-RearCase11_T)/2-WholePMT11_H-Holder11_T/2), lv_Holder11, "Holder11", lv_RearCase11, false, 98);
				}
			}

	new G4PVPlacement(0, G4ThreeVector(0,0,-WholeDet11_H/2+rear_Al_T/2), lv_rear_Al, "rear_Al", lv_Det11, false, 700);
	new G4LogicalSkinSurface("alSurface11",lv_Housing11, groundWhitepainted);
	new G4LogicalSkinSurface("pmtSurface11",lv_PMTPhotoCathode11, CathodeSurface);

	new G4LogicalSkinSurface("paintsideSurface11",lv_Paint11_1, groundWhitepainted);
	new G4LogicalSkinSurface("alSurface11",	lv_Housing11, groundWhitepainted);
	new G4LogicalSkinSurface("pmtSurface11",lv_PMTPhotoCathode11, CathodeSurface);


	//MASK
	G4VSolid *sol_Mask = new G4Box("Mask", .5 * Mask_W, .5 * Mask_W, .5 * Mask_T);
	G4LogicalVolume *lv_Mask = new G4LogicalVolume(sol_Mask, mat_Air, "Mask");
	//new G4PVPlacement(Detector_RotMat, G4ThreeVector(0., 0., M2D_Dist + .5 * Mask_T), lv_Mask, "Mask", lv_World, false, 0);

	new G4PVPlacement(Detector_RotMat, G4ThreeVector((M2D_Dist +.5*Mask_T)*mm, 0.,-200 *mm), lv_Mask, "Mask", lv_World, false, 0);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector((M2D_Dist +.5*Mask_T)*mm, 0.,200*mm), lv_Mask, "Mask", lv_World, false, 0);
	

	G4VSolid *sol_MaskPix = new G4Box("MaskPix", .5 * Mask_W / (G4double)Mask_nPix, .5 * Mask_W / (G4double)Mask_nPix, .5 * Mask_T);
	G4LogicalVolume *lv_MaskPix = new G4LogicalVolume(sol_MaskPix, mat_W, "MaskPix");
	MaskParam *para_MaskPix = new MaskParam(Mask_W, Mask_W, Mask_nPix, Mask_nPix, MaskPatternFile);
	new G4PVParameterised("MaskPix", lv_MaskPix, lv_Mask, kXAxis, Mask_nPix * Mask_nPix, para_MaskPix);

	//G4VSolid* sol_Frame_main = new G4Tubs ("Frame_m",0, 29*cm,0.3*cm,0,360.*deg);
	G4VSolid *sol_Frame_main = new G4Box("Frame_m", 0.54 * Mask_W, 0.54 * Mask_W, 0.5*Mask_T);
	G4VSolid *sol_Frame_hole = new G4Box("Frame_h", 0.54 * Mask_W, 0.54 * Mask_W, 0.5*Mask_T);
	G4VSolid *sol_Frame = new G4SubtractionSolid("Frame", sol_Frame_main, sol_Frame_hole);
	G4LogicalVolume *lv_Frame = new G4LogicalVolume(sol_Frame, mat_Stainless, "Frame");
	//new G4PVPlacement(0, G4ThreeVector(0., 0., M2D_Dist + .5 * Mask_T), lv_Frame, "Frame", lv_World, false, 0);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector((M2D_Dist +.5*Mask_T)*mm, 0., -200*mm), lv_Frame, "Frame", lv_World, false, 0);
	new G4PVPlacement(Detector_RotMat, G4ThreeVector((M2D_Dist +.5*Mask_T)*mm, 0., 200*mm), lv_Frame, "Frame", lv_World, false, 0);
	

	
/*
	//Mask frame box
	G4VSolid *sol_Frame_box = new G4Box("Frame_b", 0.6 * Mask_W + 0.5 *Mask_T, 0.6 * Mask_W + 0.5 *Mask_T, (Mask_T+ M2D_Dist + WholeDet1_H + WholeDet11_H + 50*mm)*0.5);
	G4VSolid *sol_Frame_air = new G4Box("Frame_a", 0.6 * Mask_W, 0.6 * Mask_W, (Mask_T+ M2D_Dist + WholeDet1_H + WholeDet11_H + 60*mm)*0.5);
	G4VSolid *sol_Frame_side = new G4SubtractionSolid("Frame_side", sol_Frame_box, sol_Frame_air);
	G4LogicalVolume *lv_Frame_side = new G4LogicalVolume(sol_Frame_side, mat_Stainless, "Frame_side");
	new G4PVPlacement(0, G4ThreeVector(0., 0., M2D_Dist + .5 * Mask_T-(M2D_Dist + WholeDet1_H + WholeDet11_H + 50*mm)/2), lv_Frame_side, "Frame_side", lv_World, false, 0);
*/



	// -- Visualization -- //
	// Colors
	VisColor();

	// Visualization
	G4VisAttributes* va_World = new G4VisAttributes(*white);
	va_World->SetForceWireframe(true);
	lv_World->SetVisAttributes(va_World);

	G4VisAttributes* va_Inv = new G4VisAttributes(*white);
	va_Inv->SetVisibility(false);
	lv_Det1->SetVisAttributes(va_Inv);
	lv_PMT1->SetVisAttributes(va_Inv);
	lv_PMTFrontBase1->SetVisAttributes(va_Inv);
	lv_Det11->SetVisAttributes(va_Inv);
	lv_PMT11->SetVisAttributes(va_Inv);
	lv_PMTFrontBase11->SetVisAttributes(va_Inv);

	G4VisAttributes* va_Al = new G4VisAttributes(*gray);
	lv_Pillar1->SetVisAttributes(va_Al);
	lv_Holder1->SetVisAttributes(va_Al);
	lv_Pillar11->SetVisAttributes(va_Al);
	lv_Holder11->SetVisAttributes(va_Al);

	G4VisAttributes* va_Plastic = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2));
	va_Plastic->SetForceSolid(true);
	lv_PMTTail1->SetVisAttributes(va_Plastic);
	lv_PMTTail2->SetVisAttributes(va_Plastic);
	lv_PMTTail11->SetVisAttributes(va_Plastic);
	lv_PMTTail11_2->SetVisAttributes(va_Plastic);

	G4VisAttributes* va_PCB = new G4VisAttributes(G4Colour(0.15, 0.4, 0.15));
	va_PCB->SetForceSolid(true);
	lv_PMTPCB1->SetVisAttributes(va_PCB);
	lv_PMTPCB11->SetVisAttributes(va_PCB);

	G4VisAttributes* va_Glass = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.4));
	va_Glass->SetForceSolid(true);
	lv_PMTFrontBase1_Shell->SetVisAttributes(va_Glass);
	lv_PMTRearBase1_Shell->SetVisAttributes(va_Glass);
	lv_PMTFrontBase11_Shell->SetVisAttributes(va_Glass);
	lv_PMTRearBase11_Shell->SetVisAttributes(va_Glass);
	//lv_innerPhantom-> SetVisAttributes(va_Glass);

	G4VisAttributes* va_Cathode = new G4VisAttributes(*darkOrange3);
	lv_PMTPhotoCathode1->SetVisAttributes(va_Cathode);
	lv_PMTPhotoCathode11->SetVisAttributes(va_Cathode);

	G4VisAttributes* va_Crystal = new G4VisAttributes(*white);
	lv_Crystal1->SetVisAttributes(va_Crystal);
	lv_Crystal11->SetVisAttributes(va_Crystal);
	lv_Phantom->SetVisAttributes(va_Crystal);

	G4VisAttributes* va_OptWin = new G4VisAttributes(*skyBlue);
	lv_OptWin1->SetVisAttributes(va_OptWin);
	lv_OptWin11->SetVisAttributes(va_OptWin);

	G4VisAttributes* va_CrossBar = new G4VisAttributes(*yellow);
	lv_CrossBar1->SetVisAttributes(va_CrossBar);
	lv_CrossBar1_H1->SetVisAttributes(va_CrossBar);
	lv_CrossBar1_H2->SetVisAttributes(va_CrossBar);
	lv_CrossBar11->SetVisAttributes(va_CrossBar);
	lv_CrossBar11_H1->SetVisAttributes(va_CrossBar);
	lv_CrossBar11_H2->SetVisAttributes(va_CrossBar);
	return pv_World;
}

void DetectorConstruction::ConstructSDandField(){
	G4VSensitiveDetector* CrystalSD = new DEPosSD("DEPos");
	G4SDManager::GetSDMpointer()->AddNewDetector(CrystalSD);

	G4VSensitiveDetector* PMTSD = new PhotCntSD("PhotCnt");
	G4SDManager::GetSDMpointer()->AddNewDetector(PMTSD);

	SetSensitiveDetector("Crystal1", CrystalSD);
	SetSensitiveDetector("PMTPhotoCathode1", PMTSD);

	if(fmode==2){
		SetSensitiveDetector("Crystal11", CrystalSD);
		SetSensitiveDetector("PMTPhotoCathode11", PMTSD);
	}
}

void DetectorConstruction::DefineDimensions(){
	World_Size = 25.*m;

	// User geometry
	// For source
	SrcCase_D = 25.4*mm, SrcCase_T = 3.2*mm;
	Src2DetDist = 100.*cm;

	// For collimator
	Collimator_W = 40.*mm;
	CollimatorBox_T = 20.*mm, CollimatorFace_T = 0.*mm;
	CollimatorHole_D = 1.*mm;

	// For CC
	//Sc2Ab_Dist = 250.*mm;

	// For detector type 1 (36 ch.)
	CrossBar_W = 8.*mm, CrossBar_T = 8.*mm, CrossBar_H = 319.55*mm;
	Crystal_W = 146.*mm, Crystal_Wf = 192.3*mm, Crystal_T = 20.*mm;
	OptGlue_T = 0.5*mm, OptWin_W = 158.*mm, OptWin1_T = 10.*mm;
	Paint_T = 2.0*mm, OptGrease_T = 1.*mm;
	Housing1_W = 434.*mm, Housing1_V = 376.*mm, Housing_T = 1.5*mm+Paint_T+Crystal_T+OptGlue_T+OptWin1_T+OptGrease_T;
	PMTFrontGlass_T = 1.5*mm;
	PMTPhotoCathode_W = 48.*mm, PMTPhotoCathode_T = 0.5*mm;
	PMTFrontBase_W = 51.*mm, PMTFrontBase_T = 40.*mm;
	PMTbox= 3 * PMTFrontBase_W;
	PMTRearBase_D = 51.*mm, PMTRearBase_H = 47.*mm;
	PMTTail_D = 18.7*mm, PMTTail_H =18.4*mm,PMTTail2_H = 15.6*mm;
	PMTPCB_D = PMTRearBase_D, PMTPCB_T = 1.5*mm;
	WholePMT_W = PMTFrontBase_W, WholePMT_T = 123.5*mm;
	Pillar1_D = 8.*mm, Pillar1_H = 123.5*mm; 	//Pillar1_D = 9.5*mm, Pillar1_H = 146.5*mm;
	RearCase_W = 394.6*mm, RearCase1_D = 337.2*mm, RearCase1_H = 164*mm, RearCase_T = 1.5*mm;
	rear_Al_T=6*mm, rear_Al_W=430.*mm,rear_Al_H=372.*mm;
	WholeDet1_W = Housing1_W, WholeDet1_D = Housing1_V, WholeDet1_H = Housing_T+RearCase1_H+rear_Al_T;
	Holder1_W = 384.61, Holder1_H = 11.09*mm, Holder1_T = 15*mm;

	// For detector type 1' (36 ch.)
	Crystal11_W = 146.*mm, Crystal11_Wf = 192.3*mm, Crystal11_T =30.*mm;
	OptGlue11_T = 0.5*mm, OptWin11_W = 158.*mm, OptWin11_T = 10.*mm;
	Paint11_T = 2.0*mm, OptGrease11_T = 1.*mm;
	Housing11_W = 434.*mm, Housing11_V = 376.*mm,Housing11_T = 1.5*mm+Paint11_T+Crystal11_T+OptGlue11_T+OptWin11_T+OptGrease11_T;
	PMTFrontGlass11_T = 1.5*mm;
	PMTPhotoCathode11_W = 48.*mm, PMTPhotoCathode11_T = 0.5*mm;
	PMTFrontBase11_W = 51.*mm, PMTFrontBase11_H = 48.*mm;
	PMTbox11= 3 * PMTFrontBase11_W;
	PMTRearBase11_D = 51.*mm, PMTRearBase11_H = 56.*mm;
	PMTTail11_D = 18.7*mm, PMTTail11_H = 18.5*mm, PMTTail12_H = 20.8*mm;
	PMTPCB11_D = PMTRearBase11_D, PMTPCB11_T = 1.5*mm;
	WholePMT11_W = PMTFrontBase11_W, WholePMT11_H = 143.5*mm;
	Pillar11_D = 8.*mm, Pillar11_H = 143.5*mm;
	RearCase11_W = 394.6*mm, RearCase11_D=  337.2*mm,  RearCase11_H = 174.*mm, RearCase11_T = 1.5*mm;
	Holder11_W = 384.61, Holder11_H = 11.09*mm, Holder11_T = 15*mm;
	WholeDet11_W = Housing11_W, WholeDet11_D = Housing11_V, WholeDet11_H = Housing11_T+RearCase11_H+rear_Al_T;

	// For CC
	Sc2Ab_Dist= 250*mm-Crystal_T/2-OptGlue_T-OptWin1_T-OptGrease_T-RearCase1_H-rear_Al_T-1.5*mm-Paint_T-Crystal11_T/2;

	phantom_S = 58.6;
	phantom_H = 88.4;

	Phantom_outerD1 = phantom_S *cm;
	Phantom_outerD2 = Phantom_outerD1/8;
	phantomholeD = Phantom_outerD1*3/16;
	Phantom_H = phantom_H*cm;
}
void DetectorConstruction::DefineMaterials(){
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	// NistElement definition
	G4Element* nistElH = nist->FindOrBuildElement("H");
	G4Element* nistElC = nist->FindOrBuildElement("C");
	G4Element* nistElO = nist->FindOrBuildElement("O");
	G4Element* nistElAl = nist->FindOrBuildElement("Al");
	G4Element* nistElSi = nist->FindOrBuildElement("Si");
	G4Element* nistElTl = nist->FindOrBuildElement("Tl");
	G4Element* nistElCs = nist->FindOrBuildElement("Cs");
	G4Element* nistElSb = nist->FindOrBuildElement("Sb");
	G4Element* nistElRb = nist->FindOrBuildElement("Rb");
	//G4Element* nistElPb = nist->FindOrBuildElement("Pb");
//	G4Element* nistElFe = nist->FindOrBuildElement("Fe");

	// NistMaterial definition
	nist->FindOrBuildMaterial("G4_Galactic");
	nist->FindOrBuildMaterial("G4_AIR");
	nist->FindOrBuildMaterial("G4_Al");
	nist->FindOrBuildMaterial("G4_Pb");
	nist->FindOrBuildMaterial("G4_Cu");
	nist->FindOrBuildMaterial("G4_W");
	nist->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE"); // paint
	nist->FindOrBuildMaterial("G4_POLYETHYLENE"); // source case
	nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
	nist->FindOrBuildMaterial("G4_Fe");
	nist->FindOrBuildMaterial("G4_STAINLESS-STEEL"); // source case
	G4Material* BSiO2 = nist->FindOrBuildMaterial("G4_Pyrex_Glass"); // PMT front glass

	// UserMaterial definition
	G4Material* NaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
	G4Material* NaITl = new G4Material("NaITl", 3.67*g/cm3, 2, kStateSolid); // crystal
	NaITl->AddMaterial(NaI, 99.6*perCent);
	NaITl->AddElement(nistElTl, 0.4*perCent);

	G4Material* Polydimethylsiloxane = new G4Material("Polydimethylsiloxane", 0.97*g/cm3, 4, kStateSolid); // BC630 material
	Polydimethylsiloxane->AddElement(nistElSi, 1);
	Polydimethylsiloxane->AddElement(nistElO, 1);
	Polydimethylsiloxane->AddElement(nistElC, 2);
	Polydimethylsiloxane->AddElement(nistElH, 6);

	G4Material* FusedSilica = new G4Material("FusedSilica", 2.201*g/cm3, 2, kStateSolid); // Optical window
	FusedSilica->AddElement(nistElSi, 1);
	FusedSilica->AddElement(nistElO, 2);

	G4Material* BC630 = new G4Material("BC630", 1.06*g/cm3, 2, kStateLiquid); // Optical grease
	BC630->AddMaterial(Polydimethylsiloxane, 0.95);
	BC630->AddMaterial(FusedSilica, 0.05);

	G4Material* BialkaliCathode = new G4Material("BialkaliCathode", 3.*g/cm3, 3, kStateSolid);
	BialkaliCathode->AddElement(nistElSb, 1);
	BialkaliCathode->AddElement(nistElRb, 1);
	BialkaliCathode->AddElement(nistElCs, 1);

	G4Material* Al_2 = nist->FindOrBuildMaterial("G4_Al");
	G4Material* SiO2 = new G4Material("SiO2", 0.9*g/cm3, 2, kStateSolid); // BC630 material
	SiO2->AddElement(nistElSi, 1);
	SiO2->AddElement(nistElO, 2);
	G4Material* Al2O3 = new G4Material("Al2O3", 0.9*g/cm3, 2, kStateSolid); // BC630 material
	Al2O3->AddElement(nistElAl, 2);
	Al2O3->AddElement(nistElO, 3);
	G4Material* Concrete = new G4Material("Concrete", 2.35*g/cm3, 2, kStateSolid);
	Concrete->AddMaterial(SiO2, 1);
	Concrete->AddMaterial(Al2O3, 1);

/*	G4double density_Al = Al_2->GetDensity();
	G4double density_Concrete = Concrete->GetDensity();*/

	G4Material* DAW = new G4Material("DAW", 0.3*g/cm3, 2, kStateSolid);
	DAW->AddElement(nistElC,3);
	DAW->AddElement(nistElH,6);
	
	//    Add properties   //

	// NaITl properties
	const G4int NaITl_NumEntries = 32;

	G4double NaITl_Energies[NaITl_NumEntries] = {1.2*eV, 2.265*eV, 2.303*eV, 2.341*eV, 2.381*eV, 2.422*eV, 2.464*eV
			, 2.508*eV, 2.554*eV, 2.601*eV, 2.650*eV, 2.701*eV, 2.754*eV, 2.810*eV, 2.867*eV, 2.927*eV,
			2.989*eV , 3.054*eV , 3.122*eV , 3.193*eV , 3.267*eV , 3.345*eV , 3.427*eV ,3.513*eV, 3.603*eV
			,3.698*eV, 3.798*eV, 3.903*eV, 4.015*eV ,4.133*eV, 4.259*eV, 4.95*eV};
	G4double NaITl_FastComponents[NaITl_NumEntries] = { 0., 0 ,0.0114, 0.0228, 0.0415 ,0.0622 ,0.0910 ,0.1272 ,0.1685 ,	0.2458 ,0.3380 ,0.4647 ,
			0.5915 ,0.7160	,0.8284	,0.9407 ,0.9870 ,0.9914 ,0.9654 ,0.8473 ,0.6527	,0.5192 ,0.4666, 0.4304, 0.3872 ,0.3080 ,0.2121 ,0.1162	,
			0.04095 ,0.0117 ,0.009 ,0 };
	G4double NaITl_RIndices[NaITl_NumEntries] = {1.74 ,1.78 ,1.78 ,1.78 ,1.78 ,1.79 ,1.79 ,
			1.79 , 1.80 ,1.80 ,1.80 ,1.80, 1.80 ,1.81 , 1.81 ,1.82 ,
			1.82 , 1.83, 1.83, 1.83, 1.84 , 1.85 , 1.86 ,1.87 , 1.88,
			1.88 ,1.88 , 1.90 , 1.92 , 1.93 , 1.95 , 2.08	};
	G4double NaITl_AbsorptionLength[NaITl_NumEntries] = { 10000.0 * mm, 10000.0* mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm,
			10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0* mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0* mm, 10000.0 * mm, 10000.0 * mm,
			10000.0 * mm, 10000.0* mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0* mm, 10000.0 * mm, 10000.0 * mm, 10000.0 * mm, 10000.0* mm, 10000.0 * mm,
			10000.0 * mm, 10000.0 * mm, 10000.0* mm, 10000.0 * mm };

	G4MaterialPropertiesTable* NaITl_ProP = new G4MaterialPropertiesTable();
	NaITl_ProP->AddProperty("FASTCOMPONENT", NaITl_Energies, NaITl_FastComponents, NaITl_NumEntries);
	NaITl_ProP->AddProperty("RINDEX", NaITl_Energies, NaITl_RIndices, NaITl_NumEntries);
	NaITl_ProP->AddProperty("ABSLENGTH", NaITl_Energies, NaITl_AbsorptionLength, NaITl_NumEntries);
	NaITl_ProP->AddConstProperty("SCINTILLATIONYIELD", 38000./MeV);
	NaITl_ProP->AddConstProperty("RESOLUTIONSCALE", 1.0);
	NaITl_ProP->AddConstProperty("FASTTIMECONSTANT", 230.*ns);
	NaITl->SetMaterialPropertiesTable(NaITl_ProP);

	// Fused silica (optical window)
	G4MaterialPropertiesTable* FusedSilica_Prop = new G4MaterialPropertiesTable();
	G4double FusedSilica_Energies[5] = {1.*eV, 2.*eV, 3.*eV, 4.*eV, 5.*eV};
	G4double FusedSilica_RIndices[5] = {1.56, 1.56, 1.56, 1.56, 1.56};
	FusedSilica_Prop->AddProperty("RINDEX", FusedSilica_Energies, FusedSilica_RIndices, 5);
	FusedSilica->SetMaterialPropertiesTable(FusedSilica_Prop);

	// BSiO2 (borosilicated glass)
	G4MaterialPropertiesTable* BSiO2_Prop = new G4MaterialPropertiesTable();
	G4double BSiO2_Energies[4] = {1.*eV, 2.*eV, 3.*eV, 4.*eV};
	G4double BSiO2_RIndices[4] = {1.48, 1.48, 1.48, 1.48};
	BSiO2_Prop->AddProperty("RINDEX", BSiO2_Energies, BSiO2_RIndices, 4);
	BSiO2->SetMaterialPropertiesTable(BSiO2_Prop);

	// BC630 (optical grease)
	G4MaterialPropertiesTable* BC630_Prop = new G4MaterialPropertiesTable();
	G4double BC630_Energies[2] = {1.7748*eV, 4.4371*eV};
	G4double BC630_RIndices[2] = {1.465, 1.465};
	BC630_Prop->AddProperty("RINDEX", BC630_Energies, BC630_RIndices, 2);
	BC630->SetMaterialPropertiesTable(BC630_Prop);
}
void DetectorConstruction::VisColor(){
	white = new G4VisAttributes(G4Colour());
	white -> SetVisibility(true);
	white -> SetForceSolid(true);

	blue = new G4VisAttributes(G4Colour(0., 0., 1.));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);

	gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	gray-> SetVisibility(true);
	gray-> SetForceSolid(true);

	red = new G4VisAttributes(G4Colour(1., 0., 0.));
	red-> SetVisibility(true);
	red-> SetForceSolid(true);

	yellow = new G4VisAttributes(G4Colour(1., 1., 0.));
	yellow-> SetVisibility(true);
	yellow-> SetForceSolid(true);

	green = new G4VisAttributes(G4Colour(25/255., 255/255., 25/255.));
	green -> SetVisibility(true);
	green -> SetForceSolid(true);

	darkGreen = new G4VisAttributes(G4Colour(0/255., 100/255., 0/255.));
	darkGreen -> SetVisibility(true);
	darkGreen -> SetForceSolid(true);

	darkOrange3 = new G4VisAttributes(G4Colour(205/255., 102/255., 0/255.));
	darkOrange3 -> SetVisibility(true);
	darkOrange3 -> SetForceSolid(true);

	skyBlue = new G4VisAttributes(G4Colour(135/255., 206/255., 235/255.));
	skyBlue -> SetVisibility(true);
	skyBlue -> SetForceSolid(true);

}

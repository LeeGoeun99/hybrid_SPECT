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

#ifndef MaskParam_hh_
#define MaskParam_hh_

#include "globals.hh"
#include "G4VNestedParameterisation.hh"
// #include "G4VPVParameterisation.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include <fstream>
#include <vector>

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

class MaskParam: public G4VNestedParameterisation{
public:
    MaskParam(G4double _sizeX, G4double _sizeY, G4int _nX, G4int _nY, G4String fileName);
    ~MaskParam();
    void ComputeTransformation(const G4int CpNo,
        G4VPhysicalVolume* physVol) const;
    G4Material* ComputeMaterial(G4VPhysicalVolume* physVol,
        const G4int CpNo,
        const G4VTouchable* parentTouch=0);
    G4int GetNumberOfMaterials() const;
    G4Material* GetMaterial(G4int CpNo) const;

private:  // Dummy declarations to get rid of warnings ...
    void ComputeDimensions (G4Box&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
        const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
        const G4VPhysicalVolume*) const {}

private:
    G4Material* mat_Air;
    G4Material* mat_W;
    G4VisAttributes* graySol;
    G4VisAttributes* grayWired;
    G4double sizeX, sizeY;
    G4int nX, nY;
    std::vector < std::vector < G4bool > > mask_Pattern;
};

#endif
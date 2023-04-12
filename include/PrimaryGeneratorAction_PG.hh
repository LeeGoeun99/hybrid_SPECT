#ifndef PrimaryGeneratorAction_PG_hh_
#define PrimaryGeneratorAction_PG_hh_

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "DetectorConstruction.hh"
#include "Randomize.hh"
class PrimaryGeneratorAction_PG: public G4VUserPrimaryGeneratorAction{
public:
	PrimaryGeneratorAction_PG(DetectorConstruction* );
	virtual ~PrimaryGeneratorAction_PG();
	virtual void GeneratePrimaries(G4Event*);

private:

	G4ParticleGun* fPrimary;
	std::vector < G4double > vx;
	std::vector < G4double > vy;
	std::vector < G4double > vz;
	//G4double xPos[12];
	//std::vector < G4double > xPos;
	//std::vector < G4double > Energy;
	DetectorConstruction* fDetectorConstruction;

	G4double const PI = 3.141592;

	G4double innersouce_length ;
	G4double ceta ;
	G4double ceta_degree;
	G4double degree_pos;
	G4double xpos;
	G4double zpos;

inline G4ThreeVector G4RandomDirection()
{
  G4double u, v, b;
  do {
    u = 2.*G4UniformRand() - 1.;
    v = 2.*G4UniformRand() - 1.;
    b = u*u + v*v;
  } while (b > 1.);
  G4double a = 2.*std::sqrt(1. - b);
  return G4ThreeVector(a*u, a*v, 2.*b - 1.);
}

inline G4ThreeVector G4RandomDirection(G4double cosTheta)
{
  G4double z   = (1. - cosTheta)*G4UniformRand() + cosTheta;
  G4double rho = std::sqrt((1.+z)*(1.-z));
  G4double phi = CLHEP::twopi*G4UniformRand();
  return G4ThreeVector((-1)*z,rho*std::cos(phi), rho*std::sin(phi));
}
};

#endif



#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"
#include "G4GeneralParticleSource.hh" 
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
using namespace std;
using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	
	// Default values  
	
	particleGun = new G4GeneralParticleSource();
	particleGun->SetCurrentSourceIntensity (1);
	particleGun->SetParticlePosition(G4ThreeVector());
       
	
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//create vertex
	
//        G4int NumberofParticles=RandGauss::shoot(10, 3); 
        particleGun->SetNumberOfParticles(1);
        //particleGun->SetNumberOfParticles(2);
//        G4double ParticleTime=RandGauss::shoot(5, 1.274); 
        
        //cout << NumberofParticles<<"::"<<ParticleTime<<endl;
        particleGun->SetParticleTime(1);

        //G4double ParticleEnergy=RandExponential::shoot(5);
        //cout<< particleGun->GetParticleEnergy()<<endl;
	particleGun->GeneratePrimaryVertex(anEvent);
}

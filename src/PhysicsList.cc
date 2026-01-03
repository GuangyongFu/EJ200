/////////////////////////////////////////////////////////
//
// Apr/2015  E. Nacher -> PhysicsList.cc
//
// Physics List is a simplification of the LXePhysicsList 
// ($G4INSTALL/examples/extended/optical/LXe). EM physics 
// just registering G4EmStandardPhysics and no Decay Physics.
//
////////////////////////////////////////////////////////

#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"

#include "G4NeutronHPElastic.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "QGSP_BIC_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4ParticleTypes.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayPhysics.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
	// default cut value  (0.1 mm)
	defaultCutValue = 0.001*mm;
	G4int ver=1;
	// EM Physics
	//RegisterPhysics(new G4EmStandardPhysics(verb));
//	RegisterPhysics(new G4EmPenelopePhysics(verb));
//	RegisterPhysics(new QGSP_BIC_HP(verb));

	//neutron Physics
//	G4NeutronHPElastic* neutronHPElastic = new G4NeutronHPElastic();
//	RegisterPhysics(neutronHPElastic);
//	RegisterPhysics( new G4HadronElasticPhysicsHP(verb) );
//	RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));

	// Ion Elastic scattering
//	RegisterPhysics( new G4IonElasticPhysics(verb));
      
	// Ion Inelastic physics
//	RegisterPhysics( new G4IonPhysics(verb));

  // EM Physics
	this->RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
	this->RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
	this->RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
	this->RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );

   // Hadron Physics
	this->RegisterPhysics(  new G4HadronPhysicsQGSP_BIC_HP(ver));

  // Stopping Physics
	this->RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
	this->RegisterPhysics( new G4IonPhysics(ver));
  

	// Optical Physics
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
	RegisterPhysics(opticalPhysics);
	
	//opticalPhysics->SetScintillationYieldFactor(1.0);
	//opticalPhysics->SetScintillationExcitationRatio(0.);
	
	//opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
	
}


PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts(){
	
	SetCutsWithDefault();
//  if (verboseLevel >0)
  //{
    //G4cout << "PhysicsList::SetCuts:";
  //  G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
//  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
//  SetCutValue(cutForGamma, "gamma");
 // SetCutValue(cutForElectron, "e-");
//  SetCutValue(cutForPositron, "e+");
//  SetCutValue(cutForProton, "proton");
//  SetCutValue(cutForneutron, "neutron");
//  if (verboseLevel>0) { DumpCutValuesTable(); }
}

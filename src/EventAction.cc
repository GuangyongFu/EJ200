#include "EventAction.hh"
#include "Analysis.hh"
//#include "Randomize.hh" // do we really need this?
#include <iomanip>

#include "RunAction.hh"
#include "G4GeneralParticleSource.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh"

#include "globals.hh"
#include "CLHEP/Random/RandGauss.h"

//#include "iomanip"

using namespace std;
using namespace CLHEP;

EventAction::EventAction(RunAction*)
:G4UserEventAction()//, EdepInCrystal(0.), nAbsPhotons(0.), absTime(0.)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
	// initialisation per event
	EdepInCrystal = 0.;
	nAbsPhotons = 0.;
	absTime = 0.;
        totEnergyDepCathod = 0.;
        t_min = 100;
        t_max = 0; 

    // Get energy of the primary particles (created with GPS)
   const PrimaryGeneratorAction* generatorAction
   = static_cast<const PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
    k_primary = particleGun->GetParticleEnergy();
    // G4cout << "\n And the energy is: BeginOfEventAction" << k_primary << "MeV" << G4endl;

}

void EventAction::EndOfEventAction(const G4Event* evt)
{

	  // Accumulate statistics
	  //

	  // get analysis manager
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	  G4double x= RandGauss::shoot(); 
 	  G4double FWHM = 0.001 + 0.05086 * sqrt(EdepInCrystal+0.030486*EdepInCrystal);
  	  G4double sigma = FWHM/(2*sqrt(2*log(2)));
          G4double EdepFlu = EdepInCrystal + sigma * x;
	  G4double MeVee=nAbsPhotons/10000.;
	  // fill histograms
	  analysisManager->FillH1(0, EdepFlu);
	  analysisManager->FillH1(1, MeVee);
	  analysisManager->FillH1(2, t_min);
          analysisManager->FillH1(3, t_max);
	  analysisManager->FillH1(4, k_primary);
          analysisManager->FillH1(5, totEnergyDepCathod);

	  // fill ntuple
	  analysisManager->FillNtupleDColumn(0, EdepInCrystal);
	  analysisManager->FillNtupleDColumn(1, MeVee);
	  analysisManager->FillNtupleDColumn(2, t_min);
          analysisManager->FillNtupleDColumn(3, t_max);
	  analysisManager->FillNtupleDColumn(4, k_primary);
          analysisManager->FillNtupleDColumn(5, totEnergyDepCathod);
	  analysisManager->AddNtupleRow();

  // Print per event (modulo n)
  //
  G4int eventID = 1 + evt->GetEventID();
  //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  G4int printModulo = 100;
  //cout <<"photons" <<MeVee << endl;
  //cout << "Cathod :"<< totEnergyDepCathod<<"::  " <<nAbsPhotons<< endl;
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
  {
    G4cout << "---> End of event: " << eventID << G4endl;
  }
}

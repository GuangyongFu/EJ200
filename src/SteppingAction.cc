#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"

#include "G4RunManager.hh"

#include "fstream"
#include "iomanip"

using namespace std;	 

SteppingAction::SteppingAction(EventAction* EvAct)
:eventAction(EvAct)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	const G4String currentPhysicalName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	
	const G4String particleName
	= aStep->GetTrack()->GetDefinition()->GetParticleName();
	
	if (currentPhysicalName == "Crystal"){
		
		G4double EdepStep = aStep->GetTotalEnergyDeposit();
		
		if (EdepStep > 0.) eventAction->EdepInCrystal = eventAction->EdepInCrystal + EdepStep;

		//count scintillating photons and kill the photons after the first step
		if (particleName == "opticalphoton"){
		 	eventAction->nAbsPhotons++;
		 	eventAction->absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
		 	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		 }


	}
	
	// check if the photon is absorbed in the sensitive volume
	if (currentPhysicalName == "CathodeL"){
		const G4String ProcessName = 
		aStep -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();
		
                
                if (ProcessName == "OpAbsorption"){ 
			
                        eventAction->totEnergyDepCathod += aStep->GetTotalEnergyDeposit();

			//eventAction->nAbsPhotons++;
			
			eventAction->absTime = aStep -> GetPreStepPoint() -> GetGlobalTime();
                        //cout << eventAction->absTime<<endl;
                        if (eventAction->t_max < aStep -> GetPreStepPoint() -> GetGlobalTime()){
                        	eventAction->t_max = aStep -> GetPreStepPoint() -> GetGlobalTime();
                        }
                    
                        if (eventAction->t_min > aStep -> GetPreStepPoint() -> GetGlobalTime()){
                        	eventAction->t_min = aStep -> GetPreStepPoint() -> GetGlobalTime();
                        }
                       // cout << eventAction->t_min<<"::" <<eventAction->t_max<< endl;
		} 
	}
}


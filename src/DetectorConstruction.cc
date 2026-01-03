#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include <complex>
using namespace CLHEP;

DetectorConstruction::DetectorConstruction()
{ }

DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	
	//----------------------------------------------------
	// Material definitions
	//----------------------------------------------------
	
	G4String name, symbol;             //a=mass of a mole;
	G4double a, z, density;            //z=mean number of protons;  
	
	G4int ncomponents, natoms;         
	
	G4double pressure    = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	density     = 1.e-25*g/cm3;
	
	G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, 
										density,kStateGas,temperature,pressure);
	
	
	//
	// define simple Elements
	//
	
	// H
	a = 1.0079*g/mole;
	G4Element* H  = new G4Element(name="Hadrogen" ,symbol="H" , z= 1., a);
	
	// C
	a = 12.011*g/mole;
	G4Element* C  = new G4Element(name="Carbon" ,symbol="C" , z= 6., a);
	
	// O
	a = 15.999*g/mole;
	G4Element* O  = new G4Element(name="Oxygen" ,symbol="O" , z= 8., a);
		
	// Mg
	a = 24.305*g/mole;
	G4Element* Mg  = new G4Element(name="Magnesium"  ,symbol="Mg" , z= 12., a);
	
	// Si
	a = 28.0855*g/mole;
	G4Element* Si  = new G4Element(name="Silicon"  ,symbol="Si" , z= 14., a);
	
	// Cl
	a = 35.453*g/mole;
	G4Element* Cl  = new G4Element(name="Chlorine"  ,symbol="Cl" , z= 17., a);
	
	// K
	a = 39.0983*g/mole;
	G4Element* K  = new G4Element(name="Potassium"  ,symbol="K" , z= 19., a);
	
	// Br
	a = 79.904*g/mole;
	G4Element* Br  = new G4Element(name="Bromine"  ,symbol="Br" , z= 35., a);
	
	// Sb
	a = 121.76*g/mole;
	G4Element* Sb = new G4Element(name="Antimony",symbol="Sb" , z= 51., a);
	
	// Cs
	a = 132.905*g/mole;
	G4Element* Cs = new G4Element(name="Cesium",symbol="Cs" , z= 55., a);
	
	// La
	a = 138.905*g/mole;
	G4Element* La = new G4Element(name="Lanthanum",symbol="La" , z= 57., a);
	
	// Na
	a = 22.989769*g/mole;
	G4Element* Na = new G4Element(name="Sodium",symbol="Na" , z= 11., a);

	// I
	a = 126.90447*g/mole;
	G4Element* I = new G4Element(name="iodine",symbol="I" , z= 53., a);
	
	//
	// define simple materials
	//
	
	// Al reflector
	density = 2.700*g/cm3;
	a = 26.98*g/mole;
	G4Material* AluR = new G4Material(name="AluRef", z=13., a, density);

	
	// EJ200 
//	density = 1.023*g/cm3;
//	G4Material* EJ200 = new G4Material(name="EJ200", density, ncomponents=2);
//	EJ200->AddElement(C, natoms=4.69);
//	EJ200->AddElement(H, natoms=5.17);
	
	density = 1.023*g/cm3;
	G4Material* EJ200 = new G4Material(name="EJ200", density, ncomponents=1);
	EJ200->AddElement(H, natoms=4.69);
        
	// MgO reflector
	density = 2.0*g/cm3;
	G4Material* MgO = new G4Material(name="MgO", density, ncomponents=2);
	MgO->AddElement(Mg, natoms=1);
	MgO->AddElement(O, natoms=1);
		
	// LaCl3
	density = 3.85*g/cm3;
	G4Material* LaCl3 = new G4Material(name="LaCl3", density, ncomponents=2);
	LaCl3->AddElement(Cl, natoms=3);
	LaCl3->AddElement(La, natoms=1);
	
	// LaBr3
	density = 5.08*g/cm3;
	G4Material* LaBr3 = new G4Material(name="LaBr3", density, ncomponents=2);
	LaBr3->AddElement(Br, natoms=3);
	LaBr3->AddElement(La, natoms=1);

	// NaI
	density = 3.67*g/cm3;
	G4Material* NaI = new G4Material(name="NaI", density, ncomponents=2);
	NaI->AddElement(Na, natoms=1);
	NaI->AddElement(I, natoms=1);
	
	// Quartz (fused silica)
	density = 2.20*g/cm3;
	G4Material* Quartz = new G4Material(name="Quartz", density, ncomponents=2);
	Quartz->AddElement(Si, natoms=1);
	Quartz->AddElement(O, natoms=2);
	
	// Photocathode Material 
	// (Bialkali K2CsSb,  Density=?, Thickness=25.*nm?)
	density = 2.00*g/cm3;
	G4Material* K2CsSb = new G4Material(name="K2CsSb", density, ncomponents=3);
	K2CsSb->AddElement(K, natoms=2);
	K2CsSb->AddElement(Cs, natoms=1);
	K2CsSb->AddElement(Sb, natoms=1);
	
	//------------------------------------------------------
	// Optical properties
	//------------------------------------------------------
	
	const G4int nEntries = 2;
	
	G4double PhotonEnergy[nEntries] = {2.38*eV,3.1*eV};
	
	// Al reflector
	
	G4double AlRefractionIndex[nEntries] = {0.89183,0.48774};
	
	G4double AlAbsorptionLength[nEntries] = {6.6E-9*m,6.6E-9*m};
	
	G4MaterialPropertiesTable* AlMPT = new G4MaterialPropertiesTable();
	
//	AlMPT->AddProperty("RINDEX",PhotonEnergy,AlRefractionIndex,
//						nEntries);
	AlMPT->AddProperty("ABSLENGTH",PhotonEnergy,AlAbsorptionLength,
						nEntries);
	
	AluR->SetMaterialPropertiesTable(AlMPT);
	
	// MgO reflector
	
	G4double MgORefractionIndex[nEntries] = {1.7355,1.7355};
	
	G4double MgOAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};
	
	G4MaterialPropertiesTable* MgOMPT = new G4MaterialPropertiesTable();
	
	MgOMPT->AddProperty("RINDEX",PhotonEnergy,MgORefractionIndex,
						nEntries);
	MgOMPT->AddProperty("ABSLENGTH",PhotonEnergy,MgOAbsorptionLength,
						nEntries);
	
	MgO->SetMaterialPropertiesTable(MgOMPT);
	
	
	// NaI
	
//	G4double NaIRefractionIndex[nEntries] = {1.85,1.85};
	
//	G4double NaIAbsorptionLength[nEntries] = {1.*m,1.*m};
	
//	G4MaterialPropertiesTable* NaIMPT = new G4MaterialPropertiesTable();
	
//	NaIMPT->AddProperty("RINDEX",PhotonEnergy,NaIRefractionIndex,
//						  nEntries);
//	NaIMPT->AddProperty("ABSLENGTH",PhotonEnergy,NaIAbsorptionLength,
//						  nEntries);
	
//	G4double ScintEnergy[nEntries] = {3.26*eV,3.44*eV};
//	G4double ScintFast[nEntries] = {1.0,1.0};
//	
//	NaIMPT->AddProperty("FASTCOMPONENT",ScintEnergy,ScintFast,nEntries);
//	
//	NaIMPT->AddConstProperty("SCINTILLATIONYIELD",38./keV); 
//	NaIMPT->AddConstProperty("RESOLUTIONSCALE",1.);
//	NaIMPT->AddConstProperty("FASTTIMECONSTANT",240.*ns);
//	NaIMPT->AddConstProperty("YIELDRATIO",1.);
	
//	NaI->SetMaterialPropertiesTable(NaIMPT);

	// EJ200
	
	EJ200->GetIonisation()->SetBirksConstant(0.1564*mm/MeV);
	G4double EJRefractionIndex[nEntries] = {1.58,1.58};
	
	G4double EJAbsorptionLength[nEntries] = {3.8*m,3.8*m};
	
	G4MaterialPropertiesTable* EJMPT = new G4MaterialPropertiesTable();
	
	EJMPT->AddProperty("RINDEX",PhotonEnergy,EJRefractionIndex,
						  nEntries);
	EJMPT->AddProperty("ABSLENGTH",PhotonEnergy,EJAbsorptionLength,
						  nEntries);
	const G4int num=17;
	G4double ScintEnergy[num] = {3.1*eV,3.05*eV,3.03*eV,3.01*eV,3.00*eV,2.98*eV,2.96*eV,2.92*eV,2.88*eV,2.86*eV,2.81*eV,2.76*eV,2.71*eV,2.64*eV,2.58*eV,2.48*eV,2.38*eV};
	G4double ScintFast[num] = {0.056851,0.123584,0.197078,0.401008,0.598285,0.782153,0.885785,1,0.795144,0.708104,0.691156,0.590623,0.396386,0.282365,0.201784,0.114272,0.0635452};
	G4double ScintSlow[num] = {0.056851,0.123584,0.197078,0.401008,0.598285,0.782153,0.885785,1,0.795144,0.708104,0.691156,0.590623,0.396386,0.282365,0.201784,0.114272,0.0635452};
	EJMPT->AddProperty("SCINTILLATIONCOMPONENT1",ScintEnergy,ScintFast,num);
	EJMPT->AddConstProperty("SCINTILLATIONYIELD",10./keV); 
	EJMPT->AddConstProperty("RESOLUTIONSCALE",2.);
	EJMPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",2.1*ns);
	EJMPT->AddConstProperty("SCINTILLATIONYIELD1",1);
	
	EJ200->SetMaterialPropertiesTable(EJMPT);
		
	

        // Quartz 
	G4double QuartzRefractionIndex[nEntries] = {1.47,1.53};
	
	G4double QuartzAbsorptionLength[nEntries] = {500.0*cm,500.0*cm};
	
	G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
	
	QuartzMPT->AddProperty("RINDEX",PhotonEnergy,QuartzRefractionIndex,
						   nEntries);
	QuartzMPT->AddProperty("ABSLENGTH",PhotonEnergy,QuartzAbsorptionLength,
						   nEntries);
	
	Quartz->SetMaterialPropertiesTable(QuartzMPT);
	
	// K2CsSb (Bialcali Photocathode)
	
	G4double K2CsSbRefractionIndex[nEntries] = {1.47,1.47};
	
	G4double K2CsSbAbsorptionLength[nEntries] = {1.0E-9*m,1.0E-9*m};
	
	G4MaterialPropertiesTable* K2CsSbMPT = new G4MaterialPropertiesTable();
	
	K2CsSbMPT->AddProperty("RINDEX",PhotonEnergy,K2CsSbRefractionIndex,
						   nEntries);
	K2CsSbMPT->AddProperty("ABSLENGTH",PhotonEnergy,K2CsSbAbsorptionLength,
						   nEntries);
	
	K2CsSb->SetMaterialPropertiesTable(K2CsSbMPT);
	
	// Vacuum
	
	G4double vacRefractionIndex[nEntries] = {1.0,1.0};
	
	
	G4MaterialPropertiesTable* vacMPT = new G4MaterialPropertiesTable();
	vacMPT->AddProperty("RINDEX",PhotonEnergy,vacRefractionIndex,
						nEntries);
	
	
	Vacuum->SetMaterialPropertiesTable(vacMPT);
	
	//------------------------------------------------------
	// Detector geometry
	//------------------------------------------------------
	
	//     
	// World
	//
	
	G4double WorldSize= 500.*cm;
	
	G4Box* 
    solidWorld = new G4Box("World",		       	                  //its name
						   WorldSize/2,WorldSize/2,WorldSize/2);  //its size
	
	G4LogicalVolume* 
    logicWorld = new G4LogicalVolume(solidWorld,      	//its solid
									 Vacuum,	        //its material
									 "World");		    //its name
	
	G4VPhysicalVolume* 
    physiWorld = new G4PVPlacement(0,			    //no rotation
								   G4ThreeVector(),	//at (0,0,0)
								   "World",		    //its name
								   logicWorld,		//its logical volume
								   NULL,		    //its mother  volume
								   false,	       	//no boolean operation
								   0);			    //copy number
	
	
	//
	// Detector 
	//
	
	G4double ScintHalfLength =15*cm;
	G4double ScintFrontX = 1.25*cm;
	G4double ScintFrontY = 2.25*cm;
	
    	G4double ReflectorThickness = 0.5*mm;
	G4double ReflectorHalfLength = ScintHalfLength;
	G4double ReflectorFrontX = ScintFrontX+ReflectorThickness;
	G4double ReflectorFrontY = ScintFrontY+ReflectorThickness;
	
	G4double PMTWindowHalfLength = 1.0*mm;
	G4double PMTWindowRadius = 3*cm;
	
	G4double CathodeHalfLength = 0.005*mm;
	G4double CathodeRadius =3.0*cm;
	
	G4double StartPhi = 0.*deg;
	G4double DeltaPhi = 360.*deg;
	
	
	// Reflector
	
	G4Box* solidReflector = new G4Box("Reflector",
										ReflectorFrontX,ReflectorFrontY,ReflectorHalfLength);
	
	G4LogicalVolume* logicReflector = new G4LogicalVolume(solidReflector,AluR,
														  "Reflector");
	
	G4ThreeVector positionReflector = G4ThreeVector(0.*cm,0.*cm,0.*cm);
	
	G4VPhysicalVolume* physiReflector = new G4PVPlacement(0,positionReflector,
														  "Reflector",logicReflector,
														  physiWorld,false,0);
	
	//Crystal
	
	G4Box* solidCrystal = new G4Box("Crystal", 
									  ScintFrontX,ScintFrontY,ScintHalfLength);
	
	G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal,EJ200,
														"Crystal");
	
	G4ThreeVector positionCrystal = G4ThreeVector(0.*cm,0.*cm,
                                                  0.*cm);
	
	G4VPhysicalVolume* physiCrystal = new G4PVPlacement(0,positionCrystal,
														"Crystal",logicCrystal,
														physiReflector,false,0);
	
	
	// PMT window
	
	G4Tubs* solidPMTWindowL = new G4Tubs("PMTWindowL",0.*cm,PMTWindowRadius,
										PMTWindowHalfLength,StartPhi,DeltaPhi);
	
	G4LogicalVolume* logicPMTWindowL = new G4LogicalVolume(solidPMTWindowL,
														  Quartz,"PMTWindowL");
	
	G4ThreeVector positionPMTWindowL = G4ThreeVector(0.*cm,0.*cm,
													-(ReflectorHalfLength+PMTWindowHalfLength));
	
	G4VPhysicalVolume* physiPMTWindowL = new G4PVPlacement(0,positionPMTWindowL,
														  "PMTWindowL",logicPMTWindowL,
														  physiWorld,false,0);
	
	// Photocathode
	
	G4Tubs* solidCathodeL = new G4Tubs("CathodeL",0.*cm,CathodeRadius,
									  CathodeHalfLength,StartPhi,DeltaPhi);
	
	G4LogicalVolume* logicCathodeL = new G4LogicalVolume(solidCathodeL,
														K2CsSb,"CathodeL");
	
	G4ThreeVector positionCathodeL = G4ThreeVector(0.*cm,0.*cm,
												  -(ReflectorHalfLength+2.*PMTWindowHalfLength
												  +CathodeHalfLength));
	
	G4VPhysicalVolume* physiCathodeL = new G4PVPlacement(0,positionCathodeL,
														"CathodeL",logicCathodeL,
														physiWorld,false,0);
	// PMT window
	
	G4Tubs* solidPMTWindowR = new G4Tubs("PMTWindowR",0.*cm,PMTWindowRadius,
										PMTWindowHalfLength,StartPhi,DeltaPhi);
	
	G4LogicalVolume* logicPMTWindowR = new G4LogicalVolume(solidPMTWindowR,
														  Quartz,"PMTWindowR");
	
	G4ThreeVector positionPMTWindowR = G4ThreeVector(0.*cm,0.*cm,
													ReflectorHalfLength+PMTWindowHalfLength);
	
	G4VPhysicalVolume* physiPMTWindowR = new G4PVPlacement(0,positionPMTWindowR,
														  "PMTWindowR",logicPMTWindowR,
														  physiWorld,false,0);
	
	// Photocathode
	
	G4Tubs* solidCathodeR = new G4Tubs("CathodeR",0.*cm,CathodeRadius,
									  CathodeHalfLength,StartPhi,DeltaPhi);
	
	G4LogicalVolume* logicCathodeR = new G4LogicalVolume(solidCathodeR,
														K2CsSb,"CathodeR");
	
	G4ThreeVector positionCathodeR = G4ThreeVector(0.*cm,0.*cm,
												  ReflectorHalfLength+2.*PMTWindowHalfLength
												  +CathodeHalfLength);
	
	G4VPhysicalVolume* physiCathodeR = new G4PVPlacement(0,positionCathodeR,
														"CathodeR",logicCathodeR,
														physiWorld,false,0);
	
	
	//------------------------------------------------------
	// Surfaces and boundary processes
	//------------------------------------------------------
	
	// Reflector - sintillator surface 
	
	G4OpticalSurface* OpRefCrySurface = 
	new G4OpticalSurface("RefCrySurface");
	
	OpRefCrySurface->SetType(dielectric_metal);
	OpRefCrySurface->SetModel(glisur);
	OpRefCrySurface->SetFinish(polished);
	
	G4LogicalBorderSurface* RefCrySurface = 
    new G4LogicalBorderSurface("RefCrySurface",physiCrystal,
							   physiReflector,OpRefCrySurface);
	
	
	// Scintillator - PMT window surface 
	
	G4OpticalSurface* OpCryPMTWinSurface = 
	new G4OpticalSurface("CryPMTWinSurface");
	
	OpCryPMTWinSurface->SetType(dielectric_dielectric);
	OpCryPMTWinSurface->SetModel(glisur);
	OpCryPMTWinSurface->SetFinish(polished);
	
	G4LogicalBorderSurface* CryPMTWinSurfaceL = 
    new G4LogicalBorderSurface("CryPMTWinSurfaceL",physiCrystal,physiPMTWindowL,
							   OpCryPMTWinSurface);
        G4LogicalBorderSurface* CryPMTWinSurfaceR =
    new G4LogicalBorderSurface("CryPMTWinSurfaceR",physiCrystal,physiPMTWindowR,
                                                            OpCryPMTWinSurface);

	
	// PMT window - photocathode surface
	
	G4OpticalSurface* OpPMTWinCathSurface = new G4OpticalSurface("PMTWinCathSurface");
	
	OpPMTWinCathSurface->SetType(dielectric_dielectric);
	OpPMTWinCathSurface->SetModel(glisur);
	OpPMTWinCathSurface->SetFinish(polished);
	
	G4LogicalBorderSurface* PMTWinCathSurfaceL = 
    new G4LogicalBorderSurface("CathodeSurface",physiPMTWindowL,physiCathodeL,
							   OpPMTWinCathSurface);
	G4LogicalBorderSurface* PMTWinCathSurfaceR = 
    new G4LogicalBorderSurface("CathodeSurface",physiPMTWindowR,physiCathodeR,
							   OpPMTWinCathSurface);
	
	
	//------------------------------------------------------
	// visualization attributes
	//------------------------------------------------------
	
//	logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
	
	//Green color for scintillator crystal
	G4VisAttributes* Att1= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
	logicCrystal->SetVisAttributes(Att1);
	
	//Yellow color for reflector
	G4VisAttributes* Att2= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	logicReflector->SetVisAttributes(Att2);
	
	//Blue color for PMT window
	G4VisAttributes* Att3= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	logicPMTWindowL->SetVisAttributes(Att3);
	logicPMTWindowR->SetVisAttributes(Att3);
	
	//White color for the absorber photocathode
	G4VisAttributes* Att4= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	logicCathodeL->SetVisAttributes(Att4);
	logicCathodeR->SetVisAttributes(Att4);
    
	//
	// always return the physical World
	//
	
	return physiWorld;
}

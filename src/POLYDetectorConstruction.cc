//
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
// POLYDetectorConstruction.cc
// \file   MRCP_GEANT4/External/src/POLYDetectorConstruction.cc
// \author Haegin Han
//

#include "POLYDetectorConstruction.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4PSEnergyDeposit.hh"

POLYDetectorConstruction::POLYDetectorConstruction(POLYModelImport* _POLYData)
:worldPhysical(0), container_logic(0), POLYData(_POLYData), POLYLogic(0)
{
	// initialisation of the variables for phantom information
	phantomSize     = POLYData -> GetPhantomSize();
	phantomBoxMin   = POLYData -> GetPhantomBoxMin();
	phantomBoxMax   = POLYData -> GetPhantomBoxMax();
	nOfPOLYrahedrons = POLYData -> GetNumPOLYrahedron();
}

POLYDetectorConstruction::~POLYDetectorConstruction()
{
	delete POLYData;
}

G4VPhysicalVolume* POLYDetectorConstruction::Construct()
{
	SetupWorldGeometry();
	ConstructPhantom();
	PrintPhantomInformation();
	return worldPhysical;
}

void POLYDetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10*10*10 m3)
	//
	G4double worldXYZ = 10. * m;
	G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

	G4VSolid* worldSolid
	  = new G4Box("worldSolid", worldXYZ/2, worldXYZ/2, worldXYZ/2);

	G4LogicalVolume* worldLogical
	  = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");

	worldPhysical
	  = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);

	// Define the phantom container (10-cm margins from the bounding box of phantom)
	//
	G4Box* containerSolid = new G4Box("phantomBox", phantomSize.x()/2 + 10.*cm,
										            phantomSize.y()/2 + 10.*cm,
										            phantomSize.z()/2 + 10.*cm);

	container_logic = new G4LogicalVolume(containerSolid, vacuum, "phantomLogical");

	new G4PVPlacement(0, G4ThreeVector(), container_logic, "PhantomPhysical",
			          worldLogical, false, 0);
	container_logic->SetOptimisation(TRUE);
	container_logic->SetSmartless( 0.5 ); // for optimization (default=2)
}

void POLYDetectorConstruction::ConstructPhantom()
{
	// Define the POLYrahedral mesh phantom as a parameterised geometry
	//
	// solid and logical volume to be used for parameterised geometry
	for(G4int i=0;i<nOfPOLYrahedrons;i++){
		auto mat = POLYData->GetMaterial(POLYData->GetMaterialIndex(i));
		logVec.push_back(new G4LogicalVolume(POLYData->GetTess(i), mat, mat->GetName()+"_log"));
	}

	auto tree = POLYData->GetTessTree();
	for(auto aSet:tree){
		if(aSet.first==-1){
			for(auto sub:aSet.second){
				new G4PVPlacement(0, G4ThreeVector(), logVec[sub],
						logVec[sub]->GetMaterial()->GetName()+"_phy",container_logic, false, POLYData->GetMaterialIndex(sub));
			}
		}
		else{
			for(auto sub:aSet.second){
				new G4PVPlacement(0, G4ThreeVector(), logVec[sub],
					logVec[sub]->GetMaterial()->GetName()+"_phy",logVec[aSet.first], false, POLYData->GetMaterialIndex(sub) );
			}
		}
	}
}

void POLYDetectorConstruction::ConstructSDandField()
{
	// Define detector (Phantom SD) and scorer (eDep)
	//
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	G4String phantomSDname = "PhantomSD";

	// MultiFunctional detector
	G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(phantomSDname);
	pSDman->AddNewDetector( MFDet );

	// scorer for energy depositon in each organ
	MFDet->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));

	// attach the detector to logical volume for parameterised geometry (phantom geometry)
	for(auto logV:logVec){
		SetSensitiveDetector(logV, MFDet);
	}
}

void POLYDetectorConstruction::PrintPhantomInformation()
{
	// print brief information on the imported phantom
	G4cout<< G4endl;
	G4cout.precision(3);
	G4cout<<"   Phantom name               "<<POLYData->GetPhantomName() << " POLY phantom"<<G4endl;
	G4cout<<"   Phantom size               "<<phantomSize.x()<<" * "<<phantomSize.y()<<" * "<<phantomSize.z()<<" mm3"<<G4endl;
	G4cout<<"   Phantom box position (min) "<<phantomBoxMin.x()<<" mm, "<<phantomBoxMin.y()<<" mm, "<<phantomBoxMin.z()<<" mm"<<G4endl;
	G4cout<<"   Phantom box position (max) "<<phantomBoxMax.x()<<" mm, "<<phantomBoxMax.y()<<" mm, "<<phantomBoxMax.z()<<" mm"<<G4endl;
	G4cout<<"   Number of POLYrahedrons     "<<nOfPOLYrahedrons<<G4endl<<G4endl;
}

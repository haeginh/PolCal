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
// POLYActionInitialization.hh
// \file   MRCP_GEANT4/External/include/POLYActionInitialization.hh
// \author Haegin Han
//

#ifndef POLYActionInitialization_h
#define POLYActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "G4String.hh"
#include "POLYRunAction.hh"
#include "POLYPrimaryGeneratorAction.hh"

class POLYModelImport;

// *********************************************************************
// This UserActionInitialization class initializes UserAction classes
// -- BuildForMaster: instantiate UserRunAction class to be used by
//                    G4MTRunManager. This method is not invoked in
//                    the sequential mode.
// -- Build: instantiate UserPrimaryGeneratorAction, UserRunAction, and
//           UserSteppingAction classes.
// *********************************************************************

class POLYActionInitialization : public G4VUserActionInitialization
{
public:
	POLYActionInitialization(POLYModelImport* POLYData,
			                G4String          outputFileName,
							G4Timer*          initTimer);
	virtual ~POLYActionInitialization();

	virtual void BuildForMaster() const;
	virtual void Build() const;

private:
	POLYModelImport* POLYData;
	G4String output;
	G4Timer* initTimer;
};

#endif

    

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
// POLYRunAction.hh
// \file   MRCP_GEANT4/External/include/POLYRunAction.hh
// \author Haegin Han
//

#ifndef POLYRunAction_h
#define POLYRunAction_h 1

#include <ostream>
#include <fstream>
#include <map>

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4UserRunAction.hh"
#include "G4SystemOfUnits.hh"

#include "POLYRun.hh"
#include "POLYPrimaryGeneratorAction.hh"
#include "POLYModelImport.hh"

// *********************************************************************
// The main function of this UserRunAction class is to produce the result
// data and print them.
// -- GenerateRun: Generate POLYRun class which will calculate the sum of
//                  energy deposition.
// -- BeginOfRunAction: Set the RunManager to print the progress at the
//                      interval of 10%.
// -- EndOfRunAction: Print the run result by G4cout and std::ofstream.
//  └-- PrintResult: Method to print the result.
// *********************************************************************

class POLYRunAction : public G4UserRunAction
{
public:
	POLYRunAction(POLYModelImport* POLYData, G4String output, G4Timer* initTimer);
	virtual ~POLYRunAction();

public:
	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

	void PrintResult(std::ostream &out);
	void PrintLine(std::ostream &out);
  
private:
	POLYModelImport* POLYData;
	POLYRun*         fRun;
	G4int            numOfEvent;
	G4int            runID;
	G4String         outputFile;
	G4Timer*         initTimer;
	G4Timer*         runTimer;
};

#endif






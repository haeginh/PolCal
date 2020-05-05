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
// POLYModelImport.hh
// \file   MRCP_GEANT4/External/include/POLYModelImport.hh
// \author Haegin Han
//

#ifndef POLYModelImport_h
#define POLYModelImport_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4TessellatedSolid.hh"
// *********************************************************************
// This class is to import the phantom data from *.ele, *.node, and
// *.material files.
// -- DataRead: Construct G4POLY by reading data from *.ele and *.node
//              files
// -- MaterialRead: Construct G4Material by reading material data from
//                  *.material file
// -- ColourRead: Construct std::map that contains G4Colour data
//                according to data in colour.dat file
// -- PrintMaterialInformation: Print a table that contains organ ID,
//                              number of POLYrahedrons, volume, density,
//                              mass, and name for each organ
// *********************************************************************

class POLYModelImport
{
public:
	POLYModelImport(G4bool isAF, G4UIExecutive* ui, G4String phantomName);
	virtual ~POLYModelImport() {};

	// get methods
	G4String      GetPhantomName()           { return phantomName; };
	G4Material*   GetMaterial(G4int idx)     { return materialMap[idx];}
	G4int         GetNumPOLYrahedron()       { return tessVec.size();}
	G4int         GetMaterialIndex(G4int idx){ return tessVec[idx].second; }
	G4TessellatedSolid* GetTess(G4int idx)   { return tessVec[idx].first; }
	std::map<G4int, std::vector<G4int>> GetTessTree() {return tree;}
	G4double      GetVolume(G4int idx)       { return volumeMap[idx]; }
	std::map<G4int, G4double> GetMassMap()   { return massMap; }
	std::map<G4int, G4Colour> GetColourMap() { return colourMap; }
	G4ThreeVector GetPhantomSize()           { return phantomSize; }
	G4ThreeVector GetPhantomBoxMin()         { return boundingBox_Min; }
	G4ThreeVector GetPhantomBoxMax()         { return boundingBox_Max; }

private:
	// private methods
	void DataRead(G4String objFileN);
	void MaterialRead(G4String);
	void ColourRead();
	void PrintMaterialInfomation();
	std::vector<std::vector<std::vector<G4int>>>
	     SeparateShell(std::vector<std::vector<G4int>> facePool);
	void ArrangeTess();

	G4String StringSplitterFirst(G4String str, G4String del) {
		size_t pos = 0;
		G4String token = str;
		while ((pos = str.find(del)) != std::string::npos) {
			token = str.substr(0, pos);
			break;
		}
		return token;
	}

	static bool TessCompare(const std::pair<G4TessellatedSolid*, G4int> &tess1, const std::pair<G4TessellatedSolid*, G4int> &tess2){
		return (tess1.first->Inside(tess2.first->GetPointOnSurface())==kInside);
	}

	std::vector<G4int> OutMostTess(std::vector<G4int> tessPool){
		std::vector<G4int> outMostTessVec;
		for(auto t1:tessPool){
			G4bool outMost(true);
			for(auto t2:tessPool){
				if(t1==t2) continue;
				if(TessCompare(tessVec[t2], tessVec[t1])){
					outMost = false;
					break;
				}
			}
			if(outMost) outMostTessVec.push_back(t1);
		}
		return outMostTessVec;
	}

	G4String phantomDataPath;
	G4String phantomName;

	G4ThreeVector boundingBox_Min;
	G4ThreeVector boundingBox_Max;
	G4ThreeVector phantomSize;

	std::vector<G4ThreeVector> vertexVector;
	std::vector<std::pair<G4TessellatedSolid*, G4int>> tessVec;
	std::map<G4int, std::vector<G4int>> tree;
	std::map<G4int, G4int>     numPOLYMap;
	std::map<G4int, G4double>  volumeMap;
	std::map<G4int, G4Colour>  colourMap;

	std::map<G4int, std::vector<std::pair<G4int, G4double>>> materialIndexMap;
	std::vector<G4int>                                       materialIndex;
	std::map<G4int, G4Material*>                             materialMap;
	std::map<G4int, G4double>                                densityMap;
	std::map<G4int, G4double>                                massMap;
	std::map<G4int, G4String>                                organNameMap;
};


#endif

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
// POLYModelImport.cc
// \file   MRCP_GEANT4/External/src/POLYModelImport.cc
// \author Haegin Han
//

#include "POLYModelImport.hh"
#include "G4UIcommand.hh"
#include "G4TriangularFacet.hh"
#include "G4Timer.hh"
#include <set>

POLYModelImport::POLYModelImport(G4bool isAF, G4UIExecutive* ui, G4String _phantomName)
:phantomName(_phantomName)
{
	// set path for phantom data
	char* pPATH = getenv("PHANTOM_PATH");
	if( pPATH == 0 ){
		// exception for the case when PHANTOM_PATH environment variable was not set
		G4Exception("POLYModelImport::POLYModelImport","",JustWarning,
				G4String("PHANTOM_PATH environment variable was not set.").c_str());
		// default path for phantom data
		phantomDataPath = ".";
	}
	else {
		// set path for phantom data as PHANTOM_PATH
		phantomDataPath = pPATH;
	}

	G4cout << "================================================================================"<<G4endl;
	G4cout << "\t" << phantomName << " was implemented in this CODE!!   "<< G4endl;
	G4cout << "================================================================================"<<G4endl;


	G4String objFile      =  phantomName + ".obj";
	G4String materialFile =  phantomName + ".material";

	// read phantom data files (*. ele, *.node)
	DataRead(objFile);
	// read material file (*.material)
	MaterialRead(materialFile);
	// read colour data file (colour.dat) if this is interactive mode
	if(ui) ColourRead();
	// print the summary of phantom information
	PrintMaterialInfomation();
}

void POLYModelImport::DataRead(G4String objFileN)
{
	std::ifstream ifs(phantomDataPath+"/"+objFileN);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.node file
		G4Exception("POLYModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + objFileN ).c_str());
	}
	G4cout << "Read " << objFileN <<"..."<< G4endl;

	G4String dump, firstStr;
	G4double x, y, z;
	std::vector<G4ThreeVector> pointVec;
	std::vector<std::vector<G4int>> faceVec;
    std::map<G4int, std::vector<std::vector<std::vector<G4int>>>> shells;
    G4int a, b, c, id(0), shellCount(0);
    G4Timer* timer = new G4Timer;
    G4double sepTime(0);
	while (getline(ifs, dump)) {
		std::stringstream ss(dump);
		if(!(ss >> firstStr)) continue;
		if (firstStr == "v") {
			ss >> x >> y >> z;
			pointVec.push_back(G4ThreeVector(x, y, z)*cm);
		}
		else if (firstStr == "f") {
			ss >> a >> b >> c;
			faceVec.push_back({ a-1, b-1, c-1});
		}
		else if (firstStr == "g") {
            ss >> dump;
            id = atoi(StringSplitterFirst(dump, "_").c_str());
            G4cout<<"Reading..ID "<<id<<" -- "<<std::flush;
            if(!faceVec.empty()){
                timer->Start();
                auto sep = SeparateShell(faceVec);
                timer->Stop();
                sepTime += timer->GetRealElapsed();
                shells[id] = sep;
				faceVec.clear();
                G4cout<<shells[id].size()<<" subshells"<<G4endl;
			}
			shellCount++;
		}
	}

    boundingBox_Max = pointVec[0];
    boundingBox_Min = pointVec[0];
    for(auto p:pointVec){
        if(p.getX()>boundingBox_Max.getX())
            boundingBox_Max.setX(p.getX());
        else if(p.getX()<boundingBox_Min.getX())
            boundingBox_Min.setX(p.getX());

        if(p.getY()>boundingBox_Max.getY())
            boundingBox_Max.setY(p.getY());
        else if(p.getY()<boundingBox_Min.getY())
            boundingBox_Min.setY(p.getY());

        if(p.getZ()>boundingBox_Max.getZ())
            boundingBox_Max.setZ(p.getZ());
        else if(p.getZ()<boundingBox_Min.getZ())
            boundingBox_Min.setZ(p.getZ());
    }
    phantomSize = boundingBox_Max-boundingBox_Min;
    G4ThreeVector center  = (boundingBox_Max + boundingBox_Min)*0.5;
    G4cout<< "Phantom shift by "<<center/cm<<" cm"<<G4endl;
    G4cout<<"************SEPATRAION TIME: "<<sepTime<<"************"<<G4endl;
    boundingBox_Max -= center;
    boundingBox_Min -= center;
    for(auto &point:pointVec) point -= center;

    for(auto shell:shells){
        for(size_t i=0;i<shell.second.size();i++){
            id = shell.first;
            G4TessellatedSolid* tess = new G4TessellatedSolid
                    (G4UIcommand::ConvertToString(id)+"_"+G4UIcommand::ConvertToString(G4int(i)));
            tessVec.push_back(std::make_pair(tess, id));
            for(auto f:shell.second[i]){
                tess->AddFacet(new G4TriangularFacet(pointVec[f[0]],pointVec[f[1]],pointVec[f[2]],ABSOLUTE));
            }
            tess->SetSolidClosed(true);
        }
    }

	G4cout<<"  => "<<tessVec.size()<<" tessellated solids were generated"<<std::flush;
	ArrangeTess();
	G4cout<<" and arranged"<<G4endl;
}

void POLYModelImport::ArrangeTess(){
	//std::sort(tessVec.begin(), tessVec.end(), TessCompare);
	tree[-1] = {};
	for(size_t i=0;i<tessVec.size();i++){
		G4bool outer(true);
		for(size_t j=0;j<tessVec.size();j++){
			if(i==j) continue;
			if(TessCompare(tessVec[j], tessVec[i])){
				outer = false;
				break;
			}
		}
		if(outer){
			tree[-1].push_back(i);
		}
	}

	for(size_t i=0;i<tessVec.size();i++){
		std::vector<G4int> innerTess;
		for(size_t j=0;j<tessVec.size();j++){
			if(i==j) continue;
			if(TessCompare(tessVec[i], tessVec[j])) innerTess.push_back(j);
		}
		tree[i] = OutMostTess(innerTess);
	}

	std::vector<G4double> volVec;
	for(size_t i=0;i<tessVec.size();i++){
		volVec.push_back(tessVec[i].first->GetCubicVolume());
	}

	for(size_t i=0;i<tessVec.size();i++){
		volumeMap[tessVec[i].second] += volVec[i];
		for(G4int sub:tree[i]){
			volumeMap[tessVec[i].second] -= volVec[sub];
		}
	}
}


std::vector<std::vector<std::vector<G4int>>>
     POLYModelImport::SeparateShell(std::vector<std::vector<G4int>> facePool)
{
	std::set<G4int> vToSearch, vToSearchNext;
	vToSearch.insert(facePool.back().begin(), facePool.back().end());
	std::vector<std::vector<G4int>> rest = facePool;

	std::vector<std::vector<G4int>> collected;
	collected.push_back(rest.back());
	rest.pop_back();

	std::vector<std::vector<G4int>> restNext;

	std::vector<std::vector<std::vector<G4int>>> separated;
	while (1) { // loop to detect all sepearated shells
		while (1) { // loop for finding all faces for each separated shell
			for (std::vector<G4int> f : rest) {
				if (vToSearch.find(f[0]) != vToSearch.end()) {
					collected.push_back(f);
					vToSearchNext.insert(f[1]);
					vToSearchNext.insert(f[2]);
				}
				else if (vToSearch.find(f[1]) != vToSearch.end()) {
					collected.push_back(f);
					vToSearchNext.insert(f[0]);
					vToSearchNext.insert(f[2]);
				}
				else if (vToSearch.find(f[2]) != vToSearch.end()) {
					collected.push_back(f);
					vToSearchNext.insert(f[0]);
					vToSearchNext.insert(f[1]);
				}
				else restNext.push_back(f);
			}
			if (vToSearchNext.size() == 0) break;

			vToSearch = vToSearchNext;
			vToSearchNext.clear();
			rest = restNext;
			restNext.clear();
		}
		separated.push_back(collected);
		if (restNext.size() == 0) break;
		rest = restNext; restNext.clear();

		vToSearch.clear();
		vToSearch.insert(rest.back().begin(), rest.back().end());
		collected.clear();
		collected.push_back(rest.back());
		rest.pop_back();
	}
	return separated;
}


void POLYModelImport::MaterialRead(G4String materialFile)
{
	// Read material file (*.material)
	//
	std::ifstream ifpMat;

	ifpMat.open((phantomDataPath + "/" + materialFile).c_str());
	if(!ifpMat.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("POLYModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + materialFile ).c_str());
	}

	G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String MaterialName;
	G4double density;

	while(!ifpMat.eof())
	{
		ifpMat >> read_data;                   //ex) 'C' RBM
		ifpMat >> MaterialName;                //ex)  C 'RBM'
		ifpMat >> read_data;
		density = std::atof(read_data);        //ex) 1.30
		ifpMat >> read_data;                   //ex) g/cm3
		ifpMat >> read_data;
		token = std::strtok(read_data,"m");
		G4int matID = std::atoi(token);        //ex) m'10'
		materialIndex.push_back(matID);
		organNameMap[matID]= MaterialName;
		densityMap[matID] = density*g/cm3;

		for(G4int i=0 ;  ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
			ifpMat >> read_data;
			fraction = -1.0 * std::atof(read_data);
			materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
		}
	}
	ifpMat.close();

	// Construct materials for each organ
	//
	   G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
	    G4NistManager* nistManager = G4NistManager::Instance();

		for(G4int i=0;i<(G4int)materialIndex.size();i++){
			G4int idx = materialIndex[i];
			G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
			for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
				if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
				else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
			}
		materialMap[idx]=mat;
		massMap[idx]=densityMap[idx]*volumeMap[idx];
	}
}

void POLYModelImport::ColourRead()
{
	// Read colour data file (colour.dat)
	//
	std::ifstream ifpColour;

	ifpColour.open((phantomDataPath + "/" + "colour.dat").c_str());
	if(!ifpColour.is_open()) {
		// exception for the case when there is no colour.dat file
		G4Exception("POLYModelImport::DataRead","",FatalErrorInArgument,
				G4String("Colour data file was not found ").c_str());
	}

//	G4cout << "  Opening colour data file 'colour.dat'" <<G4endl;

	G4int organID;
	G4double red, green, blue, alpha;
	while( ifpColour >> organID >> red >> green >> blue >> alpha )
		colourMap[organID] = G4Colour(red, green, blue, alpha);

	ifpColour.close();
}

void POLYModelImport::PrintMaterialInfomation()
{
	// Print the overall information for each organ
	//
	G4cout << G4endl
		   << std::setw(9)  << "Organ ID"
		   << std::setw(11) << "# of POLY"
		   << std::setw(11) << "vol [cm3]"
		   << std::setw(11) << "d [g/cm3]"
		   << std::setw(11) << "mass [g]"
		   << "\t" << "organ/tissue"<< G4endl ;
	G4cout << "--------------------------------------------------------------------------------"<<G4endl;

	std::map<G4int, G4Material*>::iterator matIter;
	G4cout<<std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
	{
		G4int idx = matIter->first;

		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << numPOLYMap[idx]              // # of POLYrahedrons
			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
			   << std::setw(11) << materialMap[idx]
			                       ->GetDensity()/(g/cm3)      // organ density
			   << std::setw(11) << massMap[idx]/g              // organ mass
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}
}

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
// $Id: G4CrossSectionDataStore.cc 68720 2013-04-05 09:18:58Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CrossSectionDataStore
//
// Modifications:
// 23.01.2009 V.Ivanchenko add destruction of data sets
// 29.04.2010 G.Folger     modifictaions for integer A & Z
// 14.03.2011 V.Ivanchenko fixed DumpPhysicsTable
// 15.08.2011 G.Folger, V.Ivanchenko, T.Koi, D.Wright redesign the class
// 07.03.2013 M.Maire cosmetic in DumpPhysicsTable
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4CrossSectionDataStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4HadronicException.hh"
#include "G4HadTmpUtil.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh"

#include "G4DynamicParticle.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include <iostream>

#include <time.h>

//PRUTH:  include Rob Fowler's code
#include "simplify.h"

using namespace std;


//PRUTH
G4int G4CrossSectionDataStore::objCnt = 0;
//std::unordered_map<std::pair<G4ParticleDefinition*,const G4Material*>,struct cycleCountEntry*>* G4CrossSectionDataStore::cycleCountMap = NULL;
//G4CrossSectionDataStore_CycleCountMap* G4CrossSectionDataStore::cycleCountMap = NULL;
//std::map<std::string,struct cycleCountEntry*>* G4CrossSectionDataStore::cycleCountMap = NULL;
//std::map<std::pair<G4ParticleDefinition*,const G4Material*>,struct fastPathEntry*>* G4CrossSectionDataStore::fastPathMap = NULL;
//std::map<G4int,struct fastPathEntry*>* G4CrossSectionDataStore::fastPathMap = NULL;
G4long G4CrossSectionDataStore::getCrossSectionCount = 0;
G4long G4CrossSectionDataStore::getCrossSectionCount_fastpath = 0;
G4long G4CrossSectionDataStore::getCrossSectionCount_slowpath = 0;
G4long G4CrossSectionDataStore::sampleZandACount = 0;
G4long G4CrossSectionDataStore::getCrossSectionCount_hitOneLineCache = 0;
/*
struct cycleCountEntry{
  G4DynamicParticle* particle;
  G4Material* material;

  uint64_t invocationCount;
  uint64_t totalCycles;
}

struct fastPathEntry{
G4DynamicParticle* particle;
  G4Material* material;

  G4double coefficient[9];
  G4double knot[9];
  G4double intercept;
  G4int knot_cnt;
};

*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include <iostream>
#include <fstream>
using namespace std;
#include "G4ParticleTable.hh"
G4CrossSectionDataStore::G4CrossSectionDataStore(const G4String& pname) :
  nDataSetList(0), verboseLevel(0),processName(pname)
{
  
  nist = G4NistManager::Instance();
  currentMaterial = elmMaterial = 0;
  currentElement = 0;  //ALB 14-Aug-2012 Coverity fix.
  matParticle = elmParticle = 0;
  matKinEnergy = elmKinEnergy = matCrossSection = elmCrossSection = 0.0;
  
  G4cout << "pruth: G4CrossSectionDataStore: " << pname << G4endl;

  G4ParticleDefinition *partDef;
  G4Material* mat;
  if(objCnt <= 0){
    G4cout << "pruth: Found first G4CrossSectionDataStore" << G4endl;

    partDef = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    if (partDef != NULL){
      G4cout << "pruth: partDef : " << partDef->GetParticleName()  << G4endl;
      
    } else {
      G4cout << "pruth: partDef is NULL" << G4endl;
    }
    
    /*
    sampleZandACount = 0;
    getCrossSectionCount = 0;
    getCrossSectionCount_fastpath = 0;
    getCrossSectionCount_slowpath = 0;
    getCrossSectionCount_hitOneLineCache = 0;
    */
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nmat = theMaterialTable->size();
  size_t i;
  for(i = 0; i < nmat; i++){
    mat = (*theMaterialTable)[i];
    //G4cout << "material (" << i << "): " << mat->GetName() << G4endl;                                                                                                           
  }

  objID = objCnt;
  //open results file if needed                                                                                                                                              
  //ofstream resultsFile ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel/RAW_XC_DATA/"  + std::to_string(objID) + "_" + processName + "_RESULTS");
  resultsFile.open("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/RAW_XC_DATA/"  + std::to_string(objID) + "_" + processName + "_RESULTS");

  //Sample and Surragate model vars
  queryMax = 10000;
  sampleMin = 0.0001;
  sampleMax = 10000;
  sampleCount = 200000;
  dpTol = 0.01;

  std::string   line;
  ifstream samplingConfigfile ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/samplingConfig.txt");
  if(samplingConfigfile.is_open()){
    G4cout << "confige file open" << G4endl;
    while (getline(samplingConfigfile,line)){
      G4cout << "pruth: line: " << line << G4endl;                                                                                                                                                                                  

      std::stringstream lineStream(line);
      std::string queryMaxStr;
      std::string sampleMinStr;
      std::string sampleMaxStr;
      std::string sampleCountStr;
      std::string dpTolStr;

      lineStream >> queryMaxStr;
      
      if(!queryMaxStr.compare(0,1,"#")){
	G4cout << "Found comment: " << line;
	continue;
      }

      G4cout << "Setting sampling config values "<< G4endl;

      lineStream >> sampleMinStr;
      lineStream >> sampleMaxStr;
      lineStream >> sampleCountStr;
      lineStream >> dpTolStr;

      

      queryMax = std::stod(queryMaxStr);
      sampleMin = std::stod(sampleMinStr);
      sampleMax = std::stod(sampleMaxStr);
      sampleCount = std::stoi(sampleCountStr);
      dpTol = std::stod(dpTolStr);

      break;
    }
  } 
  samplingConfigfile.close();
  
  G4cout << "queryMax: " << queryMax << G4endl;
  G4cout << "sampleMin: "<<  sampleMin <<  G4endl;
  G4cout << "sampleMax: "<< sampleMax<<  G4endl;
  G4cout << "sampleCount: "<< sampleCount <<  G4endl;
  G4cout << "dpTol: "<< dpTol <<  G4endl;



  //Read interseting triples
  //G4String line;
  //std::string   line;
  ifstream file ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/Triples.txt");
  if(file.is_open()){
    while (getline(file,line)){
      //G4cout << "pruth: line: " << line << G4endl;
      
      std::stringstream lineStream(line);
      std::string process;
      std::string particle;
      std::string material;
      std::string cutoff;
      
      lineStream >> process;
      lineStream >> particle;
      lineStream >> material;
      lineStream >> cutoff;
      //G4cout << "READ Triple: "  << process << ", " << particle << ", " << material << G4endl;
      
      if(processName.compare(process) == 0){
	theMaterialTable = G4Material::GetMaterialTable();
	nmat = theMaterialTable->size();
	for(i = 0; i < nmat; i++){
	  mat = (*theMaterialTable)[i];
	  //G4cout << "material (" << i << "): " << mat->GetName() << G4endl;                                                                                                                  
	  if (mat->GetName().compare(material) == 0) break;
	}
	//G4cout << "found material (" << i << "): " << mat->GetName() << G4endl;
	
	std::pair<G4ParticleDefinition*,const G4Material*> searchkey = std::make_pair(G4ParticleTable::GetParticleTable()->FindParticle(particle),mat);
	
	
	//struct fastPathEntry* entry =  (fastPathMap)[searchkey] ;
	struct fastPathEntry* entry = new struct fastPathEntry();
	
	entry->particle =  G4ParticleTable::GetParticleTable()->FindParticle(particle);
	entry->material =  mat;
	entry->min_cutoff = std::stod(cutoff);
	entry->physicsVector = NULL;
       

	//fastPathEntry debug stats
	/*
	entry->count = 0;
	entry->slowpath_sum = 0.0;
	entry->max_delta = 0.0;
	entry->min_delta = 0.0;
	entry->sum_delta = 0.0;
	entry->sum_delta_square= 0.0;
	*/

	G4cout << "ADDING Triple: "  << process << ", " << particle << ", " << material << G4endl;
	//(fastPathMap)[searchkey] = entry;
	
	//add entry                                                                                                                                                                              
	struct cycleCountEntry* cycleCountEntry; // = (cycleCountMap)[searchkey];
        cycleCountEntry = new struct cycleCountEntry();
	//G4cout << "PRUTH: 1101 GetCrossSection for matParticle: " << matParticle <<  G4endl;                                                                                                   
	cycleCountEntry->particle = particle; //part;                                                                                                                                
	cycleCountEntry->material = mat;
	cycleCountEntry->dataset = "unused"; //dataSetList[nDataSetList-1]->GetName();

	//comment out for speed
	cycleCountEntry->initCyclesFastPath = 0;
	cycleCountEntry->invocationCountSlowPath=0;
	cycleCountEntry->totalCyclesSlowPath=0;
	cycleCountEntry->invocationCountFastPath=0;
	cycleCountEntry->totalCyclesFastPath=0.0;
	cycleCountEntry->invocationCountOneLineCache=0;
	

	cycleCountEntry->fastPath = entry;

	cycleCountEntry->energy = -1.0;
	cycleCountEntry->crossSection = 0.0;
	cycleCountEntry->cacheHitCount = 0;
	(cycleCountMap)[searchkey] = cycleCountEntry;


      }
      

    }
    file.close();
  }

  objCnt++;
  G4cout << "pruth: End G4CrossSectionDataStore Constructor> " <<  objCnt <<  G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::~G4CrossSectionDataStore() 
{
  objCnt--;
  G4cout << "pruth: G4CrossSectionDataStore Destructor> " << processName << ", " <<  objCnt <<  G4endl;

  
  if(objCnt == 0){
    G4cout << "pruth: Cleaning up last G4CrossSectionDataStore"<< G4endl;
    G4cout << "pruth: sampleZandACount: " << sampleZandACount << G4endl;
    G4cout << "pruth: getCrossSectionCount: " << getCrossSectionCount << G4endl;
    G4cout << "pruth: getCrossSectionCount_fastpath: " << getCrossSectionCount_fastpath << G4endl;
    G4cout << "pruth: getCrossSectionCount_slowpath: " << getCrossSectionCount_slowpath << G4endl;
    G4cout << "pruth: getCrossSectionCount_hitOneLineCache: " << getCrossSectionCount_hitOneLineCache << G4endl;

  }

    //Read interseting triples                                                                                                                                                                                        
    //G4String line;                                                                                                                                                                                               
  G4ParticleDefinition *partDef;
  G4Material* mat;
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    size_t nmat = theMaterialTable->size();
    size_t i;
    for(i = 0; i < nmat; i++){
      mat = (*theMaterialTable)[i];
      //G4cout << "material (" << i << "): " << mat->GetName() << G4endl;                                                                                                                                             
    }

    std::string   line;
    ifstream file ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/Triples.txt");
    if(file.is_open()){
      while (getline(file,line)){

	std::stringstream lineStream(line);
	std::string process;
	std::string particle;
	std::string material;

	lineStream >> process;
	lineStream >> particle;
	lineStream >> material;

	//G4cout << "Process: " << processName <<  ", Writing Output for Triple: "  << process << ", " << particle << ", " << material << G4endl;

	if(processName.compare(process) == 0){
	  theMaterialTable = G4Material::GetMaterialTable();
	  nmat = theMaterialTable->size();
	  for(i = 0; i < nmat; i++){
	    mat = (*theMaterialTable)[i];
	    //G4cout << "material (" << i << "): " << mat->GetName() << G4endl;                                                                                                                                       
	    if (mat->GetName().compare(material) == 0) break;
	  }
	  //G4cout << "found material (" << i << "): " << mat->GetName() << G4endl;                                                                                                                                   
	  std::pair<G4ParticleDefinition*,const G4Material*> searchkey = std::make_pair(G4ParticleTable::GetParticleTable()->FindParticle(particle),mat);

	  //struct fastPathEntry* entry =  (fastPathMap)[searchkey] ;                                                                                                

    
	  struct cycleCountEntry* countMap = (cycleCountMap)[searchkey];
	  struct fastPathEntry* entry = countMap->fastPath;

	  if(entry != NULL){
	    G4cout << "PRUTH: TRIPLE SUMMARY: ";
	    G4cout << ", processNum " << objCnt;  
            G4cout << ", processName " << processName;
            G4cout << ", particle " << particle;
	    G4cout << ", material " <<  mat->GetName();
	    
	    //fastPathEntry debug stats 
	    G4cout << ", initCyclesFastPath " << countMap->initCyclesFastPath; 

	    G4cout << ", invocationCountSlowPath " << countMap->invocationCountSlowPath;
	    G4cout << ", totalCyclesSlowPath " << countMap->totalCyclesSlowPath;
	    
	    G4cout << ", invocationCountFastPath " << countMap->invocationCountFastPath;
	    G4cout << ", totalCyclesFastPath " << countMap->totalCyclesFastPath;

	    G4cout << ", invocationCountTriedOneLineCache " << countMap->invocationCountTriedOneLineCache;
	    G4cout << ", invocationCountOneLineCache " <<  countMap->invocationCountOneLineCache;
	    
                                                                                                               
	    G4cout << ", count " << entry->count;
	    G4cout << ", slowpath_sum " << entry->slowpath_sum;
	    G4cout << ", max_delta " << entry->max_delta;
	    G4cout << ", min_delta " << entry->min_delta;
	    G4cout << ", sum_delta " << entry->sum_delta;
	    G4cout << ", sum_delta_square " << entry->sum_delta_square;
	    
	    G4cout << G4endl;
	  }
	}
    //Clean up memory for maps...
      }
      
    }      
    resultsFile.close();
   

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
#include <inttypes.h>
static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}


//PRUTH
/*
void G4CrossSectionDataStore::RemoveBias(vector<G4double> data1, vector<G4double> data2, vector<G4double>& result){
  
  //Create original data arrays
  G4int originalSize = data1.size()/2;
  vector<G4double> originalx; //[originalSize];
  vector<G4double> originaly; //[originalSize];
  int j = 0;
  
  for (int i = 0; i < data1.size(); i += 2, j++){
    originalx.push_back(data1[i]);
    originaly.push_back(data1[i+1]);
  }

  //Create new data arrays
  G4int newSize =  data2.size()/2;
  G4double newx [newSize];
  G4double newy [newSize];
  j = 0;
  for (int i = 0; i < data2.size(); i += 2, j++){
    newx[j] = data2[i];
    newy[j] = data2[i+1];
  }

  //Create index mapping array
  G4int xindex[newSize];
  j = 0;
  for (int k = 0; k <newSize; k++) {
    for (int i = 0; i < originalSize; i++) {
      if (originalx[i] == newx[k]) {
	xindex[j++] = i;
      }
    }
  }

  //Convert to Sarah's names
  G4int m = newSize;

  G4double GArea [m-1];
  G4double GAreatotal = 0;
  G4double GAreatemp = 0;
  G4double trap;
  //Area of oversampled curve                                                                                                                                
  for(int i = 0; i < m-1; i++){
    for(j = xindex[i]; j< xindex[i+1]; j++){
      trap = (originaly[j+1]+originaly[j])*(originalx[j+1]-originalx[j])/2.0;
      GAreatemp = GAreatemp + trap;
    }
    GArea[i] = GAreatemp;
    GAreatotal = GAreatotal + GAreatemp;
    GAreatemp = 0;
  }

  //aleph                                                                                                                                                    
  G4double aleph [m-1];
  for(int i = 0; i< m-1; i++){
    aleph[i] = (newx[i+1]-newx[i])/2.0;
  }
  
  //solve for f                                                                                                                                              
  G4double adjustedy [m];
  adjustedy[m-1] = newy[m-1];
  for(int i = 2; i < m+1; i++) {
    adjustedy[m-i] = (GArea[m-i]/aleph[m-i]) - adjustedy[m-i+1];
  }
  
  //error and difference tracking                                                                                                                            
  G4double difference [m];
  G4double maxdiff = 0;
  G4double adjustedarea = 0;
  G4double newarea = 0;
  for(int i = 0; i < m-1; i++){
    trap = (adjustedy[i+1]+adjustedy[i])*(newx[i+1]-newx[i])/2.0;
    adjustedarea = adjustedarea+trap;
    trap = (newy[i+1]+newy[i])*(newx[i+1]-newx[i])/2.0;
    newarea = newarea + trap;
  }
  for(int i = 0; i <m; i++) {
    difference[i] = newy[i]-adjustedy[i];
  }
  for(int i = 0; i <m; i++){
    if(abs(difference[i]) > maxdiff) {
      maxdiff = abs(difference[i]);
    }
  }

  //vector<G4double>* unbiasData = new vector<G4double>(); 
  for(int i = 0; i < newSize; i++){
    result.push_back(newx[i]);
    result.push_back(adjustedy[i]);
  }

  //return unbiasData;
}
*/

G4double G4CrossSectionDataStore::SampleCrossSectionValue(const G4DynamicParticle* part, const G4Material* mat, G4double xval)
{
  G4int i, j;
  G4double acc;
  G4double ret_value = 0.0;
  G4ParticleDefinition* matParticleTemp;
  G4double matKinEnergyTemp;

  G4Material* savedMaterial = (G4Material*)mat;
  G4double savedEnergy = part->GetKineticEnergy();

  ((G4DynamicParticle*)(part))->SetKineticEnergy(xval);
  // G4cout << " (3) Energy " << part->GetKineticEnergy() << G4endl; 

  matParticle   = matParticleTemp;
  matKinEnergy  = matKinEnergyTemp;
  currentMaterial = mat;

  //matCrossSection = 0;

  const G4int nElements = mat->GetNumberOfElements();
  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

  //if(G4int(xsecelm.size()) < nElements) {
  //  xsecelm.resize(nElements);
  //}
  acc = 0.0;
  for(i=0; i < nElements; ++i) {
    ret_value = GetCrossSection(part, (*mat->GetElementVector())[i] , mat);

    acc += nAtomsPerVolume[i] * ret_value;
    //xsecelm[i] = acc;
  }

  //for(i=0; i < nElements; ++i) {
  //  xsecelm[i] = 0.0;
  //}
  
  ((G4DynamicParticle*)(part))->SetKineticEnergy(savedEnergy);
  matKinEnergyTemp = savedEnergy;
  currentMaterial = savedMaterial;
  //matCrossSection = 0;

  return acc;
}




//#include "psimpl.h"

void G4CrossSectionDataStore::writeXSVector(G4String label, vector <G4double> vect){
  //comment out to write
  //return;

  ofstream myfile ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/RAW_XC_DATA/"  + std::to_string(objID) + "_" + processName + "_" + label);

  if (myfile.is_open())
    {
      myfile << "#" << label << G4endl;
      
      for(int i = 0; i < vect.size(); i+=2){
	myfile << vect[i] << " " << vect[i+1] << G4endl;
      }
      myfile.close();
    }
  else{
    G4cout << "PRUTH: Could not open file for writing vector: " << label <<  G4endl;   
  }
}


void G4CrossSectionDataStore::writeLinearXSSample(G4String label, G4double min, G4double max, G4double step, const G4DynamicParticle* part, const G4Material* mat){
  ofstream myfile ("/hpc_shared/home/pruth/GEANT-4/run-blanco-geant4.10.00_SurrogateModel_v2/RAW_XC_DATA/"  + std::to_string(objID) + "_" + processName + "_" + label);


  G4double xs;
  if (myfile.is_open())
    {
      myfile << "#" << label << G4endl;

      G4double currEnergy = 0.0;
      G4int i = 0;
      for(currEnergy = min; currEnergy < max; currEnergy += step){
	xs = SampleCrossSectionValue(part, mat, currEnergy);

	if(xs != 0.0){
	  myfile << currEnergy << " " << xs << G4endl;
	}

	if(i++%1000000 == 0){ 
	  G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl;
	}
      } // --- end of loop i  


      myfile.close();
    }
  else{
    G4cout << "PRUTH: Could not open file for writing vector: " << label <<  G4endl;
  }


}


G4PhysicsVector* G4CrossSectionDataStore::InitializeFastPathPhysicsVector(const G4DynamicParticle* part, const G4Material* mat, G4double cutoff){
  //G4cout << "PRUTH: Begin Initializing Physics Vector for process: " << processName <<  G4endl; 
 
  using namespace std;

  vector <Point> data_in;
  vector <Point> decimated_data;
  vector <Point> debiased_data;

  /*
  G4double max_query = 1000000;
  G4double min = 0.0001;
  G4double max = 10000;
  G4int count = 200000;
  G4double tol = 0.00001;
  */

  G4double max_query = queryMax;
  G4double min = sampleMin;
  G4double max = sampleMax;
  G4int count = sampleCount;
  G4double tol = dpTol;


  G4double xs;

  /*
  G4double linear_min = 0.0;
  G4double linear_max = 21;
  G4double linear_step = 0.001;
  */
  //Write a linear sample for debug
  //writeLinearXSSample("LINEAR_RAW_" + part->GetDefinition()->GetParticleName() + "_" +  mat->GetName(), linear_min, linear_max, linear_step, part, mat);
 
 
  /*
  G4double currEnergy;
  G4int i = 0;
  for(currEnergy = linear_min; currEnergy < linear_max; currEnergy += linear_step){
    xs = SampleCrossSectionValue(part, mat, currEnergy);
    //G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl;                                                                                                     
    polyline.push_back(currEnergy);
    polyline.push_back(xs);
   
    //if(xs > 0){
    //  G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl;
    //}
 
    if(i++%1000000 == 0){
      G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl; 
    }
  } // --- end of loop i                                                                                                                                                                   
  xs = SampleCrossSectionValue(part, mat, max);
  polyline.push_back(max);
  polyline.push_back(xs);

  writeXSVector("LINEAR-RESULTS_RAW-" + processName + "-" + part->GetDefinition()->GetParticleName() + "-" +  mat->GetName(), polyline);
  */
  
  //polyline.clear();
  


  
  //Shift so max and min are >= 1. 
  //Don't forget to shift back before computing XS
  G4double shift = 0.0;
  if(min < 1.0){
    shift = 1.0 - min;
  }
  min += shift;
  max += shift;

  G4double log_max = log10(max);
  G4double log_min = log10(min);
  G4double log_step = (log_max-log_min)/(1.0*count);
  
  G4int resultCount=10;
  
  G4double max_xs = 0.0;

  //G4cout << "PRUTH: energy value " << ", max: " << max << ", min: " << min << ", count: " << count <<  G4endl;
  //G4cout << "PRUTH: energy value " << ", log_max: " << log_max << ", log_min: " << log_min <<  ", log_step: " << log_step <<  G4endl;
  
  //add the cutoff energy  
  xs = SampleCrossSectionValue(part, mat, cutoff); 
  data_in.push_back({max_query,xs});  

  G4double currEnergy = 0.0;
  G4int i = 0;
  //log results
  for(G4double log_currEnergy = log_min; log_currEnergy < log_max; log_currEnergy += log_step){
    currEnergy = exp10(log_currEnergy) - shift;
    if (currEnergy <  cutoff) continue;

    xs = SampleCrossSectionValue(part, mat, currEnergy);
    //G4cout << "PRUTH: energy value " << currEnergy << ", XS value " << xs << G4endl;
    if (xs > max_xs) max_xs = xs;
    data_in.push_back({currEnergy,xs});
    i++;
  } // --- end of loop i
  xs = SampleCrossSectionValue(part, mat, max - shift);
  data_in.push_back({max-shift,xs});
  
  G4cout << "PRUTH: init fast path took " << i << " samples. size of data_in = "  << data_in.size() << G4endl;  

  //add the highest possible query
  //xs = SampleCrossSectionValue(part, mat, max_query);
  //data_in.push_back({max_query,xs});                

  
  //writeXSVector("RAW_" + part->GetDefinition()->GetParticleName() + "_" +  mat->GetName(), decimated_data);

  // NOT THIS ONE ///psimpl::simplify_douglas_peucker_n <2> (polyline.begin(), polyline.end(), resultCount, back_inserter(result));
  //psimpl::simplify_douglas_peucker <2> (polyline.begin(), polyline.end(), tol, back_inserter(result));
  //num_points_simplified = simplify_function(tolerance,  data_in,  decimated_data);
  tol = max_xs * 0.01;
  simplify_function(tol,  data_in,  decimated_data);
 
  //G4cout << "PRUTH: BEFORE zero removal result size: " << result.size() << G4endl; 

  //test to write output
  //writeXSVector("LOG_DP_RAW_" + part->GetDefinition()->GetParticleName() + "_" +  mat->GetName(), result);

  
  /*
  G4int zeroCount = 0;
  for(i = 0; i < result.size(); i+=2){
    if(result[i+1] == 0.0){
      zeroCount++;
    } else {
      zeroCount = 0;
    }

    //if series of zeros then erase the middle ones
    if(zeroCount >= 3) result.erase(result.begin()+(i-2),result.begin()+(i));
  }
 
  G4cout << "PRUTH: AFTER zero removal result size: " << result.size() << G4endl;
  */  

  /*
  writeXSVector("LOG_DP_NOZEROS_" + part->GetDefinition()->GetParticleName() + "_" +  mat->GetName(), result);
  */

  /*
  for(i = 0; i < result.size(); i+=2){
    if (result[i+1] < 0){
      G4cout << "PRUTH: i: NEGATIVE result: " << i/2 << ", e: " << result[i] << ", xs: " << result[i+1] << G4endl;
    }
  }
  */

  RemoveBias( data_in,  decimated_data,  debiased_data);



  //writeXSVector("LOG_DP_UNBIAS_" + part->GetDefinition()->GetParticleName() + "_" +  mat->GetName(), unbiasData);

  /*
  for(i = 0; i < unbiasData.size(); i+=2){
    if (unbiasData[i+1] < 0){
      G4cout << "PRUTH: i: NEGATIVE unbiasData: " << i/2 << ", e: " << unbiasData[i] << ", xs: " << unbiasData[i+1] << G4endl;
    }
  

  //  G4cout << "PRUTH: i: unbiasData: " << i/2 << ", e: " << unbiasData[i] << ", xs: " << unbiasData[i+1]; // <<  G4endl;
  //  G4cout << ",  results: "  << ", e: " << result[i] << ", xs: " << result[i+1];  
  //  G4cout << ",  diff: " << unbiasData[i+1] - result[i+1] << G4endl;
  }
  */
  
  //G4cout << "physics vector size: " << unbiasData.size()/2 << G4endl;  
  
  
  G4PhysicsFreeVector *physicsVector = new G4PhysicsFreeVector(decimated_data.size());
  G4int physicsVectorIndex = 0;
  for(i = 0; i < decimated_data.size(); i++){
    physicsVector->PutValue(physicsVectorIndex++, decimated_data[i].e, decimated_data[i].xs);
  }
  

  /*
  G4PhysicsFreeVector *physicsVector = new G4PhysicsFreeVector(polyline.size()/2);
  G4int physicsVectorIndex = 0;
  for(i = 0; i < polyline.size(); i+=2){
    physicsVector->PutValue(physicsVectorIndex++, polyline[i], polyline[i+1]);
  }
  */
  
  /*
  G4PhysicsFreeVector *physicsVector = new G4PhysicsFreeVector(result.size()/2);
  G4int physicsVectorIndex = 0;
  for(i = 0; i < result.size(); i+=2){
    physicsVector->PutValue(physicsVectorIndex++, result[i], result[i+1]);
  }
  */
  
  

  

  //G4cout << "PRUTH: End Initializing Physics Vector size: " << physicsVectorIndex <<  G4endl;
  return (G4PhysicsVector*) physicsVector;
}


G4double G4CrossSectionDataStore::GetCrossSectonFastPath(struct fastPathEntry* fast_entry, const G4DynamicParticle* part){
  
  if(fast_entry->physicsVector == NULL){
    //G4cout << "PRUTH: Init FastPath:  "  << processName << ", " << part->GetDefinition()->GetParticleName() << ", " << fast_entry->material->GetName() << G4endl;      
    fast_entry->physicsVector = InitializeFastPathPhysicsVector(part,fast_entry->material, fast_entry->min_cutoff);
    //G4cout << "PRUTH: AFTER Init FastPath:  "  << processName << ", " << part->GetDefinition()->GetParticleName() << ", " << fast_entry->material->GetName() << ", physicsVector " << fast_entry->physicsVector << G4endl;
  }

 
  //With physics vector
  //G4cout << "FastPath used BEFORE get Value: " << processName << ". energy: " << part->GetKineticEnergy() << G4endl;
  G4double retVal_pv = fast_entry->physicsVector->Value(part->GetKineticEnergy());
  
  //prevent negative cross sections
  if(retVal_pv < 0){
    retVal_pv = 0.0;
  }
  

  //G4double realVal = SampleCrossSectionValue(part, fast_entry->material, part->GetKineticEnergy() );
  //Log results compared to "correct" result
  //if(abs(realVal - retVal_pv) > 0.00001 && abs(realVal - retVal_pv) > abs(realVal * 0.05)){
  //  resultsFile  << "PRUTH DELTA GT 5%: " << processName << " " << part->GetDefinition()->GetParticleName() << " " << fast_entry->material->GetName()  << " " << part->GetKineticEnergy()  << " " << realVal  << " " << retVal_pv << " " <<  (realVal-retVal_pv)<< " " << (abs(retVal_pv - realVal)/realVal)*100 << "% "  << G4endl;
  //} 
  //else {
  //  G4cout << "PRUTH DELTA LT 5%: " << processName << " " << part->GetDefinition()->GetParticleName() << " " << fast_entry->material->GetName()  << " " << part->GetKineticEnergy()   << " " << realVal  <<      " " << retVal_pv << " " << (realVal-retVal_pv)<< " " << (abs(retVal_pv - realVal)/realVal)*100 << "% "  << G4endl;
  //}



  //IF DEBUG
  /*
    G4cout << "PRUTH: particle " << part->GetDefinition()->GetParticleName();
    G4cout << ", material " <<  fast_entry->material->GetName();

    G4cout << ", count " << fast_entry->count;
    G4cout << ", slowpath_sum " << fast_entry->slowpath_sum;
    G4cout << ", max_delta " << fast_entry->max_delta;
    G4cout << ", min_delta " << fast_entry->min_delta;
    G4cout << ", sum_delta " << fast_entry->sum_delta;
    G4cout << ", sum_delta_square " << fast_entry->sum_delta_square;
    G4cout << G4endl;
  */

  
  
  fast_entry->count++;
  
  /*
  fast_entry->slowpath_sum += realVal;
  
  G4double delta = retVal_pv-realVal;
  
  if(delta > fast_entry->max_delta)  
    fast_entry->max_delta = delta;

  if(delta < fast_entry->min_delta) 
    fast_entry->min_delta = delta;

  fast_entry->sum_delta += abs(delta);
  fast_entry->sum_delta_square += abs(delta)*abs(delta);
  */


  
  

  //return fast path value  
  return retVal_pv;

  //return slow path value  
  //return realVal;

}


G4bool prevCalcUsedFastPath=false;

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Material* mat)
{
  //set 3rd arg to true if you want to force slowpath for all calls
  return GetCrossSection(part, mat,false);
}

G4double 
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Material* mat, G4bool requireSlowPath){
  getCrossSectionCount++;

  uint64_t t_start_total;
  uint64_t t_end_total;
  uint64_t t_diff_total;
  
  //G4cout << "PRUTH: BEGIN GetCrossSection for process: " << processName << ",  particle: " << part->GetDefinition()->GetParticleName() << ", material: " << mat->GetName() <<  ", energy: " << part->GetKineticEnergy() <<  G4endl;


  t_start_total=rdtsc();

  G4ParticleDefinition* tmp_matParticle = part->GetDefinition();
  //matParticle = part->GetDefinition();
  pair<G4ParticleDefinition*,const G4Material*> searchkey = std::make_pair(tmp_matParticle,mat);
  //std::string entryName= tmp_matParticle->GetParticleName() + ":" + mat->GetName() + ":" + dataSetList[nDataSetList-1]->GetName();
  struct cycleCountEntry* entry = (cycleCountMap)[searchkey];

  //actually might be a bug... i think that if the dataset is different this gives a wrong answer
  if(mat == currentMaterial && part->GetDefinition() == matParticle
     && part->GetKineticEnergy() == matKinEnergy) 
    {
      if(entry!=NULL){  entry->invocationCountTriedOneLineCache++; }
      if(!prevCalcUsedFastPath && !requireSlowPath){
	getCrossSectionCount_hitOneLineCache++;

	if(entry!=NULL){ entry->invocationCountOneLineCache++; }

	return matCrossSection; 
      } else {
	requireSlowPath=true;
      }
    }

  if(entry != NULL && entry->energy == part->GetKineticEnergy()){
    entry->cacheHitCount++;
    if(!requireSlowPath){
      return entry->crossSection;
    }
  }

  //G4int searchkey=0;
  //if(entryName == "neutron:materials_StainlessSteel:G4NeutronInelasticXS"){
  //  searchkey=1;
  //}

  //pair<G4ParticleDefinition*,const G4Material*> searchkey = std::make_pair(tmp_matParticle,mat);
  
  
  //G4cout << "pruth: GetCrossSectionDataStore (1000) " << G4endl; 
  

  //t_start_total=rdtsc();

  currentMaterial = mat;
  matParticle = part->GetDefinition();
  matKinEnergy = part->GetKineticEnergy();
  matCrossSection = 0;
  
  //PRUTH map entry key
  //entryName= matParticle->GetParticleName() + ":" + mat->GetName() + ":" + dataSetList[nDataSetList-1]->GetName();
  struct fastPathEntry* fast_entry = NULL;
  if(!requireSlowPath && entry != NULL){
    //fast_entry = (fastPathMap)[searchkey];
    fast_entry = entry->fastPath;
  } 

  if(fast_entry != NULL && part->GetKineticEnergy() < fast_entry->min_cutoff){
    requireSlowPath = true;
  }

  if (!requireSlowPath && fast_entry != NULL){
    //G4cout << "pruth: Fast Path Found: " << entryName << G4endl; 
    getCrossSectionCount_fastpath++;
    matCrossSection = GetCrossSectonFastPath(fast_entry,part);
    prevCalcUsedFastPath = true;

  
  } else {
    
    getCrossSectionCount_slowpath++;
    prevCalcUsedFastPath=false;

    //G4cout << "=" << elmParticle->GetParticleName() << ":" << mat->GetName();

    G4int nElements = mat->GetNumberOfElements();
    const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();
    
    if(G4int(xsecelm.size()) < nElements) { xsecelm.resize(nElements); }
    
    for(G4int i=0; i<nElements; ++i) {
      matCrossSection += nAtomsPerVolume[i] * 
	GetCrossSection(part, (*mat->GetElementVector())[i], mat);
      xsecelm[i] = matCrossSection;
    }
    
  }
  
  t_end_total=rdtsc();
  t_diff_total=t_end_total-t_start_total;

  /*  
  if (fast_entry != NULL){
    G4cout << "pruth: Fast Path Result: " << entryName << " -- ";
    G4cout << "energy: " << part->GetKineticEnergy();
    G4cout << ",fastpath: " <<  fastPathXS;
    G4cout << ", slowpath: " <<  matCrossSection;
    G4cout << ", diff: " <<  fastPathXS - matCrossSection;
    
    G4double diff = 0.0;
    if(fastPathXS > matCrossSection)
    diff = fastPathXS - matCrossSection;
    else
      diff = matCrossSection - fastPathXS;

    G4double percentDiff=diff/matCrossSection;
    G4cout << ", percentDiff: " << percentDiff;
    if(percentDiff < 0.0001)
      G4cout << "XXXXX"; 

    G4cout << G4endl;
  }
  */

  
  //G4cout << "pruth: GetCrossSectionDataStore (2000) " << G4endl;

  // G4cout << "PRUTH: 1000 GetCrossSection for process: " << processName <<  G4endl;
  //G4cout << ":" << t_diff << G4endl;
  //std::string entryName= matParticle->GetParticleName() + ":" + mat->GetName();
  //struct cycleCountEntry* entry = (*cycleCountMap)[entryName];
  matParticle = part->GetDefinition();
  if (entry == NULL){
    //add entry
    //matParticle = part->GetDefinition();

    //G4cout << "pruth: Adding: " << entryName << G4endl; 
    entry = new struct cycleCountEntry();
    //G4cout << "PRUTH: 1101 GetCrossSection for matParticle: " << matParticle <<  G4endl;
    entry->particle = matParticle->GetParticleName(); //part;
    entry->material = mat;
    entry->dataset = dataSetList[nDataSetList-1]->GetName();
     
    entry->invocationCountSlowPath=0;
    entry->totalCyclesSlowPath=0;
    
    entry->invocationCountFastPath=0;
    entry->totalCyclesFastPath=0.0;
    
    entry->invocationCountOneLineCache=0;
    
    entry->fastPath = NULL;

    entry->energy = 0.0;
    entry->crossSection = 0.0;
    entry->cacheHitCount = 0;
    (cycleCountMap)[searchkey] = entry;
  } 

  entry->energy = part->GetKineticEnergy();
  entry->crossSection = matCrossSection;
  
  if(fast_entry != NULL){
    if(entry->invocationCountFastPath == 0){
      G4cout << "PRUTH: initCyclesFastPath = " << t_diff_total << G4endl;
      entry->initCyclesFastPath = t_diff_total;
      entry->invocationCountFastPath++;
    } else {
      //the first one includes the initialization... don't count it for now
      entry->invocationCountFastPath++;
      entry->totalCyclesFastPath+=t_diff_total;
    }
  } else {
    entry->invocationCountSlowPath++;
    entry->totalCyclesSlowPath+=t_diff_total;
  }
  
  //G4cout << "pruth: END GetCrossSectionDataStore (end): CrossSection:  " << matCrossSection << G4endl;
  
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  if(mat == elmMaterial && elm == currentElement &&
     part->GetDefinition() == elmParticle &&
     part->GetKineticEnergy() == elmKinEnergy) 
    { return elmCrossSection; }

  elmMaterial = mat;
  currentElement = elm;
  elmParticle = part->GetDefinition();
  elmKinEnergy = part->GetKineticEnergy();
  elmCrossSection = 0.0;

  G4int i = nDataSetList-1;  
  G4int Z = G4lrint(elm->GetZ());
  if (elm->GetNaturalAbundanceFlag() &&
      dataSetList[i]->IsElementApplicable(part, Z, mat)) {

    // element wise cross section
    //PRUTH
    /*
    std::string entryName= matParticle->GetParticleName() + ":" + mat->GetName();
    struct fastPathEntry* fast_entry = (*fastPathMap)[entryName];
    if (fast_entry != NULL){
      G4cout << "pruth: Fast Path: Element wise cross section: Z = " << Z << ", i = " << i << ", dataSetList[i]->GetName() = " << dataSetList[i]->GetName()  << G4endl ;
    }
    */
    //

    elmCrossSection = dataSetList[i]->GetElementCrossSection(part, Z, mat);

    //G4cout << "Element wise " << elmParticle->GetParticleName() 
    //	   << " xsec(barn)= " <<  elmCrossSection/barn 
    //	   << "  E(MeV)= " << elmKinEnergy/MeV 
    //	   << " Z= " << Z << " AbundFlag= " << elm->GetNaturalAbandancesFlag()
    //	   <<G4endl;

  } else {
    // isotope wise cross section
    G4int nIso = elm->GetNumberOfIsotopes();    
    G4Isotope* iso = 0;

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = elm->GetIsotopeVector();
    G4double* abundVector = elm->GetRelativeAbundanceVector();

    //PRUTH
    /*
    std::string entryName= matParticle->GetParticleName() + ":" + mat->GetName();
    struct fastPathEntry* fast_entry = (*fastPathMap)[entryName];
    if (fast_entry != NULL){                                         
      G4cout << "pruth: Fast Path: Isotope wise cross section: Z = " << Z << ", i = " << i  << ", dataSetList[i]->GetName() = " << dataSetList[i]->GetName()  << G4endl;
    }
    */
    //

    for (G4int j = 0; j<nIso; ++j) {
      if(abundVector[j] > 0.0) {
	iso = (*isoVector)[j];
	elmCrossSection += abundVector[j]*
	  GetIsoCrossSection(part, Z, iso->GetN(), iso, elm, mat, i);
	//G4cout << "Isotope wise " << elmParticle->GetParticleName() 
	//       << " xsec(barn)= " <<  elmCrossSection/barn 
	//       << "  E(MeV)= " << elmKinEnergy/MeV 
	//       << " Z= " << Z << "  A= " << iso->GetN() << "  j= " << j << G4endl;
      }
    }
  }
  //G4cout << "  E(MeV)= " << elmKinEnergy/MeV 
  //	 << "xsec(barn)= " <<  elmCrossSection/barn <<G4endl;
  return elmCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetIsoCrossSection(const G4DynamicParticle* part,
					    G4int Z, G4int A, 
					    const G4Isotope* iso,
					    const G4Element* elm,
					    const G4Material* mat, 
					    G4int idx)
{
  // this methods is called after the check that dataSetList[idx] 
  // depend on isotopes, so for this DataSet only isotopes are checked

  // isotope-wise cross section does exist
  if(dataSetList[idx]->IsIsoApplicable(part, Z, A, elm, mat) ) {
    return dataSetList[idx]->GetIsoCrossSection(part, Z, A, iso, elm, mat);

  } else {
    // seach for other dataSet
    for (G4int j = nDataSetList-1; j >= 0; --j) { 
      if (dataSetList[j]->IsElementApplicable(part, Z, mat)) {
	return dataSetList[j]->GetElementCrossSection(part, Z, mat);
      } else if (dataSetList[j]->IsIsoApplicable(part, Z, A, elm, mat)) {
	return dataSetList[j]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }
  }
  G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	 << " no isotope cross section found"
	 << G4endl;
  G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	 << " off Element " << elm->GetName()
         << "  in " << mat->GetName() 
	 << " Z= " << Z << " A= " << A
	 << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  throw G4HadronicException(__FILE__, __LINE__, 
                      " no applicable data set found for the isotope");
  return 0.0;
  //return dataSetList[idx]->ComputeCrossSection(part, elm, mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         G4int Z, G4int A,
					 const G4Isotope* iso,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  for (G4int i = nDataSetList-1; i >= 0; --i) {
    if (dataSetList[i]->IsIsoApplicable(part, Z, A, elm, mat) ) {
      return dataSetList[i]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
    }
  }
  G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	 << " no isotope cross section found"
	 << G4endl;
  G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	 << " off Element " << elm->GetName()
         << "  in " << mat->GetName() 
	 << " Z= " << Z << " A= " << A
	 << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  throw G4HadronicException(__FILE__, __LINE__, 
                      " no applicable data set found for the isotope");
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Element*
G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* part, 
                                     const G4Material* mat,
				     G4Nucleus& target)
{

  sampleZandACount++;
  G4int nElements = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  G4Element* anElement = (*theElementVector)[0];

  G4double cross = GetCrossSection(part, mat, true);

  // zero cross section case
  if(0.0 >= cross) { return anElement; }

  // select element from a compound 
  if(1 < nElements) {
    cross *= G4UniformRand();
    for(G4int i=0; i<nElements; ++i) {
      if(cross <= xsecelm[i]) {
	anElement = (*theElementVector)[i];
        break;
      }
    }
  }

  G4int Z = G4lrint(anElement->GetZ());
  G4Isotope* iso = 0;

  G4int i = nDataSetList-1; 
  if (dataSetList[i]->IsElementApplicable(part, Z, mat)) {

    //----------------------------------------------------------------
    // element-wise cross section
    // isotope cross section is not computed
    //----------------------------------------------------------------
    G4int nIso = anElement->GetNumberOfIsotopes();
    if (0 >= nIso) { 
      G4cout << " Element " << anElement->GetName() << " Z= " << Z 
	     << " has no isotopes " << G4endl; 
      throw G4HadronicException(__FILE__, __LINE__, 
                      " Isotope vector is not defined");
      return anElement;
    }
    // isotope abundances        
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    iso = (*isoVector)[0];

    // more than 1 isotope
    if(1 < nIso) { 
      iso = dataSetList[i]->SelectIsotope(anElement, part->GetKineticEnergy());
    }

  } else {

    //----------------------------------------------------------------
    // isotope-wise cross section
    // isotope cross section is computed
    //----------------------------------------------------------------
    G4int nIso = anElement->GetNumberOfIsotopes();
    cross = 0.0;

    if (0 >= nIso) { 
      G4cout << " Element " << anElement->GetName() << " Z= " << Z 
	     << " has no isotopes " << G4endl; 
      throw G4HadronicException(__FILE__, __LINE__, 
                      " Isotope vector is not defined");
      return anElement;
    }

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    iso = (*isoVector)[0];

    // more than 1 isotope
    if(1 < nIso) {
      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      if(G4int(xseciso.size()) < nIso) { xseciso.resize(nIso); }

      for (G4int j = 0; j<nIso; ++j) {
	G4double xsec = 0.0;
	if(abundVector[j] > 0.0) {
	  iso = (*isoVector)[j];
	  xsec = abundVector[j]*
	    GetIsoCrossSection(part, Z, iso->GetN(), iso, anElement, mat, i);
	}
	cross += xsec;
	xseciso[j] = cross;
      }
      cross *= G4UniformRand();
      for (G4int j = 0; j<nIso; ++j) {
	if(cross <= xseciso[j]) {
	  iso = (*isoVector)[j];
	  break;
	}
      }
    }
  }
  target.SetIsotope(iso);
  return anElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4CrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if (nDataSetList == 0) 
    {
      throw G4HadronicException(__FILE__, __LINE__, 
				"G4CrossSectionDataStore: no data sets registered");
      return;
    }
  for (G4int i=0; i<nDataSetList; ++i) {
    dataSetList[i]->BuildPhysicsTable(aParticleType);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4CrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  // Print out all cross section data sets used and the energies at
  // which they apply

  if (nDataSetList == 0) {
    G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: "
	   << " no data sets registered" << G4endl;
    return;
  }

  for (G4int i = nDataSetList-1; i >= 0; --i) {
    G4double e1 = dataSetList[i]->GetMinKinEnergy();
    G4double e2 = dataSetList[i]->GetMaxKinEnergy();
     G4cout 
      << "     Cr_sctns: " << std::setw(25) << dataSetList[i]->GetName() << ": "
      <<  G4BestUnit(e1, "Energy")
      << " ---> "
      <<  G4BestUnit(e2, "Energy") << "\n";
      if (dataSetList[i]->GetName() == "G4CrossSectionPairGG") {
        dataSetList[i]->DumpPhysicsTable(aParticleType);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition&,
                                       std::ofstream& outFile)
{
  // Write cross section data set info to html physics list
  // documentation page

  G4double ehi = 0;
  G4double elo = 0;
  for (G4int i = nDataSetList-1; i > 0; i--) {
    elo = dataSetList[i]->GetMinKinEnergy()/GeV;
    ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
    outFile << "      <li><b><a href=\"" << dataSetList[i]->GetName() << ".html\"> "
            << dataSetList[i]->GetName() << "</a> from "
            << elo << " GeV to " << ehi << " GeV </b></li>\n";
  }

  G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
  if (ehi < defaultHi) {
    outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() << ".html\"> "
            << dataSetList[0]->GetName() << "</a> from "
            << ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....



// Rob Fowler's simplify code

//  This is a curve simplification routine based on the Douglas-Peucker
//  algorithm.
//  Simplifying assumptions are that the input polyline is a piecewise
//  function with the x values monotonically increasing,  that the function
//  reaches an asymptote at the right (high energy) end.
//  Also, the correct error measure is the difference in y between the original
//  curve and the result.
//  In GEANT4 use, the assumption is that the calling program has identified
//  low- and high-energy cutoffs and that the vector passed in is restricted
//  to the region between the cutoffs.


// The raw_data vector comes in ordered left to right (small energy to large).
//  The simplified_data vector is initially empty.


int  simplify_function(G4double tolerance,
                       vector <Point> & raw_data,
                       vector<Point>  & simplified_data)
{

  int debug = 3;

  int  gap_left, gap_right;  // indices of the current region

  int stack_top, here;  
 
  G4double tolsq = tolerance*tolerance;  // Alternative to working with absolute values.

  vector <int> working_stack; 
  //A stack of the points to the right of the current interval that 
  // are known to be selected.

  here = 0;
  gap_right = raw_data.size() - 1;  // index of the last element.

  gap_left = 0;

  if (debug>2)  cout << "First and last elements  " << gap_left <<"  " <<gap_right <<endl;

  simplified_data.push_back(raw_data[0]); //copy first element over.

  if (debug > 2) cout << "first point ( 0  " 
		      <<simplified_data[0].e  <<",  "<<simplified_data[0].xs <<" )" <<endl;
  working_stack.push_back(gap_right); // 0th element on the stack.

  while ( !working_stack.empty() )
    {  G4double   a, slope, delta;
      G4double deltasq_max= tolsq;  
      int i_max;  

      gap_right = working_stack.back();  //get current TOS
      i_max = gap_right;      
      

      if ( (gap_left +1) < gap_right )  // At least three points in the range.
        {
          // co-efficients for the left to right affine line segment
          slope =  (raw_data[gap_right].xs - raw_data[gap_left].xs) /
            (raw_data[gap_right].e - raw_data[gap_left].e);
          a = raw_data[gap_left].xs - slope * raw_data[gap_left].e;
     
          for ( int i = gap_left +1; i <gap_right;  i++) {
            delta = raw_data[i].xs - a - slope * raw_data[i].e;
            if ( delta * delta > deltasq_max){  
              deltasq_max = delta * delta;
              i_max = i;
            }
          }
        }  else
        { if (debug >3) cout << "      Less than 3 point interval at [ "
			     << gap_left <<", " <<gap_right<< " ]" <<endl;
	}

      if(i_max < gap_right) {  //  Found a new point, push it on the stack
	working_stack.push_back(i_max);
	if (debug >3) cout << "         pushing point " << i_max <<endl;
	gap_right = i_max;
      }
      else { // didn't find a new point betweek gap_left and gap_right.
	simplified_data.push_back(raw_data[gap_right]);
	if (debug >1) cout << "inserting point (" 
			   <<gap_right  <<",  "<<raw_data[gap_right].e <<", "
			   << raw_data[gap_right].xs <<" )" <<endl;
	gap_left = gap_right;
	working_stack.pop_back(); 
	gap_right = working_stack.back();
	if (debug >3) cout << "      new gap_right " << gap_right <<endl;
      }
          
    }
  if (debug > 2) cout << "Simplified curve size "<< simplified_data.size() <<endl;
  return (simplified_data.size());
}          

// Rob Fowler's debias code   

//  This is a de-biasing routine applied after using a curve simplification 
//  routine based on the Douglas-Peucker
//  algorithm.
//  Simplifying assumptions are that the input polyline is a piecewise
//  function with the x values monotonically increasing, and
//  The right error measure is the difference in y between the original
//  curve and the result.



void RemoveBias(vector<Point> &  original, vector<Point> & simplified, 
                vector<Point> & result){
  

  G4int originalSize = original.size();
  G4int simplifiedSize = simplified.size();

  //Create index mapping array
  G4int xindex[simplifiedSize];
  G4int lastmatch = 0;
  G4int  j = 0;

  G4int debug = 3;

  if( debug >2)  cout  << "   original and simplified vector sizes  " << originalSize <<"  "
		       <<simplifiedSize <<endl;

  for (int k = 0; k <simplifiedSize; k++) { 
    for (int i = lastmatch; i < originalSize; i++) {
      if (original[i].e == simplified[k].e) {
        xindex[j++] = i;
        lastmatch = i;  
      }
    }
  }

  if (debug > 2 ) cout << "Matched  " << j << " values of the simplified vector" <<endl;
  // Use short names here.
  G4int m = simplifiedSize;

  G4double GArea [m-1];
  G4double GAreatotal = 0;

  //Area of original simplified curve                                                                                                                                
  for(int i = 0; i < m-1; i++){
    G4double GAreatemp = 0;

    for(j = xindex[i]; j< xindex[i+1]; j++){
      G4double trap = (original[j+1].xs + original[j].xs) * (original[j+1].e - original[j].e)/2.0;
      GAreatemp = GAreatemp + trap;
    }
    GArea[i] = GAreatemp;
    GAreatotal = GAreatotal + GAreatemp;
  }

  cout << "   Area under the original curve " << GAreatotal << endl;

  //aleph        Why is this not alpha?                                                                                                                                            


  G4double aleph [m-1];
  for(int i = 0; i< m-1; i++){
    aleph[i] = (simplified[i+1].e - simplified[i].e)/2.0;
  }
  
  //solve for f                                                                                                                                              
  G4double adjustedy [m];
  adjustedy[m-1] = simplified[m-1].xs;
  for(int i = 2; i < m+1; i++) {
    adjustedy[m-i] = (GArea[m-i]/aleph[m-i]) - adjustedy[m-i+1];
    if (adjustedy[m-i] <0.0) {
      adjustedy[m-i] = 0.0;
      cout << "   Fixing negative cross section at index " << (m-i) <<endl;
    }
  }
   
  //error and difference tracking                                                                                                                            
  G4double difference [m];
  G4double maxdiff = 0;
  G4double adjustedarea = 0;
  G4double simplifiedarea = 0;
  for(int i = 0; i < m-1; i++){  
    G4double trap;
    trap = (adjustedy[i+1]+adjustedy[i])*(simplified[i+1].e-simplified[i].e)/2.0;
    adjustedarea = adjustedarea+trap;
    trap = (simplified[i+1].xs+simplified[i].xs)*(simplified[i+1].e-simplified[i].e)/2.0;
    simplifiedarea = simplifiedarea + trap;
  }

  cout <<  "   Area:  Simplified curve = " <<simplifiedarea <<endl;
  cout << "    Area:  Debiased curve  = " << adjustedarea <<endl ;

  for(int i = 0; i <m; i++) {
    difference[i] = simplified[i].xs-adjustedy[i];
  }
  for(int i = 0; i <m; i++){
    if(fabs(difference[i]) > maxdiff) {
      maxdiff = fabs(difference[i]);
    }
  }
  //  what is the significance of the loops above ?
 
  for(int i = 0; i < simplifiedSize; i++){
    result.push_back( {simplified[i].e , adjustedy[i] } );
  }

}

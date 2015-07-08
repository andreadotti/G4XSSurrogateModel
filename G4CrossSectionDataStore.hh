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
// $Id: $
//
// -------------------------------------------------------------------
// File name:     G4CrossSectionDataStore
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source,
//                         use STL vector instead of C-array
//
// August 2011  Re-designed
//              by G. Folger, V. Ivantchenko, T. Koi and D.H. Wright

// Class Description
// This is the class to which cross section data sets may be registered. 
// An instance of it is contained in each hadronic process, allowing
// the use of the AddDataSet() method to tailor the cross sections to
// your application.
// Class Description - End

#ifndef G4CrossSectionDataStore_h
#define G4CrossSectionDataStore_h 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

//#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsFreeVector.hh"

#include <vector>

//PRUTH
#include <inttypes.h>
//#include <map>
#include <unordered_map>
//#include "psimpl.h"

#include "simplify.h"

struct fastPathEntry{
  //std::string particle;
  //std::string material;
  //G4String dataset;
  
  //G4double coefficient[9];
  //G4double knot[9];
  //G4double intercept;
  //G4int knot_cnt;

  G4ParticleDefinition *particle;
  G4Material *material;
  G4double min_cutoff;

  G4PhysicsVector *physicsVector;

  
  //stats for debug
  G4int count;
  G4double slowpath_sum; //sum of all slowpath xs
  G4double max_delta;
  G4double min_delta;
  G4double sum_delta;
  G4double sum_delta_square;
  

};

struct cycleCountEntry{
  //G4DynamicParticle* particle; 

  G4String particle;
  const G4Material* material;
  G4String dataset;

  // Comment out for speed 
  uint64_t initCyclesFastPath;
  uint64_t invocationCountSlowPath;
  uint64_t totalCyclesSlowPath;
  uint64_t invocationCountFastPath;
  uint64_t totalCyclesFastPath;
  uint64_t invocationCountTriedOneLineCache;
  uint64_t invocationCountOneLineCache;
 

  //optional fastPathEntry
  struct fastPathEntry* fastPath;

  //cache per element of material test 
  G4double energy;
  G4double crossSection;
  uint64_t cacheHitCount;
};



//

class G4Nucleus;
class G4DynamicParticle;
class G4ParticleDefinition;
class G4Isotope;
class G4Element;
class G4Material;
class G4NistManager;


//PRUTH
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <stdlib.h>
//#include "psimpl.h"

using namespace std;

namespace std {
  template <>
  struct hash<std::pair<G4ParticleDefinition*,const G4Material*>> {
  public:
    std::hash<uint64_t> hash_uint64_t;
    inline size_t operator()(std::pair<G4ParticleDefinition*,const G4Material*> x) const throw() {
      return hash_uint64_t(hash_uint64_t( ((uint64_t)(x.first)) ) +  hash_uint64_t(   ((uint64_t)(x.second))));
    }

  };


}


typedef std::pair<G4ParticleDefinition*,const G4Material*> G4CrossSectionDataStore_Key;
typedef std::unordered_map<G4CrossSectionDataStore_Key, struct cycleCountEntry*> G4CrossSectionDataStore_CycleCountMap;
//typedef std::unordered_map<G4CrossSectionDataStore_Key, struct fastPathEntry*>  G4CrossSectionDataStore_FastPathEntryMap;



class G4CrossSectionDataStore
{
public:
  //PRUTH
  static G4int objCnt;
 
  //static std::map<std::pair<G4ParticleDefinition*,const G4Material*>,struct cycleCountEntry*>* cycleCountMap;  
  //static std::unordered_map<std::pair<G4ParticleDefinition*,const G4Material*>,struct cycleCountEntry*>* cycleCountMap;
  G4CrossSectionDataStore_CycleCountMap cycleCountMap;
  ofstream resultsFile;
  G4int objID;

  //PRUTH vars for sampling and surragate model
  G4double queryMax;
  G4double sampleMin;
  G4double sampleMax;
  G4int sampleCount;
  G4double dpTol;

  //std::map<std::pair<G4ParticleDefinition*,const G4Material*>,struct fastPathEntry*>* fastPathMap;
  //std::unordered_map<std::pair<G4ParticleDefinition*,const G4Material*>,struct fastPathEntry*>* fastPathMap; 
  //G4CrossSectionDataStore_FastPathEntryMap fastPathMap;

  static G4long sampleZandACount;
  static G4long getCrossSectionCount;
  static G4long getCrossSectionCount_fastpath;
  static G4long getCrossSectionCount_slowpath;
  static G4long getCrossSectionCount_hitOneLineCache;

//

  const G4String processName;

  G4CrossSectionDataStore(const G4String& pname);

  ~G4CrossSectionDataStore();

  //PRUTH
  G4PhysicsVector* InitializeFastPathPhysicsVector(const G4DynamicParticle* part, const G4Material* mat, G4double cuttoff);
  G4double GetCrossSectonFastPath(struct fastPathEntry* fast_entry, const G4DynamicParticle* part);
  G4double SampleCrossSectionValue(const G4DynamicParticle* part, const G4Material* mat, G4double xval);
  //void RemoveBias(vector<G4double> data1, vector<G4double> data2, vector<G4double>& result);
  //void RemoveBias(vector<Point> &  original, vector<Point> & simplified, vector<Point> & result);
  void writeLinearXSSample(G4String label, G4double min, G4double max, G4double step, const G4DynamicParticle* part, const G4Material* mat);
  void writeXSVector(G4String label, vector <G4double> vect);

  // Cross section per unit volume is computed (inverse mean free path)
  G4double GetCrossSection(const G4DynamicParticle*, const G4Material*);


  G4double GetCrossSection(const G4DynamicParticle* part,
			   const G4Material* mat, G4bool requireSlowPath);

  // Cross section per element is computed
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, const G4Material*);

  // Cross section per isotope is computed
  G4double GetCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Isotope*,
			   const G4Element*, const G4Material*);

  // Sample Z and A of a target nucleus and upload into G4Nucleus
  G4Element* SampleZandA(const G4DynamicParticle*, const G4Material*,
			 G4Nucleus& target);

  // Initialisation before run
  void BuildPhysicsTable(const G4ParticleDefinition&);

  // Dump store to G4cout
  void DumpPhysicsTable(const G4ParticleDefinition&);

  // Dump store as html
  void DumpHtml(const G4ParticleDefinition&, std::ofstream&);

  inline void AddDataSet(G4VCrossSectionDataSet*);

  inline void SetVerboseLevel(G4int value);

private:

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
			      const G4Isotope*,
			      const G4Element*, const G4Material* aMaterial,
			      G4int index);

  G4CrossSectionDataStore & operator=(const G4CrossSectionDataStore &right);
  G4CrossSectionDataStore(const G4CrossSectionDataStore&);

  G4NistManager* nist;

  std::vector<G4VCrossSectionDataSet*> dataSetList;
  std::vector<G4double> xsecelm;
  std::vector<G4double> xseciso;

  const G4Material* currentMaterial;
  const G4ParticleDefinition* matParticle;
  G4double matKinEnergy;
  G4double matCrossSection;

  const G4Material* elmMaterial;
  const G4Element* currentElement;
  const G4ParticleDefinition* elmParticle;
  G4double elmKinEnergy;
  G4double elmCrossSection;

  G4int nDataSetList;
  G4int verboseLevel;
};

inline void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p)
{
  dataSetList.push_back(p);
  ++nDataSetList;
}

inline void G4CrossSectionDataStore::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

#endif

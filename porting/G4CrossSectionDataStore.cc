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
// $Id: G4CrossSectionDataStore.cc 90276 2015-05-22 10:30:24Z gunter $
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::G4CrossSectionDataStore() :
  nDataSetList(0), verboseLevel(0)
{
  nist = G4NistManager::Instance();
  currentMaterial = elmMaterial = 0;
  currentElement = 0;  //ALB 14-Aug-2012 Coverity fix.
  matParticle = elmParticle = 0;
  matKinEnergy = elmKinEnergy = matCrossSection = elmCrossSection = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::~G4CrossSectionDataStore()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Material* mat )
{
	//Reminder:
	//fastPathFlags contains control flags for the fast-path algorithm:
	//             .requiresSlowPath == true => Use slow path always
	// 	           .prevCalcUsedFastPath == true => Previous call to GetCrossSection used the fast-path
	//												it is used in the decision to assess if xsecelem is
	//												correctly set-up
	//			   .useFastPathIfAvailable == true => User requested the use of fast-path algorithm
	//			   .initializationPhase == true => If true we are in Geant4 Init phase before the event-loop
	//											   used in the decision to build the fast-path table
	//Starting from a fresh call, in general we want to use, if possible fast-path algorithm
	fastPathFlags.requiresSlowPath = false;
	//Check user-request, does he want fast-path?
	//if so, are we in initialization phase?
	if ( !fastPathFlags.useFastPathIfAvailable
		||	(fastPathFlags.useFastPathIfAvailable&&fastPathFlags.initializationPhase) ) {
		fastPathFlags.requiresSlowPath=true;
	}

	//Logging for performance calculations and counter, active only in FPDEBUG mode
	counters.MethodCalled();
	//Measure number of cycles
	G4FastPathHadronicCrossSection::logStartCountCycles(timing);

	//This is the cache entry of the fast-path cross-section parametrization
	G4FastPathHadronicCrossSection::cycleCountEntry* entry = nullptr;
	//Did user requrest fast-path in first place and are we not in the initialization phase
	if ( fastPathFlags.useFastPathIfAvailable && !fastPathFlags.initializationPhase ) {
		//Important: if it is in initialization phase we should NOT use fast path: we are going to build it
		//G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Key searchkey = {part->GetParticleDefinition(),mat};
		entry = fastPathCache[{part->GetParticleDefinition(),mat}];
	}
	assert( fastPathFlags.initializationPhase && entry == nullptr );
	assert( fastPathFlags.useFastPathIfAvailable && entry == nullptr );

  //Super fast check: are we calling again this method for exactly the same interaction?
  if(mat == currentMaterial && part->GetDefinition() == matParticle
     && part->GetKineticEnergy() == matKinEnergy) 
    {
	  G4FastPathHadronicCrossSection::logInvocationTriedOneLine(entry);
	  //Check that the previous time we called this method we used the slow
	  //path: we need the data-member xsecelem to be the one for the current
	  //interaction, this is ensured only if: we will do the slow path right now or we
	  //did it exactly for the same conditions the last call.
	  if ( !fastPathFlags.prevCalcUsedFastPath && !fastPathFlags.requiresSlowPath ) {
		  counters.HitOneLine();
		  G4FastPathHadronicCrossSection::logInvocationOneLine(entry);
		  //Good everything is setup correctly, exit!
		  return matCrossSection;
	  } else {
		  fastPathFlags.requiresSlowPath = true;
	  }
    }
  
  //Ok, now check if we have cached for this {particle,material,energy} the cross-section
  //in this case we did, let's return immediately, if we are not forced to take the slow path
  //(e.g. as before if the xsecelem is not up-to-date we need to take the slow-path).
  //Note that this is not equivalent to the previous ultra-fast check: we now have a map that
  //we are using.
  //TODO: I think the previous if is contained in this case, so we could merge the
  //      two ifs. The previous one is the algorithm in G4Ver<10.2
  if ( entry != nullptr && entry->energy == part->GetKineticEnergy() ) {
	  G4FastPathHadronicCrossSection::logHit(entry);
	  if ( !fastPathFlags.requiresSlowPath ) {
		  return entry->crossSection;
	  }
  }

  currentMaterial = mat;
  matParticle = part->GetDefinition();
  matKinEnergy = part->GetKineticEnergy();
  matCrossSection = 0;

  //Now check if the cache entry has a fast-path cross-section calculation available
  G4FastPathHadronicCrossSection::fastPathEntry* fast_entry = nullptr;
  if ( entry != nullptr && !fastPathFlags.requiresSlowPath ) {
	  fast_entry = entry->fastPath;
  }

  //Each fast-path cross-section has a minimum value of validity, if energy is below
  //that skip
  if ( fast_entry != nullptr && part->GetKineticEnergy() < fast_entry->min_cutoff )
  {
	  assert(fastPathFlags.requiresSlowPath==false);
	  fastPathFlags.requiresSlowPath = true;
  }

  //Ready to use the fast-path calculation
  if ( !fastPathFlags.requiresSlowPath && fast_entry != nullptr ) {
	  counters.FastPath();
	  //Retrieve cross-section from fast-path cache
	  matCrossSection = fast_entry->GetCrossSection(part->GetKineticEnergy());
	  fastPathFlags.prevCalcUsedFastPath=true;
  } else {
	  counters.SlowPath();
	  //Remember that we are now doing the full calculation: xsecelem will
	  //be made valid
	  fastPathFlags.prevCalcUsedFastPath=false;

	  G4int nElements = mat->GetNumberOfElements();
	  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

	  if(G4int(xsecelm.size()) < nElements) { xsecelm.resize(nElements); }

	  for(G4int i=0; i<nElements; ++i) {
		  matCrossSection += nAtomsPerVolume[i] *
				  GetCrossSection(part, (*mat->GetElementVector())[i], mat);
		  xsecelm[i] = matCrossSection;
	  }
  }
  //Stop measurement of cpu cycles
  G4FastPathHadronicCrossSection::logStopCountCycles(timing);

  //We are in initialization phase, we want to use the fast-path
  if ( fastPathFlags.useFastPathIfAvailable && fastPathFlags.initializationPhase ) {
	  //Check if this particular {particle,material} combination has never been seen before
	  G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Key searchkey = {part->GetParticleDefinition(),mat};
	  entry = fastPathCache[searchkey];
	  if (entry == nullptr){
		  //add entry
		  entry = new G4FastPathHadronicCrossSection::cycleCountEntry();
		  entry->particle = matParticle->GetParticleName(); //TODO: needed?
		  entry->material = mat; //TODO: needed?
		  entry->dataset = dataSetList[nDataSetList-1]->GetName();//TODO: needed?
		  fastPathCache[searchkey] = entry;
	  }
  }

  if ( entry != nullptr ) {
	  entry->energy = part->GetKineticEnergy();
	  entry->crossSection = matCrossSection;
  }
  //Some logging of timing
  G4FastPathHadronicCrossSection::logTiming(entry,fast_entry,timing);
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4CrossSectionDataStore::DumpFastPath(const G4ParticleDefinition* pd, const G4Material* mat,std::ostream& os)
{
	const G4FastPathHadronicCrossSection::cycleCountEntry* entry = fastPathCache[{pd,mat}];
	if ( entry != nullptr ) {
		if ( entry->fastPath != nullptr ) {
			os<<*entry->fastPath;
		} else {
			os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
			os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} found, but no fast path defined";
		}
	} else {
		os<<"#Cache entry for {"<<(pd!=nullptr?pd->GetParticleName():"UNDEFINED")<<",";
		os<<(mat!=nullptr?mat->GetName():"UNDEFINED")<<"} not found.";
	}
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
  G4int nElements = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  G4Element* anElement = (*theElementVector)[0];

  G4double cross = GetCrossSection(part, mat);

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
#include <typeinfo>
void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition& /* pD */,
                                       std::ofstream& outFile) const
{
  // Write cross section data set info to html physics list
  // documentation page

  G4double ehi = 0;
  G4double elo = 0;
  G4String physListName(getenv("G4PhysListName"));
  for (G4int i = nDataSetList-1; i > 0; i--) {
    elo = dataSetList[i]->GetMinKinEnergy()/GeV;
    ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
    outFile << "      <li><b><a href=\"" << physListName << "_"
	         << dataSetList[i]->GetName() << ".html\"> "
            << dataSetList[i]->GetName() << "</a> from "
            << elo << " GeV to " << ehi << " GeV </b></li>\n";
	 //G4cerr << i << ": XS for " << pD.GetParticleName() << " : " << dataSetList[i]->GetName() 
	 //       << " typeid : " << typeid(dataSetList[i]).name()<< G4endl;			
	 PrintCrossSectionHtml(dataSetList[i]);			
  }

  G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
  if (ehi < defaultHi) {
    outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() << ".html\"> "
            << dataSetList[0]->GetName() << "</a> from "
            << ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
	 PrintCrossSectionHtml(dataSetList[0]);			
  }
}

void G4CrossSectionDataStore::PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const
{
   G4String dirName(getenv("G4PhysListDocDir"));
	G4String physListName(getenv("G4PhysListName"));

	G4String pathName = dirName + "/" + physListName + "_" + HtmlFileName(cs->GetName());
	std::ofstream outCS;
	outCS.open(pathName);
	outCS << "<html>\n";
	outCS << "<head>\n";
	outCS << "<title>Description of " << cs->GetName() 
		 << "</title>\n";
	outCS << "</head>\n";
	outCS << "<body>\n";

	cs->CrossSectionDescription(outCS);

	outCS << "</body>\n";
	outCS << "</html>\n";

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//private 
G4String G4CrossSectionDataStore::HtmlFileName(const G4String & in) const
{
   G4String str(in);
    // replace blanks by _  C++11 version:
#ifdef G4USE_STD11
	std::transform(str.begin(), str.end(), str.begin(), [](char ch) {
     return ch == ' ' ? '_' : ch;
   });
#else	
	  // and now in ancient language
	   for(std::string::iterator it = str.begin(); it != str.end(); ++it) {
        if(*it == ' ') *it = '_';
      }
#endif
   str=str + ".html";		
   return str;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

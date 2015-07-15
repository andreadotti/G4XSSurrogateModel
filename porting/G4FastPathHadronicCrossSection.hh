#ifndef G4FastPathHadronicCrossSection_hh
#define G4FastPathHadronicCrossSection_hh

#include "G4PhysicsFreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include <functional>
#include <utility>
#include <unordered_map>

#define FPDEBUG

namespace G4FastPathHadronicCrossSection {
	//This data type contains the simplified representation of the
	//cross-section, by default it is a G4PhysicsVector type
	using XSParam=G4PhysicsFreeVector;
	//The key used to search in the cache.
	using G4CrossSectionDataStore_Key=std::pair<G4ParticleDefinition*,const G4Material*>;
	//This represents the fast XS implementation.
	struct fastPathEntry{
		fastPathEntry();
		~fastPathEntry();
		G4ParticleDefinition *particle;
		G4Material *material;
		G4double min_cutoff;

		XSParam *physicsVector;

#       ifdef FPDEBUG
		//stats for debug
		G4int count;
		G4double slowpath_sum; //sum of all slowpath xs
		G4double max_delta;
		G4double min_delta;
		G4double sum_delta;
		G4double sum_delta_square;
#		endif
	};

	//A cache entry.
	struct cycleCountEntry{
		cycleCountEntry();
		~cycleCountEntry();
		G4String particle;
		const G4Material* material;
		G4String dataset;

		//optional fastPathEntry
		fastPathEntry* fastPath;

		//cache per element of material test
		G4double energy;
		G4double crossSection;
#	  ifdef FPDEBUG
		uint64_t cacheHitCount;//
		uint64_t initCyclesFastPath;
		uint64_t invocationCountSlowPath;
		uint64_t totalCyclesSlowPath;
		uint64_t invocationCountFastPath;
		uint64_t totalCyclesFastPath;
		uint64_t invocationCountTriedOneLineCache;//
		uint64_t invocationCountOneLineCache;//
		unsigned long long rdtsc_start;
		unsigned long long rdtsc_stop;
#	  endif
	};

	struct getCrossSectionCount {
		getCrossSectionCount();
		inline void MethodCalled();
		inline void HitOneLine();
		inline void FastPath();
		inline void SlowPath();
#ifdef FPDEBUG
		uint64_t methodCalled;
		uint64_t hitOneLineCache;
		uint64_t fastPath;
		uint64_t slowPath;
#endif
	};

	//Hashing the key
	struct G4CrossSectionDataStore_Key_Hash {
		std::hash<uint64_t> hash_uint64_t;
		inline size_t operator()(std::pair<G4ParticleDefinition*,const G4Material*> x) const throw() {
			return hash_uint64_t(hash_uint64_t( ((uint64_t)(x.first)) ) +  hash_uint64_t(   ((uint64_t)(x.second))));
		}
	};
	//Equality for two key elements
	struct G4CrossSectionDataStore_Key_EqualTo {
		inline bool operator() (	const G4CrossSectionDataStore_Key& lhs, const G4CrossSectionDataStore_Key& rhs ) {
			//TODO: Verify this: particles are singletons, materials use operator==
			return (lhs.first==rhs.first)&&(*lhs.second == *rhs.second);
		}
	};

	//The cache itself
	using G4CrossSectionDataStore_Cache=std::unordered_map<G4CrossSectionDataStore_Key,cycleCountEntry,G4CrossSectionDataStore_Key_Hash,G4CrossSectionDataStore_Key_EqualTo>;

	//Configure the caching mechanism
	struct controlFlag {
		G4bool prevCalcUsedFastPath;
	};

	//Logging functionalities, disabled if not in FPDEBUG mode
	static inline void logInvocationTriedOneLine( cycleCountEntry* );
	static inline void logInvocationOneLine( cycleCountEntry* );
	static inline void logHit(cycleCountEntry*);
	static inline void logInvocationCountFastPath( cycleCountEntry* );
	static inline void logInvocationCountSlowPAth( cycleCountEntry* );
	static inline void logStartCountCycles( cycleCountEntry* );
	static inline void logStopCountCycles( cycleCountEntry* );
	static inline void logInitCyclesFastPath( cycleCountEntry* );
	static inline void logTotalCyclesFastPath( cycleCountEntry* );
	static inline void logTotalCyclesSlowPath( cycleCountEntry* );
}

#endif //G4FastPathHadronicCrossSection_hh

#ifndef G4FastPathHadronicCrossSection_hh
#define G4FastPathHadronicCrossSection_hh

#include "G4PhysicsFreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include <functional>
#include <utility>
#include <unordered_map>

#define FPDEBUG
//TODO: Move all logging and debug functionality to separate header
namespace G4FastPathHadronicCrossSection {
	//This data type contains the simplified representation of the
	//cross-section, by default it is a G4PhysicsVector type
	using XSParam=G4PhysicsFreeVector;
	//The key used to search in the cache.
	using G4CrossSectionDataStore_Key=std::pair<const G4ParticleDefinition*,const G4Material*>;
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
#	  endif
	};

	struct timing {
		unsigned long long rdtsc_start;
		unsigned long long rdtsc_stop;
	};

	struct getCrossSectionCount {
		getCrossSectionCount();
		inline void MethodCalled();
		inline void HitOneLine();
		inline void FastPath();
		inline void SlowPath();
		inline void SampleZandA();
#ifdef FPDEBUG
		uint64_t methodCalled;
		uint64_t hitOneLineCache;
		uint64_t fastPath;
		uint64_t slowPath;
		uint64_t sampleZandA;
#endif
	};

	//Hashing the key
	struct G4CrossSectionDataStore_Key_Hash {
		std::hash<uint64_t> hash_uint64_t;
		inline size_t operator()(const G4CrossSectionDataStore_Key& x) const throw() {
			return hash_uint64_t(hash_uint64_t( ((uint64_t)(x.first)) ) +  hash_uint64_t(   ((uint64_t)(x.second))));
		}
	};
	//Equality for two key elements
	struct G4CrossSectionDataStore_Key_EqualTo {
		inline bool operator()(const G4CrossSectionDataStore_Key& lhs, const G4CrossSectionDataStore_Key& rhs ) const {
			//TODO: Verify this: particles are singletons, materials use operator==
			return (lhs.first==rhs.first)&&(*lhs.second == *rhs.second);
		}
	};

	//The cache itself
	using G4CrossSectionDataStore_Cache=std::unordered_map<G4CrossSectionDataStore_Key,cycleCountEntry*,G4CrossSectionDataStore_Key_Hash,G4CrossSectionDataStore_Key_EqualTo>;

	//Configure the caching mechanism
	struct controlFlag {
		G4bool prevCalcUsedFastPath;
		G4bool useFastPathIfAvailable;
		G4bool initializationPhase;
		controlFlag() : prevCalcUsedFastPath(false),useFastPathIfAvailable(false),initializationPhase(false) {}
	};

	//Logging functionalities, disabled if not in FPDEBUG mode
	static inline void logInvocationTriedOneLine( cycleCountEntry* );
	static inline void logInvocationOneLine( cycleCountEntry* );
	static inline void logHit(cycleCountEntry*);
	static inline void logInvocationCountFastPath( cycleCountEntry* );
	static inline void logInvocationCountSlowPAth( cycleCountEntry* );
	void logStartCountCycles( timing& );
	void logStopCountCycles( timing& );
	static inline void logInitCyclesFastPath( cycleCountEntry* , timing& );
	static inline void logTotalCyclesFastPath( cycleCountEntry* , timing& );
	static inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& );
}


namespace G4FastPathHadronicCrossSection {

#ifdef FPDEBUG
inline void logInvocationTriedOneLine(cycleCountEntry* cl ) {
	if ( cl != nullptr ) ++(cl->invocationCountTriedOneLineCache);
}
inline void logInvocationOneLine( cycleCountEntry* cl ) {
	if ( cl != nullptr ) ++(cl->invocationCountOneLineCache);
}
inline void logHit(cycleCountEntry* cl) {
	if ( cl != nullptr ) ++(cl->cacheHitCount);
}
inline void logInvocationCountFastPath( cycleCountEntry* cl )
{
	if ( cl != nullptr ) ++(cl->invocationCountFastPath);
}
inline void logInvocationCountSlowPAth( cycleCountEntry* cl)
{
	if ( cl != nullptr ) ++(cl->invocationCountSlowPath);
}

inline void logInitCyclesFastPath(cycleCountEntry* cl,timing& tm)
{
	if ( cl != nullptr ) cl->initCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
}
inline void logTotalCyclesFastPath( cycleCountEntry* cl,timing& tm)
{
	if ( cl!=nullptr ) cl->totalCyclesFastPath = tm.rdtsc_stop - tm.rdtsc_start;
}
inline void logTotalCyclesSlowPath( cycleCountEntry* cl,timing& tm)
{
	if ( cl!=nullptr ) cl->totalCyclesSlowPath = tm.rdtsc_stop - tm.rdtsc_start;
}
#else
inline void logInvocationTriedOneLine(cycleCountEntry*){}
inline void logInvocationOneLine( cycleCountEntry*){}
inline void logHit(cycleCountEntry*){}
inline void logInvocationCountFastPath( cycleCountEntry*){}
inline void logInvocationCountSlowPAth( cycleCountEntry*){}
inline void logInitCyclesFastPath( cycleCountEntry* , timing& ){}
inline void logTotalCyclesFastPath( cycleCountEntry* , timing& ){}
inline void logTotalCyclesSlowPath( cycleCountEntry* , timing& ){}
#endif

inline void getCrossSectionCount::MethodCalled() {
#ifdef FPDEBUG
	++methodCalled;
#endif
}

inline void getCrossSectionCount::HitOneLine() {
#ifdef FPDEBUG
	++hitOneLineCache;
#endif
}

inline void getCrossSectionCount::FastPath() {
#ifdef FPDEBUG
	++fastPath;
#endif
}

inline void getCrossSectionCount::SlowPath() {
#ifdef FPDEBUG
	++slowPath;
#endif
}

inline void getCrossSectionCount::SampleZandA() {
#ifdef FPDEBUG
	++sampleZandA;
#endif
}
}

#endif //G4FastPathHadronicCrossSection_hh

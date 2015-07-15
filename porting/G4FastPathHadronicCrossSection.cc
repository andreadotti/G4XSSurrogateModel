#include "G4FastPathHadronicCrossSection.hh"
#include "G4ios.hh"

#ifdef FPDEBUG
#define DBG( msg ) G4cout<< msg <<G4endl;
#define DUMP() G4cout<< <<G4endl;
#else
#define DBG(msg)
#define DUMP()
#endif

using namespace G4FastPathHadronicCrossSection;

fastPathEntry::fastPathEntry() :
	particle(nullptr),material(nullptr),min_cutoff(-1.),
	physicsVector(nullptr)
{
	DBG("Initializing a fastPathEntry");
#ifdef FPDEBUG
	count = 0;
	slowpath_sum=0.;
	max_delta=0.;
	min_delta=0.;
	sum_delta=0.;
	sum_delta_square=0.;
#endif
}

fastPathEntry::~fastPathEntry()
{
	DBG("Deleting fastPathEntry");
	DBG("Dumping status for: "<<(particle?particle->GetParticleName():"PART_NONE")<<" "\
		  <<(material?material->GetName():"MAT_NONE")<<" min_cutoff:"<<min_cutoff<<" "\
		  <<" count:"<<count<<" slowpath_sum:"<<slowpath_sum<<" max_delta:"<<max_delta\
		  <<" min_delta"<<min_delta<<" sum_delta"<<sum_delta<<" sum_delta_squared:"<<sum_delta_square);
	delete physicsVector;
}

cycleCountEntry::cycleCountEntry() :
		particle(""),material(nullptr),dataset(""),fastPath(nullptr),
		energy(-1.),crossSection(-1.)
{
	DBG("Initializing cache entry");
#ifdef FPDEBUG
	cacheHitCount = 0;
	initCyclesFastPath=0;
	invocationCountSlowPath=0;
	totalCyclesSlowPath=0;
	invocationCountFastPath=0;
	totalCyclesFastPath=0;
	invocationCountTriedOneLineCache=0;
	invocationCountOneLineCache=0;
#endif
}

cycleCountEntry::~cycleCountEntry()
{
	DBG("Deleting cache entry");
	DBG(particle<<" "<<material<<" ("<<(material?material->GetName():"MAT_NONE")<<") "<<dataset<<" "\
			<<"fast path pointer:"<<fastPath<<" stored:"<<energy<<" "<<crossSection<<" "\
			<<cacheHitCount<<" "<<initCyclesFastPath<<" "<<invocationCountSlowPath<<" "\
			<<totalCyclesSlowPath<<" "<<invocationCountFastPath<<" "<<totalCyclesFastPath<<" "\
			<<invocationCountTriedOneLineCache<<" "<<invocationCountOneLineCache);
}

#ifdef FPDEBUG
namespace {
	static inline unsigned long long rdtsc() {
		unsigned hi=0,lo=0;
#if defined(__GNUC__) &&( defined(__i386__)|| defined(__x86_64__) )
		__asm__ __volatile__ ("rdtsc":"=a"(lo),"=d"(hi));
#endif
		return ((unsigned long long)lo) | ((unsigned long long)hi<<32 );
	}
}
void G4FastPathHadronicCrossSection::logStartCountCycles(timing& tm)
{
	tm.rdtsc_start=rdtsc();
}
void G4FastPathHadronicCrossSection::logStopCountCycles(timing& tm)
{
	tm.rdtsc_stop=rdtsc();
}
#else
void logStartCountCycles(timing&) {}
void logStopCountCycles(timing&) {}
#endif
getCrossSectionCount::getCrossSectionCount() {
#ifdef FPDEBUG
	methodCalled = 0;
	hitOneLineCache=0;
	fastPath=0;
	slowPath=0;
	sampleZandA = 0;
#endif
}



#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_xsect
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_xsect
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 86150 2014-11-07 12:27:27Z vnivanch $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_xsect
    HEADERS
	G4BGGNucleonElasticXS.hh
	G4BGGNucleonInelasticXS.hh
	G4BGGPionElasticXS.hh
	G4BGGPionInelasticXS.hh
	G4ChipsAntiBaryonElasticXS.hh
	G4ChipsAntiBaryonInelasticXS.hh
	G4ChipsComponentXS.hh
	G4ChipsHyperonElasticXS.hh
	G4ChipsHyperonInelasticXS.hh
	G4ChipsKaonMinusElasticXS.hh
	G4ChipsKaonMinusInelasticXS.hh
	G4ChipsKaonPlusElasticXS.hh
	G4ChipsKaonPlusInelasticXS.hh
	G4ChipsKaonZeroElasticXS.hh
	G4ChipsKaonZeroInelasticXS.hh
	G4ChipsNeutronElasticXS.hh
	G4ChipsNeutronInelasticXS.hh
	G4ChipsPionMinusElasticXS.hh
	G4ChipsPionMinusInelasticXS.hh
	G4ChipsPionPlusElasticXS.hh
	G4ChipsPionPlusInelasticXS.hh
	G4ChipsProtonElasticXS.hh
	G4ChipsProtonInelasticXS.hh
	G4ComponentAntiNuclNuclearXS.hh
	G4ComponentBarNucleonNucleusXsc.hh
	G4ComponentGGHadronNucleusXsc.hh
	G4ComponentGGNuclNuclXsc.hh
	G4ComponentSAIDTotalXS.hh
	G4CrossSectionDataSetRegistry.hh
	G4CrossSectionDataStore.hh
	G4CrossSectionElastic.hh
	G4CrossSectionFactory.hh
	G4CrossSectionInelastic.hh
	G4CrossSectionPairGG.hh
	G4ElectroNuclearCrossSection.hh
        G4DiffElasticRatio.hh
	G4EMDissociationCrossSection.hh
	G4EMDissociationSpectrum.hh
	G4GeneralSpaceNNCrossSection.hh
	G4GGNuclNuclCrossSection.hh
	G4GlauberGribovCrossSection.hh
	G4HadronCaptureDataSet.hh
	G4HadronCrossSections.hh
	G4HadronElasticDataSet.hh
	G4HadronFissionDataSet.hh
	G4HadronInelasticDataSet.hh
	G4HadronNucleonXsc.hh
	G4IonProtonCrossSection.hh
	G4IonsKoxCrossSection.hh
	G4IonsShenCrossSection.hh
	G4IonsSihverCrossSection.hh
	G4KokoulinMuonNuclearXS.hh
	G4NeutronCaptureXS.hh
	G4NeutronElasticXS.hh
	G4NeutronInelasticCrossSection.hh
	G4NeutronInelasticXS.hh
	G4NucleonNuclearCrossSection.hh
	G4PhotoNuclearCrossSection.hh
	G4PiData.hh
	G4PiNuclearCrossSection.hh
	G4ProjectileFragmentCrossSection.hh
	G4ProtonInelasticCrossSection.hh
	G4TripathiCrossSection.hh
	G4TripathiLightCrossSection.hh
	G4UPiNuclearCrossSection.hh
	G4VComponentCrossSection.hh
	G4VCrossSectionDataSet.hh
	G4VCrossSectionRatio.hh
	G4CrossSectionFactoryRegistry.hh
	G4FastPathHadronicCrossSection.hh
    SOURCES
	G4BGGNucleonElasticXS.cc
	G4BGGNucleonInelasticXS.cc
	G4BGGPionElasticXS.cc
	G4BGGPionInelasticXS.cc
	G4ChipsAntiBaryonElasticXS.cc
	G4ChipsAntiBaryonInelasticXS.cc
	G4ChipsComponentXS.cc
	G4ChipsHyperonElasticXS.cc
	G4ChipsHyperonInelasticXS.cc
	G4ChipsKaonMinusElasticXS.cc
	G4ChipsKaonMinusInelasticXS.cc
	G4ChipsKaonPlusElasticXS.cc
	G4ChipsKaonPlusInelasticXS.cc
	G4ChipsKaonZeroElasticXS.cc
	G4ChipsKaonZeroInelasticXS.cc
	G4ChipsNeutronElasticXS.cc
	G4ChipsNeutronInelasticXS.cc
	G4ChipsPionMinusElasticXS.cc
	G4ChipsPionMinusInelasticXS.cc
	G4ChipsPionPlusElasticXS.cc
	G4ChipsPionPlusInelasticXS.cc
	G4ChipsProtonElasticXS.cc
	G4ChipsProtonInelasticXS.cc
	G4ComponentAntiNuclNuclearXS.cc
	G4ComponentBarNucleonNucleusXsc.cc
	G4ComponentGGHadronNucleusXsc.cc
	G4ComponentGGNuclNuclXsc.cc
	G4ComponentSAIDTotalXS.cc
	G4CrossSectionDataSetRegistry.cc
	G4CrossSectionDataStore.cc
	G4CrossSectionElastic.cc
	G4CrossSectionInelastic.cc
	G4CrossSectionPairGG.cc
        G4DiffElasticRatio.cc
	G4ElectroNuclearCrossSection.cc
	G4EMDissociationCrossSection.cc
	G4EMDissociationSpectrum.cc
	G4GeneralSpaceNNCrossSection.cc
	G4GGNuclNuclCrossSection.cc
	G4GlauberGribovCrossSection.cc
	G4HadronCaptureDataSet.cc
	G4HadronCrossSections.cc
	G4HadronElasticDataSet.cc
	G4HadronFissionDataSet.cc
	G4HadronInelasticDataSet.cc
	G4HadronNucleonXsc.cc
	G4IonProtonCrossSection.cc
	G4IonsKoxCrossSection.cc
	G4IonsShenCrossSection.cc
	G4IonsSihverCrossSection.cc
	G4KokoulinMuonNuclearXS.cc
	G4NeutronCaptureXS.cc
	G4NeutronElasticXS.cc
	G4NeutronInelasticCrossSection.cc
	G4NeutronInelasticXS.cc
	G4NucleonNuclearCrossSection.cc
	G4PhotoNuclearCrossSection.cc
	G4PiData.cc
	G4PiNuclearCrossSection.cc
	G4ProtonInelasticCrossSection.cc
	G4TripathiCrossSection.cc
	G4TripathiLightCrossSection.cc
	G4UPiNuclearCrossSection.cc
	G4VComponentCrossSection.cc
	G4VCrossSectionDataSet.cc
	G4VCrossSectionRatio.cc
	G4CrossSectionFactoryRegistry.cc
	G4FastPathHadronicCrossSection.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_util
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here


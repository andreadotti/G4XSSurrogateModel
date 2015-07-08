#!/bin/bash

module add gcc/4.8.1
#module add gcc/4.6.3

#export PATH=$PATH:/hpc_shared/geant4/blanco/cmake-2.8.7-INSTALL/bin

BUILD_NAME='geant4.10.00_SurrogateModel_v2'

#-----------------------------------------------------------------------
# Default options
#-----------------------------------------------------------------------

unset G4P_USE_CLHEP ; G4P_USE_CLHEP=0
unset G4P_CXX       ; G4P_CXX=/usr/bin/c++
#G4P_CXX='/hpc_shared/apps/RHEL-6/x86_64/gcc/gcc-4.5.1/bin/g++'
#G4P_CXX='/hpc_shared/apps/RHEL-6/x86_64/gcc/gcc-4.6.3/bin/g++'
G4P_CXX='/hpc_shared/apps/RHEL-6/x86_64/gcc/gcc-4.8.1/bin/g++'
export LD_LIBRARY_PATH=/hpc_shared/apps/RHEL-6/x86_64/gcc/gcc-4.8.1/lib64 

#G4P_CXX='/hpc_shared/home/pruth/GEANT-4/gcc/install/hpc_shared/home/pruth/GEANT-4/gcc/install/usr/bin/g++'
#export LD_LIBRARY_PATH=/hpc_shared/home/pruth/GEANT-4/gcc/install/hpc_shared/home/pruth/GEANT-4/gcc/install/usr/lib64

#-----------------------------------------------------------------------
# Create the directory structure
#-----------------------------------------------------------------------
umask 0002

WORK_DIR="/hpc_shared/home/pruth/GEANT-4"
BUILD_DIR="${WORK_DIR}/${BUILD_NAME}-build"
INSTALL_DIR="${WORK_DIR}/${BUILD_NAME}-install"
SOURCE_DIR="${WORK_DIR}/${BUILD_NAME}"

mkdir -p ${INSTALL_DIR}
mkdir -p ${BUILD_DIR}

cd ${BUILD_DIR}

#-----------------------------------------------------------------------
# configure with cmake
#-----------------------------------------------------------------------
XERCESC_DIR=/hpc_shared/geant4/blanco/xerces-c-3.1.1-INSTALL/
#XERCESC_DIR=/usr/
export XERCESC_DIR

#-g -DG4clhep_EXPORTS -DG4_STORE_TRAJECTORY -DG4VERBOSE -DCLHEP_EXPORT -W -Wall -ansi -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -pipe -DG4MULTITHREADED -ftls-model=initial-exec -O2  -fPIC

#-DCMAKE_CXX_FLAGS="-std=c++11 -O3 -g -fpermissive -fno-omit-frame-pointer -ftls-model=initial-exec  -fPIC " \
#      -DCMAKE_CXX_FLAGS_RELEASE="-std=c++11  -O3 -DNDEBUG -fpermissive -ftls-model=initial-exec  -fPIC " \
#      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-std=c++11  -O3 -g -fpermissive -fno-omit-frame-pointer -DG4clhep_EXPORTS -DG4_STORE_TRAJECTORY -DG4VERBOSE -DCLHEP_EXPORT -W -Wall -ansi -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -pipe -ftls-model=initial-exec  -fPIC   " \
#      -DCMAKE_BUILD_TYPE='RELWITHDEBINFO' \

#-DXERCESC_LIBRARY=${XERCESC_DIR}/lib/libxerces-c-3.1.so \
#-DXERCESC_INCLUDE_DIR=${XERCESC_DIR}/include \

#-DXERCESC_LIBRARY=${XERCESC_DIR}/lib64/libxerces-c-3.0.so \ 


cmake -DCMAKE_CXX_COMPILER="${G4P_CXX}" \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_INSTALL_LIBDIR=lib \
-DCMAKE_CXX_FLAGS="-std=c++0x -O3 -g -fpermissive -fno-omit-frame-pointer -ftls-model=initial-exec  -fPIC " \
      -DCMAKE_CXX_FLAGS_RELEASE="-std=c++0x  -O3 -DNDEBUG -fpermissive -ftls-model=initial-exec  -fPIC " \
      -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-std=c++0x  -O3 -g -fpermissive -fno-omit-frame-pointer -DG4clhep_EXPORTS -DG4_STORE_TRAJECTORY -DG4VERBOSE -DCLHEP_EXPORT -W -Wall -ansi -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -pipe -ftls-model=initial-exec  -fPIC   " \
      -DCMAKE_BUILD_TYPE='RELWITHDEBINFO' \
     -DGEANT4_BUILD_CXXSTD='c++0x' \
      -DGEANT4_BUILD_MULTITHREADED=ON \
      -DGEANT4_USE_SYSTEM_CLHEP=0 \
      -DGEANT4_INSTALL_DATA=0 \
      -DGEANT4_USE_GDML=ON \
      -DXERCESC_LIBRARY=${XERCESC_DIR}/lib/libxerces-c-3.1.so \
      -DXERCESC_INCLUDE_DIR=${XERCESC_DIR}/include \
      ${SOURCE_DIR}/source ${SOURCE_DIR}
#-----------------------------------------------------------------------
# build and install
#-----------------------------------------------------------------------
make -j16
#make

#make install

cmake_minimum_required(VERSION 3.10)
project(myfitter)

# Usa C++17 per garantire compatibilità con ROOT e GenFit
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova ROOT con i componenti specifici
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist MathCore Geom)
include_directories(${ROOT_INCLUDE_DIRS})

# Trova Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Usa le variabili d'ambiente di GenFit
if(DEFINED ENV{GENFIT})
    set(GENFIT_DIR $ENV{GENFIT})
    message(STATUS "Using GenFit from ${GENFIT_DIR}")
else()
    message(FATAL_ERROR "GENFIT environment variable not set. Source your env.sh script.")
endif()

# Includi le directory di GenFit
include_directories(${GENFIT_DIR}/core ${GENFIT_DIR}/include)
link_directories(${GENFIT_DIR}/lib)

# Aggiungi il file sorgente del tuo eseguibile
add_executable(myfitter myfitter.cc)

# Linka tutte le librerie di ROOT e libGeom
target_link_libraries(myfitter 
    genfit2
    Eigen3::Eigen
    ${ROOT_LIBRARIES}
    ROOT::Geom
)
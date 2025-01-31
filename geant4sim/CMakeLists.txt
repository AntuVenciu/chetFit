cmake_minimum_required(VERSION 2.8)

# Project name
project(chet_sim)

# Find Geant4 package
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif() 

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

# Include Geant4 libraries and headers
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add executable for your application
add_executable(chet_sim
    chet_sim.cc
    ${sources}
    ${headers}
)

# Link Geant4 and ROOT libraries to your application
target_link_libraries(chet_sim ${Geant4_LIBRARIES})

# find_package(ROOT REQUIRED)
target_link_libraries(chet_sim ${ROOT_LIBRARIES})

# Install targets if needed
install(TARGETS chet_sim DESTINATION bin)


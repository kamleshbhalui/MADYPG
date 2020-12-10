cmake_minimum_required(VERSION 3.4)

project(MESH2YARNS CXX)

SET(CMAKE_CXX_STANDARD 17) # NOTE: c++17 e.g. for filesystem iteration
SET(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS} -Wall -pthread")
# set(CMAKE_CXX_FLAGS "-std=c++17 -lstdc++fs") 

# Add module path in case this is project root
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../modules/" ${CMAKE_MODULE_PATH})


# common libraries
find_package(Eigen3 CONFIG REQUIRED)

add_library(tinyxml2 dependencies/tinyxml2/tinyxml2.cpp)
target_include_directories(tinyxml2 SYSTEM INTERFACE dependencies/)
target_compile_options(tinyxml2 PRIVATE -Wno-all)

add_library(bitsery INTERFACE)
target_include_directories(bitsery SYSTEM INTERFACE dependencies/)
target_compile_options(bitsery INTERFACE -Wno-all)

add_library(cnpy dependencies/cnpy/cnpy.cpp)
target_include_directories(cnpy SYSTEM INTERFACE dependencies/)
target_compile_options(cnpy PRIVATE -Wno-all)

add_library(fbx dependencies/fbx/fbxdocument.cpp
                dependencies/fbx/fbxnode.cpp
                dependencies/fbx/fbxproperty.cpp
                dependencies/fbx/fbxutil.cpp)
target_include_directories(fbx SYSTEM INTERFACE dependencies/)
target_compile_options(fbx PRIVATE -Wno-all)


# NOTE: internally commented out CubicSDFCollisionDetection.h and .cpp to avoid discregrid
##### single lib PBD
file(GLOB_RECURSE glob_PBD_SRCS
  dependencies/PositionBasedDynamics/PositionBasedDynamics/*.cpp
  dependencies/PositionBasedDynamics/Simulation/*.cpp
  dependencies/PositionBasedDynamics/Utils/*.cpp
)
add_library(pbd ${glob_PBD_SRCS})
target_include_directories(pbd
SYSTEM PUBLIC dependencies/PositionBasedDynamics/ # for "Common/Common.h" etc.
SYSTEM PUBLIC dependencies/GenericParameters/GenericParameters/) 
# # target_include_directories(pbd SYSTEM INTERFACE dependencies/)
target_link_libraries(pbd PRIVATE Eigen3::Eigen)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(pbd PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(pbd PRIVATE -Wno-all)

# target: obj2bin. object sequence to binary format converter

file(GLOB_RECURSE glob_BIN_SRCS
  io/*.cpp
  utils/*.cpp
  obj2bin.cpp
)

add_executable(obj2binary ${glob_BIN_SRCS})
target_link_libraries(obj2binary PRIVATE
  Eigen3::Eigen
  tinyxml2
  bitsery
  fbx # not actually needed, only because globbing io/*.cpp ...
  stdc++fs # <filesystem> on gcc/g++ 8.4
)

# optional cpu parallelism
set (NO_PARALLEL FALSE CACHE BOOL "no parallel")
if (${NO_PARALLEL})
  target_compile_definitions(obj2binary PUBLIC NO_PARALLEL)
  target_compile_definitions(obj2binary PUBLIC EIGEN_DONT_PARALLELIZE)
else()
  find_package(TBB CONFIG REQUIRED)
  target_link_libraries(obj2binary PUBLIC TBB::tbb)
endif()

# target: mesh2yarns ###########################################

find_package(Corrade REQUIRED Main)
find_package(Magnum REQUIRED
    GL
    MeshTools
    Shaders
    SceneGraph
    Trade
    Sdl2Application
    DebugTools
    Primitives
    # AnySceneImporter
    AnyImageImporter
    )
find_package(MagnumIntegration REQUIRED ImGui)
find_package(MagnumPlugins REQUIRED PngImporter JpegImporter)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

corrade_add_resource(Render_Rcs render/shaders/resources.conf)
corrade_add_resource(Compute_Rcs yarns/shaders/resources.conf)

file(GLOB_RECURSE glob_APP_SRCS
  arcball/*.cpp
  io/*.cpp
  mesh/*.cpp
  render/*.cpp
  utils/*.cpp
  yarns/*.cpp
  MainApplication.cpp
)

add_executable(mesh2yarns ${glob_APP_SRCS} ${Render_Rcs} ${Compute_Rcs})

target_link_libraries(mesh2yarns PRIVATE
    Corrade::Main
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Trade
    Magnum::Shaders
    Magnum::DebugTools
    Magnum::Primitives
    # Magnum::AnySceneImporter
    Magnum::AnyImageImporter
    MagnumIntegration::ImGui
    MagnumPlugins::PngImporter
    MagnumPlugins::JpegImporter
    )

# optional cpu parallelism
# compile with NO_PARALLEL fo threadutils to ignore threading
set (NO_PARALLEL FALSE CACHE BOOL "no parallel")
if (${NO_PARALLEL})
  message("Compiling without multithreading support")
  target_compile_definitions(mesh2yarns PRIVATE NO_PARALLEL)
  target_compile_definitions(mesh2yarns PRIVATE EIGEN_DONT_PARALLELIZE)
else()
  message("Compiling with multithreading support")
  find_package(TBB CONFIG REQUIRED)
  target_link_libraries(mesh2yarns PRIVATE TBB::tbb)
endif()

target_link_libraries(mesh2yarns PRIVATE
Eigen3::Eigen
tinyxml2
bitsery
cnpy
fbx
pbd
stdc++fs # <filesystem> on gcc/g++ 8.4
)



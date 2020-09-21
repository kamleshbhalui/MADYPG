cmake_minimum_required(VERSION 3.4)

project(MESH2YARNS CXX)

SET(CMAKE_CXX_STANDARD 17) # NOTE: c++17 e.g. for filesystem iteration
SET(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS} -Wall -pthread")
# set(CMAKE_CXX_FLAGS "-std=c++17 -lstdc++fs") 

# Add module path in case this is project root
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../modules/" ${CMAKE_MODULE_PATH})

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

corrade_add_resource(SSAO_Rcs render/shaders/resources.conf)


file(GLOB_RECURSE glob_SRCS
  *.cpp)

add_executable(mesh2yarns
        ${glob_SRCS}
        ${SSAO_Rcs})

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

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(mesh2yarns PUBLIC Eigen3::Eigen)

# optional cpu parallelism
# compile with NO_PARALLEL fo threadutils to ignore threading
set (NO_PARALLEL FALSE CACHE BOOL "no parallel")
if (${NO_PARALLEL})
  message("Compiling without multithreading support")
  target_compile_definitions(mesh2yarns PUBLIC NO_PARALLEL)
  target_compile_definitions(mesh2yarns PUBLIC EIGEN_DONT_PARALLELIZE)
else()
  message("Compiling with multithreading support")
  find_package(TBB CONFIG REQUIRED)
  target_link_libraries(mesh2yarns PUBLIC TBB::tbb)
endif()

# <filesystem> not working on gcc/g++ 8.4
target_link_libraries(mesh2yarns PUBLIC stdc++fs)

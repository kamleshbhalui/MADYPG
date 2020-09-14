cmake_minimum_required(VERSION 3.4)

project(MagnumSSAOExample CXX)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")

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

add_executable(magnum-ssao
        # SsaoExample.cpp
        # arcball/ArcBall.cpp
        # Shaders/SsaoShader.cpp
        # Shaders/SsaoShader.h
        # Shaders/SsaoApplyShader.cpp
        # Shaders/SsaoApplyShader.h
        ${glob_SRCS}
        ${SSAO_Rcs})

target_link_libraries(magnum-ssao PRIVATE
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
target_link_libraries(magnum-ssao PUBLIC Eigen3::Eigen)

find_package(TBB CONFIG REQUIRED)
target_link_libraries(magnum-ssao PUBLIC TBB::tbb)
cmake_minimum_required(VERSION 3.4)

project(MagnumSSAOExample CXX)

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
    )
find_package(MagnumIntegration REQUIRED ImGui)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

corrade_add_resource(SSAO_Rcs Shaders/resources.conf)


file(GLOB_RECURSE glob_SRCS
  *.cpp)

add_executable(magnum-ssao
        # SsaoExample.cpp
        # arcball/ArcBall.cpp
        # Shaders/DeferredGeometryShader.cpp
        # Shaders/DeferredGeometryShader.h
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
    MagnumIntegration::ImGui)
# INSTALLATION / COMPILATION
# ensure submodules downloaded / initialized
git submodule update --init --recursive
# build vcpkg (linux) and install dependencies
pushd vcpkg
./bootstrap-vcpkg.sh
./vcpkg install eigen3 tbb magnum[gl,meshtools,shaders,scenegraph,trade,debugtools,sdl2application] magnum-integration[imgui] magnum-plugins[pngimporter,jpegimporter]
popd
# build the project without launching/running
python exec.py mesh2yarns -r 0

  
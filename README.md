./bootstrap-vcpkg.sh
./vcpkg install eigen3 tbb magnum[gl,meshtools,shaders,scenegraph,trade,debugtools,sdl2application] magnum-integration[imgui]

or

./vcpkg install eigen3 tbb magnum[gl] magnum[meshtools] magnum[shaders] magnum[scenegraph] magnum[trade] magnum[debugtools] magnum[sdl2application] magnum-integration[imgui] magnum-plugins[pngimporter] magnum-plugins[jpegimporter]

or maybe it just always installs too much from vcpkg..


NOTE: had to comment out bottom of sdl2 portfile ... line 88 onwards

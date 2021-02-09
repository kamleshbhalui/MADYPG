./bootstrap-vcpkg.sh
./vcpkg install eigen3 tbb magnum[gl,meshtools,shaders,scenegraph,trade,debugtools,sdl2application] magnum-integration[imgui]

or

./vcpkg install eigen3 tbb magnum[gl] magnum[meshtools] magnum[shaders] magnum[scenegraph] magnum[trade] magnum[debugtools] magnum[sdl2application] magnum-integration[imgui] magnum-plugins[pngimporter] magnum-plugins[jpegimporter]

or maybe it just always installs too much from vcpkg..


NOTE: had to comment out bottom of sdl2 portfile ... line 88 onwards


dummy.fbx needed for exporting fbx geometry
zlib might be required for fbx, works on my pc, but if missing somewhere else problably can download from vcpkg and link in cmake

compile and run using:
  python exec.py 
debug build: -d, non-parallel: -p 0, don't run: -r 0 -- run instead using gdb or similar
  python exec.py -d -p 0 -r 0


TODO about yarnmapper.h main entry of method
and mainapplication.h for gui/rendering messy


TODO PBD/Discregrid/GenericParams: cmake download specific tags?

TODO how to convert a folder of obj seq to bin seq (untested remeshing animations with third parameter)

TODO about folder structure (data/)
  maybe separate readme there

TODO about model structure, pix and strain ... axes.txt
  show where created and where loaded
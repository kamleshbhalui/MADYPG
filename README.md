# Mechanics-Aware Deformation of Yarn Pattern Geometry

![Teaser Image](teaser.jpg "")
This repository contains the published code for the paper 'Mechanics Aware Deformation of Yarn Pattern Geometry' by Georg Sperl, Rahul Narain, Chris Wojtan.

Project Website: [https://visualcomputing.ist.ac.at/publications/2021/MADYPG/](https://visualcomputing.ist.ac.at/publications/2021/MADYPG/)

## Automatic Compilation & Launching

We provide a script for automatic downloading of dependencies and compilation `./automatic_compilation.sh` (tested on Ubuntu 18.04).

We also provide a selection of scripts in the `launch_scripts` folder to run some of the animations used in the paper or submission video.

Alternatively, see below for more detailed instructions instead.

## Download / Setup

This repository uses submodules that need to be cloned recursively:
```sh
git clone --recurse-submodules <REPOSITORY>
```
If cloned non-recursively, you can download the submodules afterwards using
```sh
git submodule update --init --recursive
```

Next, use the submodule [vcpkg](https://github.com/microsoft/vcpkg/) to install additional dependencies from within the `vcpkg/` directory:
```sh
./bootstrap-vcpkg.sh
./vcpkg install eigen3 tbb magnum[gl,meshtools,shaders,scenegraph,trade,debugtools,sdl2application] magnum-integration[imgui] magnum-plugins[pngimporter,jpegimporter]
```

## Compile & Run

We tested our code on Ubuntu 18.04 with g++ version > 8 for the filesystem API.
To compile and run our interactive program, we use a helper python script: `python exec.py`.

To only compile a debug build without cpu-paralleism, use e.g.: `python exec.py -d -p 0 -r 0`.

With the file `obj2bin.cpp` we provide a separate script/target that converts a folder of '.obj'-sequences into a single binary file for fast preloading. Usage: `python exec.py obj2binary IN-FOLDER OUT-FILE`

## Files/Folders

### Source Code

The code itself is located in the `src/` directory. You'll find `src/yarns/YarnMapper.h` (specifically the `step()` method) to contain the main algorithm detailed in the paper.

Loop/rendering/gui related code can be found in `MainApplication.h` and the `render/` directory.

### Data

The `data/` directory contains the following subdirectories:

- `objseqs/` contains cloth mesh animations as folders of '.obj'-sequences. (See also: `ObjSeqAnimation.h`)

- `binseqs/` contains binary single-file versions of such mesh animations for fast loading. (See also: `BinSeqAnimation.h, obj2bin.cpp`)

- `pbdsock/` contains the meshes used in our real-time sock example. (See also: `PBDSimulation.h`)

- `yarnmodels/` contains the precomputed local-displacement data for the several yarn patterns we experimented with. (See also: `Model.h`)

- `textures/` contains an assortment of matcaps (for shading yarns, or cloth/obstacle meshes), cloth textures, and the base texture for twistable ply/fiber-detail.

**Missing animations and yarnmodels:** Due to filesize/bloating we have not uploaded all mesh animations to this repository. Crucially, the sweater animations (each of which is around 180mb in binary format) are missing. Similarly, due to the filesize we did not upload the 4D bending models or the 15x15x15 and 31x31x31 yarn models. Instead we provide a direct download link to these files:
[https://doi.org/10.15479/AT:ISTA:9327](https://doi.org/10.15479/AT:ISTA:9327) 
<!-- TODO write GB if direct link.. -->

The meshes/animation have been mostly created with our previous method ["Homogenized Yarn-Level Cloth"](https://visualcomputing.ist.ac.at/publications/2020/HYLC/) (sweater animations, 30x30cm^2 stretches), or Blender.

The dummy.fbx file in the `data` directory is used for a rudimentary fbx-export of yarn or cloth geometry.

### Displacement Data Generation

The `generate_yarnmodels` folder contains the scripts we used to precompute the local-displacement data. These scripts use the python-bound optimization from our previous paper ["Homogenized Yarn-Level Cloth"](https://visualcomputing.ist.ac.at/publications/2020/HYLC/) ([Code](https://git.ist.ac.at/gsperl/HYLC)).

Note that we already provide generated data in the `data/yarnmodels/` directory, so it is not necessary to download/compile the HYLC code or rerun those scripts.

## License & Citation

This code is released under the MIT license (see [LICENSE.txt](LICENSE.txt)).
Note that some files in `src/` (files from the Magnum library and its [SSAO example](https://github.com/Janos95/magnum-examples/tree/master/src/ssao), but also `threadutils.h` and other dependencies) may have a (compatible) license comment at the top instead.

If you use our code, please consider citing our work:
```bibtex
@article{sperl2021madypg,
  author    = {Sperl, Georg and Narain, Rahul and Wojtan, Chris},
  title     = {Mechanics-Aware Deformation of Yarn Pattern Geometry},
  journal   = {ACM Transactions on Graphics (TOG)},
  number    = {4},
  volume    = {40},
  year      = {2021},
  publisher = {ACM}
}
```
<!-- @article{sperl2020hylc,
  author    = {Sperl, Georg and Narain, Rahul and Wojtan, Chris},
  title     = {Homogenized Yarn-Level Cloth},
  journal   = {ACM Transactions on Graphics (TOG)},
  number    = {4},
  volume    = {39},
  year      = {2020},
  publisher = {ACM}
} -->

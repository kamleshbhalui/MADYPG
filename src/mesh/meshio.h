#ifndef __MESHIO__H__
#define __MESHIO__H__

#include <string>
#include "Mesh.h"

Mesh load_obj_mesh(const std::string& objfile, float scale=1.0f);

#endif // __MESHIO__H__
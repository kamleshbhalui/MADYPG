#ifndef __TRAFOIO__H__
#define __TRAFOIO__H__

#include <string>

#include "../mesh/Trafo.h"

void load_xml_trafo(const std::string& path, Trafo& tf);

#endif  // __TRAFOIO__H__
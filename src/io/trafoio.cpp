#include "trafoio.h"

#include <tinyxml2/tinyxml2.h>

void load_xml_trafo(const std::string& path, Trafo& tf) {
  tinyxml2::XMLDocument doc;
  doc.LoadFile(path.c_str());
  auto el = doc.FirstChildElement("rotate");
  tf.axis_angle << atof(el->Attribute("angle")), atof(el->Attribute("x")),
      atof(el->Attribute("y")), atof(el->Attribute("z"));
  tf.axis_angle.tail<3>().normalize();
  el       = doc.FirstChildElement("scale");
  tf.scale = atof(el->FirstAttribute()->Value());
  el       = doc.FirstChildElement("translate");
  tf.translation << atof(el->Attribute("x")), atof(el->Attribute("y")),
      atof(el->Attribute("z"));
}
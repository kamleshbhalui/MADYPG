#include "Model.h"
template <>
double Model::str2num<double>(const std::string& str) { return std::strtod(str.c_str(),NULL); }
template <>
float Model::str2num<float>(const std::string& str) { return std::atof(str.c_str()); }
template <>
int Model::str2num<int>(const std::string& str) { return std::atoi(str.c_str()); }
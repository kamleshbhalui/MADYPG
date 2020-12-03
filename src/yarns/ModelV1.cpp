#include "ModelV1.h"
template <>
double ModelV1::str2num<double>(const std::string& str) { return std::strtod(str.c_str(),NULL); }
template <>
float ModelV1::str2num<float>(const std::string& str) { return std::atof(str.c_str()); }
template <>
int ModelV1::str2num<int>(const std::string& str) { return std::atoi(str.c_str()); }
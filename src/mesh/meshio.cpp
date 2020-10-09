// #include "meshio.h"

// #include <fstream>
// #include <iostream>
// #include <iterator>
// #include <sstream>

// #include "../utils/debug_logging.h"
// #include "../utils/threadutils.h"

// int count_words_(const std::string &str);
// int count_words_(const std::string &str) {
//   std::istringstream iss(str);
//   return std::distance(std::istream_iterator<std::string>(iss),
//                        std::istream_iterator<std::string>());
// }

// // int count_remaining_words(std::stringstream &ss)
// // {
// //   auto spos = ss.tellg();
// //   auto state = ss.rdstate();
// //   int count = std::distance(std::istream_iterator<std::string>(ss),
// //   std::istream_iterator<std::string>()); ss.seekg(spos); ss.setstate(state);
// //   std::cout<<state<<spos<<"\n";
// //   return count;
// // }

// void load_obj_mesh(const std::string &filename, Mesh &mesh,
//                    bool with_uv_or_topology, float scale) {
//   std::deque<VectorGLi> f;      // polygonal faces (world-space)
//   std::deque<VectorGLi> fms;    // polygonal faces (material-space)
//   AlignedDeque<VectorGL3f> v;   // world coordinates
//   AlignedDeque<VectorGL2f> vt;  // texture coordinates / material space

//   std::fstream ifs(filename.c_str(), std::ios::in);
//   if (!ifs) {
//     std::cerr << "Error: failed to open file " << filename << "\n";
//     return;
//   }

//   std::string line;
//   std::string kw;
//   while (std::getline(ifs, line)) {
//     // trim(line);
//     if (line.empty())
//       continue;

//     std::stringstream linestream(line);
//     linestream >> kw;

//     if (kw == "vt" && with_uv_or_topology) {
//       VectorGL2f vec;
//       linestream >> vec[0] >> vec[1];
//       vt.push_back(vec * scale);
//     } else if (kw == "v") {
//       VectorGL3f vec;
//       linestream >> vec[0] >> vec[1] >> vec[2];
//       v.push_back(vec * scale);
//     } else if (kw == "f" && with_uv_or_topology) {
//       int nprimverts = count_words_(line) - 1;
//       VectorGLi ixs_ms(nprimverts);
//       VectorGLi ixs_ws(nprimverts);
//       std::string w;
//       int j = 0;
//       while (linestream >> w) {
//         std::stringstream wstream(w);
//         int msix, wsix;
//         char c;
//         wstream >> wsix >> c >> msix;
//         if (wstream.fail()) {  // couldn't read format 'ws/ms'. trying instead
//                                // just 'ws', ie only one set of face indices /
//                                // no seam vertices
//           wstream.clear();
//           wstream.str(w);
//           wstream >> wsix;
//           msix = wsix;
//           Debug::msgassert(
//               "OBJ: Couldn't read face indices: '" + line + "'",
//               !wstream.fail());
//         }
//         ixs_ms[j] = msix - 1;  // NOTE: obj files are 1-indexed
//         ixs_ws[j] = wsix - 1;
//         ++j;
//       }
//       assert(j == nprimverts);

//       Debug::msgassert("OBJ: polygonal face triangulation not implemented",
//                        nprimverts <= 4);

//       if (nprimverts == 3) {
//         f.push_back(ixs_ws);
//         if (with_uv_or_topology)
//           fms.push_back(ixs_ms);
//       } else if (nprimverts == 4) {
//         // triangulate
//         bool verts_loaded = true;
//         for (size_t k = 0; k < 3; ++k)
//           if (ixs_ms[k] >= vt.size()) {
//             verts_loaded = false;
//             break;
//           }
//         VectorGLi f0, f1;
//         if (verts_loaded) {  // find shortest diagonal
//           scalar diag02 = (vt[ixs_ms[0]] - vt[ixs_ms[2]]).squaredNorm();
//           scalar diag13 = (vt[ixs_ms[1]] - vt[ixs_ms[3]]).squaredNorm();
//           if (diag02 <= diag13) {
//             fms.emplace_back(3);
//             fms.back() << ixs_ms[0], ixs_ms[1], ixs_ms[2];
//             fms.emplace_back(3);
//             fms.back() << ixs_ms[0], ixs_ms[2], ixs_ms[3];
//             f.emplace_back(3);
//             f.back() << ixs_ws[0], ixs_ws[1], ixs_ws[2];
//             f.emplace_back(3);
//             f.back() << ixs_ws[0], ixs_ws[2], ixs_ws[3];
//           } else {
//             fms.emplace_back(3);
//             fms.back() << ixs_ms[0], ixs_ms[1], ixs_ms[3];
//             fms.emplace_back(3);
//             fms.back() << ixs_ms[1], ixs_ms[2], ixs_ms[3];
//             f.emplace_back(3);
//             f.back() << ixs_ws[0], ixs_ws[1], ixs_ws[3];
//             f.emplace_back(3);
//             f.back() << ixs_ws[1], ixs_ws[2], ixs_ws[3];
//           }
//         } else {  // fallback to fixed arbitrary (if vt not specified before f)
//           fms.emplace_back(3);
//           fms.back() << ixs_ms[0], ixs_ms[1], ixs_ms[2];
//           fms.emplace_back(3);
//           fms.back() << ixs_ms[0], ixs_ms[2], ixs_ms[3];
//           f.emplace_back(3);
//           f.back() << ixs_ws[0], ixs_ws[1], ixs_ws[2];
//           f.emplace_back(3);
//           f.back() << ixs_ws[0], ixs_ws[2], ixs_ws[3];
//         }
//       }
//     }
//   }

//   // allocate
//   mesh.X.resize(v.size(), 3);
//   if (with_uv_or_topology)
//     mesh.U.resize(vt.size(), 2);
//   if (f.size() > 0) {
//     mesh.F.resize(f.size(), f.back().size());
//     if (with_uv_or_topology)
//       mesh.Fms.resize(fms.size(), f.back().size());
//   }

//   // fill
//   threadutils::parallel_for(size_t(0), v.size(),
//                             [&](size_t i) { mesh.X.row(i) = v[i]; });
//   if (with_uv_or_topology) {
//     threadutils::parallel_for(size_t(0), vt.size(),
//                               [&](size_t i) { mesh.U.row(i) = vt[i]; });
//     threadutils::parallel_for(size_t(0), fms.size(),
//                               [&](size_t i) { mesh.Fms.row(i) = fms[i]; });
//     // NOTE: TODO currently this method is called also for obstacle meshes where we don't need that guarantee
//     //   have to make it so it only checks for cloth meshes :/
//     // Debug::msgassert(
//     //     "OBJ: v and vt data size inconsistent! " + std::to_string(v.size()) + " !<= " + std::to_string(vt.size()), 
//     //     v.size() <=
//     //         vt.size());  // NOTE: one worldspace vertex can have multiple uvs
//     threadutils::parallel_for(size_t(0), f.size(),
//                               [&](size_t i) { mesh.F.row(i) = f[i]; });
//   }

//   // TODO optionally: assert each index in faces is in [0,nvertices)
// }

// // Mesh load_obj_mesh(const std::string &filename, float scale) {
// //   Mesh mesh;
// //   load_obj_mesh(filename, mesh, true, scale);
// //   return mesh;
// // }
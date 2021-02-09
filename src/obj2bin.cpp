
#include <algorithm>
#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)
#include <regex>
namespace fs = std::filesystem;

#include <bitsery/adapter/stream.h>
#include <bitsery/bitsery.h>
#include <bitsery/traits/vector.h>
#include "io/framesio.h"
#include "io/objio.h"
#include "io/trafoio.h"

// This file is a script to convert obj-sequences into our custom binary format
// for fast preloading of mesh animations.
// It basically follows the same principle as the 'ObjSeqAnimation',
// with two modes for general cloth sims and ARCsim, due to naming conventions
// and obstacle data
// Compile & use:
//    python exec.py obj2binary IN-FOLDER OUT-FILE [constuv=1]

// #include "io/bitsery_eigen.h"
std::vector<Frame> load_folder(const std::string& folder, bool const_uv=true) {
  std::vector<Frame> frames;

  // heuristic check for arcsim simulation, with specific file names etc.
  bool arcsim_mode = fs::exists(fs::path(folder) / "conf.json");

  if (!arcsim_mode) {
    static const auto regex_obj = std::regex(".*\\.obj$");

    std::vector<std::string> files;
    for (auto& p : fs::directory_iterator(folder)) {
      if (fs::is_regular_file(p) &&
          std::regex_match(p.path().filename().string(), regex_obj))
        files.push_back(p.path().filename());
      else if (fs::is_directory(p) && p.path().filename() == "obs") {
        frames.resize(std::max(int(frames.size()), 1));
        for (auto& p2 : fs::directory_iterator(p)) {
          if (fs::is_regular_file(p2) &&
              std::regex_match(p2.path().filename().string(), regex_obj)) {
            frames[0].obs_V.emplace_back();
            frames[0].obs_F.emplace_back();
            load_obj(p2.path(), frames[0].obs_V.back(), frames[0].obs_F.back());
          }
        }
      }
    }
    std::sort(files.begin(), files.end());
    frames.resize(files.size());
    for (size_t i = 0; i < files.size(); ++i) {
      std::string fpath = fs::path(folder) / files[i];
      if (i == 0u || !const_uv)  // load entire mesh for first frame
        load_obj(fpath, frames[i].cloth_V, frames[i].cloth_F, frames[i].cloth_U,
                 frames[i].cloth_Fms);
      else  // load only ws vertex updates, assuming no changes to UVs or
            // topology
        load_obj<MatrixGLf,MatrixGLi>(fpath, frames[i].cloth_V);
    }

    // obstacles
    for (auto& p : fs::directory_iterator(folder)) {
      if (fs::is_directory(p) && p.path().filename() == "obs") {
        if (frames.size() < 1u)
          frames.resize(1);
        for (auto& p2 : fs::directory_iterator(p)) {
          if (fs::is_regular_file(p2) &&
              std::regex_match(p2.path().filename().string(), regex_obj)) {
            frames[0].obs_V.emplace_back();
            frames[0].obs_F.emplace_back();
            load_obj(p2.path(), frames[0].obs_V.back(), frames[0].obs_F.back());
          }
        }
      }
    }

  } else {
    static const auto regex_clothobj = std::regex("(\\d{4})_00\\.obj");
    static const auto regex_obstrafo = std::regex("(\\d{4})obs(\\d{2})\\.txt");
    static const auto regex_obsobj   = std::regex("obs_(\\d{2})\\.obj");

    for (auto& p : fs::directory_iterator(folder)) {
      if (!fs::is_regular_file(p))
        continue;

      const auto& str = p.path().filename().string();
      std::smatch match;
      if (std::regex_search(str, match, regex_clothobj)) {
        size_t i = std::stoul(match[1].str());
        if (i >= frames.size())
          frames.resize(i + 1);
        auto& frame = frames[i];
        if (i == 0u || !const_uv)  // load entire mesh for first frame
          load_obj(p.path(), frame.cloth_V, frame.cloth_F, frame.cloth_U,
                   frame.cloth_Fms);
        else  // load only ws vertex updates, assuming no changes to UVs or
              // topology
          load_obj<MatrixGLf,MatrixGLi>(p.path(), frame.cloth_V);
      } else if (std::regex_search(str, match, regex_obstrafo)) {
        size_t i = std::stoul(match[1].str());
        if (i >= frames.size())
          frames.resize(i + 1);
        size_t id = std::stoul(match[2].str());

        auto& frame = frames[i];

        if (id >= frame.obs_trafo.size())
          frame.obs_trafo.resize(id + 1);

        load_xml_trafo(p.path(), frame.obs_trafo[id]);
      } else if (std::regex_search(str, match, regex_obsobj)) {
        size_t id = std::stoul(match[1].str());
        if (frames.size() < 1u)
          frames.resize(1);
        auto& frame = frames[0];
        if (id >= frame.obs_V.size())
          frame.obs_V.resize(id + 1);
        if (id >= frame.obs_F.size())
          frame.obs_F.resize(id + 1);
        load_obj(p.path(), frame.obs_V[id], frame.obs_F[id]);
      }
    }
  }

  return frames;
}

bool convert_folder(const std::string& folder, const std::string& outfile, bool const_uv = false) {
  std::cout << "Loading folder: " << folder << "\n";
  std::vector<Frame> frames = load_folder(folder, const_uv);
  std::cout << "Loaded folder.\n";

  std::fstream s(outfile, std::ios::binary | std::ios::trunc | std::ios::out);
  if (!s.is_open()) {
    std::cerr << "Cannot open " << outfile << " for writing\n";
    return false;
  }

  std::cout << "Serializing to: " << outfile << "\n";
  bitsery::Serializer<bitsery::OutputBufferedStreamAdapter> ser{s};
  ser.object(frames);
  ser.adapter().flush();
  s.close();
  std::cout << "Done.\n";

  return true;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Not enough arguments. Expecting ./obj2binary in-folder out-file [constuv=1].\n";
    return 1;
  }

  // assuming: constant number of obstacles, data existing for each rame, ...

  bool const_uv = true; // false if it's a remeshing sim
  if (argc >= 4) {
    const_uv = bool(atof(argv[3]));
  }
  if (!const_uv) {
    std::cout<< "Assuming remeshed simulation, each frame will get updated topology.\n";
  }

  if (!convert_folder(argv[1], argv[2], const_uv))
    return 2;

  return 0;
}

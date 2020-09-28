#ifndef __OBJSEQANIMATION__H__
#define __OBJSEQANIMATION__H__

#include <string>

#include "AbstractMeshProvider.h"

// TODO CPP includes
#include <algorithm>
#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)
#include <regex>

#include "../tinyxml2/tinyxml2.h"
#include "../utils/debug_logging.h"
#include "meshio.h"
// namespace fs = std::experimental::filesystem;
namespace fs = std::filesystem;

// IN: folder, repeat, constant_material_space(=topology&positions)
// OUT: update (sets mesh and flags)

class ObjSeqAnimation : public AbstractMeshProvider {
 public:
  struct Settings {
    std::string folder;
    bool repeat                  = true;
    bool constant_material_space = false;
  } m_settings;

  // TODO PUT SOMEWHERE
  bool endsWith(const std::string& mainStr, const std::string& toMatch) {
    if (mainStr.size() >= toMatch.size() &&
        mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(),
                        toMatch) == 0)
      return true;
    else
      return false;
  }

  ObjSeqAnimation(const Settings& settings) : m_settings(settings) {
    m_indicesDirty = true;

    if (!fs::is_directory(m_settings.folder)) {  // keep uninitialized
      Debug::error("Could not load mesh sequence directory:",
                   m_settings.folder);
      return;
    }

    // heuristic check for arcsim simulation, with specific file names etc.
    arcsim_mode = fs::exists(fs::path(m_settings.folder) / "conf.json");

    // get sorted list of files
    m_iter = 0;
    if (!arcsim_mode) {
      for (auto& p : fs::directory_iterator(m_settings.folder)) {
        if (fs::is_regular_file(p))
          m_files.push_back(p.path());
      }
      std::sort(m_files.begin(), m_files.end());
    } else {
      static const auto regex_clothobj = std::regex("(\\d{4})_00\\.obj");
      static const auto regex_obstrafo =
          std::regex("(\\d{4})obs(\\d{2})\\.txt");
      static const auto regex_obsobj = std::regex("obs_(\\d{2})\\.obj");

      for (auto& p : fs::directory_iterator(m_settings.folder)) {
        if (!fs::is_regular_file(p))
          continue;

        const auto& str = p.path().filename().string();
        std::smatch match;
        if (std::regex_search(str, match, regex_clothobj)) {
          int frame = std::stoi(match[1].str());
          if (frame >= int(m_files.size()))
            m_files.resize(frame + 1);
          m_files[frame] = p.path();
        } else if (std::regex_search(str, match, regex_obstrafo)) {
          int frame = std::stoi(match[1].str());
          int id    = std::stoi(match[2].str());

          if (id >= int(m_obstacle_trafos.size()))
            m_obstacle_trafos.resize(id + 1);
          auto& data = m_obstacle_trafos[id];
          if (frame >= int(data.size()))
            data.resize(frame + 1);

          auto& tf = data[frame];

          // TODO function this
          tinyxml2::XMLDocument doc;
          doc.LoadFile(p.path().string().c_str());
          auto el = doc.FirstChildElement("rotate");
          tf.axis_angle << atof(el->Attribute("angle")),
              atof(el->Attribute("x")), atof(el->Attribute("y")),
              atof(el->Attribute("z"));
          tf.axis_angle.tail<3>().normalize();
          el       = doc.FirstChildElement("scale");
          tf.scale = atof(el->FirstAttribute()->Value());
          el       = doc.FirstChildElement("translate");
          tf.translation << atof(el->Attribute("x")), atof(el->Attribute("y")),
              atof(el->Attribute("z"));
        } else if (std::regex_search(str, match, regex_obsobj)) {
          int id = std::stoi(match[1].str());
          if (id >= int(m_obstacles.size()))
            m_obstacles.resize(id + 1);
          load_obj_mesh(p.path(), m_obstacles[id].mesh, true);
        }
      }
    }

    update();  // load first frame

    // initial flags
    m_indicesDirty = true;
  }

  ~ObjSeqAnimation() {}

  void update() {
    if (m_files.empty()) {
      m_indicesDirty = false;
      return;
    }

    // repeat
    if (m_iter >= m_files.size() && m_settings.repeat) {
      m_iter = 0;
    }

    // load next
    if (m_iter < m_files.size()) {
      // load mesh, with material space data if non-constant or uninitialized
      bool load_uv =
          !m_settings.constant_material_space || m_mesh.Fms.rows() == 0;
      m_indicesDirty = load_uv;
      load_obj_mesh(m_files[m_iter], m_mesh, load_uv);

      if (arcsim_mode) {
        // update obstacle transformations
        Debug::msgassert("Inconsistent obstacle and transformation ids",
                         m_obstacles.size() == m_obstacle_trafos.size());
        for (size_t i = 0; i < m_obstacle_trafos.size(); ++i) {
          Debug::msgassert("Inconsistent frames in cloth and ostacles",
                           m_iter < m_obstacle_trafos[i].size());
          m_obstacles[i].transformation = m_obstacle_trafos[i][m_iter];
        }
      }

      ++m_iter;
    }
  }

 private:
  std::vector<std::string> m_files;
  std::vector<std::deque<AbstractMeshProvider::Obstacle::Trafo>>
      m_obstacle_trafos;
  size_t m_iter;

  bool arcsim_mode;
};

#endif  // __OBJSEQANIMATION__H__

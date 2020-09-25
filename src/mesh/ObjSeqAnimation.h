#ifndef __OBJSEQANIMATION__H__
#define __OBJSEQANIMATION__H__

#include <string>

#include "AbstractMeshProvider.h"

// TODO CPP includes
#include <algorithm>
#include <filesystem>  // requires C++17, g++ >= 8 & linking to stdc++fs)

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

  ObjSeqAnimation(const Settings& settings) : m_settings(settings) {
    m_indicesDirty = true;

    if (!fs::is_directory(m_settings.folder)) {  // keep uninitialized
      Debug::error("Could not load mesh sequence directory:",
                   m_settings.folder);
      return;
    }

    // get sorted list of files
    m_iter = 0;
    for (auto& p : fs::directory_iterator(m_settings.folder)) {
      m_files.push_back(p.path());
    }
    std::sort(m_files.begin(), m_files.end());

    update();  // load first

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

      ++m_iter;
    }
  }

 private:
  std::vector<std::string> m_files;
  size_t m_iter;
};

#endif  // __OBJSEQANIMATION__H__

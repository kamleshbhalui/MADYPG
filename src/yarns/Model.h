#ifndef __MODEL__H__
#define __MODEL__H__

#include "PeriodicYarnPattern.h"
#include "YarnSoup.h"  // vertexmsdata

class Model {
 public:
  Model(const std::string& folder);

  Vector4s deformation(const Vector6s& strain, uint32_t pix) const;

  const PeriodicYarnPattern& getPYP() const { return m_pyp; }
  bool isInitialized() const { return m_initialized; }
  auto& getTexAxes() { return m_tex_sxsasy_axes; }
  auto& getTexData() { return m_tex_sxsasy_data; }

 private:
  bool m_initialized;
  PeriodicYarnPattern m_pyp;

  struct DeformationEntry {
    float x, y, z, th;
    Eigen::Map<Eigen::Matrix<float, 4, 1, Eigen::ColMajor>, Eigen::Unaligned>
    map() {
      return Eigen::Map<Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<float*>(this), 4);
    }
    Eigen::Map<const Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
               Eigen::Unaligned>
    map() const {
      return Eigen::Map<const Eigen::Matrix<float, 4, 1, Eigen::ColMajor>,
                        Eigen::Unaligned>(reinterpret_cast<const float*>(this),
                                          4);
    }
  };

#define AXES_MAX_LENGTH 32
  struct AxesInfo {
    uint32_t lenSX, lenSA, lenSY;
    float SX[AXES_MAX_LENGTH];
    float SA[AXES_MAX_LENGTH];
    float SY[AXES_MAX_LENGTH];
    float invSX[AXES_MAX_LENGTH];
    float invSA[AXES_MAX_LENGTH];
    float invSY[AXES_MAX_LENGTH];
  };

  struct axwrapper {  // helper struct to pass arrays in the 'AxesInfo' uniform
                      // buffer into the file loading method
    uint32_t& len;
    float (&data)[AXES_MAX_LENGTH];  // ref of fixes size array..
    float (&invdata)[AXES_MAX_LENGTH];
    axwrapper(uint32_t& lenref, float (&dataref)[AXES_MAX_LENGTH],
              float (&invdataref)[AXES_MAX_LENGTH])
        : len(lenref), data(dataref), invdata(invdataref) {}
  };

  bool load_axes(const std::string& filepath, std::vector<axwrapper>& axes);

  Vector4s sample3D(Vector3s strain, uint32_t pix) const;

  VectorBuffer<AxesInfo>
      m_tex_sxsasy_axes;  // note: axinfo exploiting "vector"buffer for just
                          // single entry
  VectorBuffer<DeformationEntry> m_tex_sxsasy_data;
};

#endif  // __MODEL__H__
#ifndef __MODEL__H__
#define __MODEL__H__

#include "PeriodicYarnPattern.h"



class Model {
 public:
  Model(const std::string& folder) {
    // Load model / pyp
    m_pyp.deserialize(folder + "/pyp");  // DEBUG hardcoded file
    m_pyp.rectangulize();  // TODO potentially assume that this is true for new
    // pyp, but it should take only a millisecond anyway

    // TODO compute periodic yarns

    // TODO compute vix2y, vix2t // TODO cache in pyp?

    // TODO compute (y,ix)->t ??? // TODO cache ?

    // TODO load sx sy data like python
  }
  
  const Vector4s deformation(const Vector6s& strain, int vix) {
    Vector4s g = Vector4s::Zero();

    // TODO vix -> y,t, xref=pyp.Q.row(vix)

    // TODO
    // int i0=0,i1=2;
    // x = lerp2D[y].blerp(strain[i0],strain[i1], t)
    // g += x - xref
    // g += 
    DECLARE_UNUSED(vix);
    DECLARE_UNUSED(strain);

    return g;
  }

  const PeriodicYarnPattern& getPYP() const { return m_pyp; }

 private:
  PeriodicYarnPattern m_pyp;

  class Lerp2D{
    // ...
  };
};

#endif  // __MODEL__H__
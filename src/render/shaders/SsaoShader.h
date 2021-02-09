#ifndef Magnum_Examples_Ssao_SsaoShader_h
#define Magnum_Examples_Ssao_SsaoShader_h

#include <Corrade/Containers/Array.h>

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/GL/MultisampleTexture.h>

#include "../render_definitions.h"

#include <Magnum/Math/TypeTraits.h>
#include <type_traits>
#include <Magnum/Tags.h>

// fmod from more up-to-date version of Magnum (than the one in vcpkg)
namespace Magnum {
namespace Math {

template <class T>
inline typename std::enable_if<IsScalar<T>::value, T>::type fmod(T a, T b) {
  return T(std::fmod(UnderlyingTypeOf<T>(a), UnderlyingTypeOf<T>(b)));
}

template <std::size_t size, class T>
inline Vector<size, T> fmod(const Vector<size, T>& a,
                            const Vector<size, T>& b) {
  Vector<size, T> out{Magnum::NoInit};
  for (std::size_t i = 0; i != size; ++i) out[i] = Math::fmod(a[i], b[i]);
  return out;
}
}
}

namespace Magnum {

class SsaoShader : public GL::AbstractShaderProgram {
public:
    enum {
        AmbientOcclusionOutput = 0,
    };

    explicit SsaoShader(UnsignedInt sampleCount = 64);

    explicit SsaoShader(Magnum::NoCreateT) : GL::AbstractShaderProgram{Magnum::NoCreate} {};

    #ifdef MSAA
    SsaoShader& bindPositionTexture(GL::MultisampleTexture2D&);
    SsaoShader& bindNormalTexture(GL::MultisampleTexture2D&);
    #else
    SsaoShader& bindPositionTexture(GL::Texture2D&);
    SsaoShader& bindNormalTexture(GL::Texture2D&);
    #endif

    SsaoShader& bindNoiseTexture(GL::Texture2D&);
    SsaoShader& bindOcclusionTexture(GL::Texture2D&);
    SsaoShader& setProjectionMatrix(Matrix4 const&);
    SsaoShader& setSampleRadius(Float);
    SsaoShader& setBias(Float);

private:

    UnsignedInt _projectionUniform = 0,
                _biasUniform = 1,
                _radiusUniform = 2,
                _samplesUniform = 3;
};

}

#endif
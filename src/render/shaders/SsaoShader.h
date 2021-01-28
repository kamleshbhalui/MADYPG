#ifndef Magnum_Examples_Ssao_SsaoShader_h
#define Magnum_Examples_Ssao_SsaoShader_h

#include "../../magnum-master-changes.h"
#include <Corrade/Containers/Array.h>

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/GL/MultisampleTexture.h>

#include "../render_definitions.h"

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
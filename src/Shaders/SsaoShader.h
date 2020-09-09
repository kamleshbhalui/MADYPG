#ifndef Magnum_Examples_Ssao_SsaoShader_h
#define Magnum_Examples_Ssao_SsaoShader_h

/*
    This file is part of Magnum.

    Original authors — credit is appreciated but not required:

        2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 —
            Vladimír Vondruš <mosra@centrum.cz>
        2020 — Janos Meny <janos.meny@gmail.com>

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or distribute
    this software, either in source code form or as a compiled binary, for any
    purpose, commercial or non-commercial, and by any means.

    In jurisdictions that recognize copyright laws, the author or authors of
    this software dedicate any and all copyright interest in the software to
    the public domain. We make this dedication for the benefit of the public
    at large and to the detriment of our heirs and successors. We intend this
    dedication to be an overt act of relinquishment in perpetuity of all
    present and future rights to this software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "../magnum-master-changes.h"
#include <Corrade/Containers/Array.h>

#include <Magnum/GL/AbstractShaderProgram.h>
#include <Magnum/Shaders/Generic.h>
#include <Magnum/Math/Vector2.h>

namespace Magnum { namespace Examples {

class SsaoShader : public GL::AbstractShaderProgram {
public:
    enum {
        AmbientOcclusionOutput = 0,
    };

    enum class Flag : UnsignedInt {
        ComputeShader,
        FragmentShader
    };

    explicit SsaoShader(Flag flag = Flag::FragmentShader, UnsignedInt sampleCount = 64);

    explicit SsaoShader(Magnum::NoCreateT) : GL::AbstractShaderProgram{Magnum::NoCreate} {};

    SsaoShader& bindPositionTexture(GL::Texture2D&);

    SsaoShader& bindNormalTexture(GL::Texture2D&);

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

}}

#endif
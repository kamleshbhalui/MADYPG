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

in vec2 textureCoordinate;

uniform sampler2D positions;
uniform sampler2D normals;
uniform sampler2D noise;

uniform mat4 projection;
uniform float bias = 0.025;//0.2; // TODO UNIFORM
uniform float radius = 0.5;//1.; // TODO UNIFORM 
uniform vec3 samples[SAMPLE_COUNT];

out float ambientOcclusion;

vec2 noiseScale = vec2(textureSize(positions,0))/vec2(textureSize(noise,0));

void main()
{
    vec3 position = texture(positions, textureCoordinate).xyz;
    if(position.z == 0) {
        ambientOcclusion = 1;
        return;
    }
    vec3 normal = normalize(texture(normals, textureCoordinate).xyz);
    vec3 randomVector = normalize(texture(noise, textureCoordinate*noiseScale).xyz);

    /* tangent-space to view-space */
    vec3 tangent = normalize(randomVector - normal*dot(randomVector, normal));
    vec3 bitangent = cross(normal, tangent);
    mat3 TBN = mat3(tangent, bitangent, normal);

    /* iterate over the sample kernel and calculate occlusion factor */
    float occlusion = 0.0;
    for(int i = 0; i < SAMPLE_COUNT; ++i) {
        vec3 randomSample = TBN*samples[i];
        randomSample = position + randomSample * radius;

        vec4 sampleNDC = projection*vec4(randomSample, 1.);
        sampleNDC.xyz /= sampleNDC.w; /* perspective division */
        vec2 sampleCoords = sampleNDC.xy*0.5 + 0.5;
        float sampleDepth = texture(positions, sampleCoords).z;

        /* range check & accumulate */
        float rangeCheck = smoothstep(0.0, 1.0, radius/abs(position.z - sampleDepth));
        occlusion += (sampleDepth >= randomSample.z + bias ? 1.0 : 0.0)*rangeCheck;
    }

    ambientOcclusion = 1.0 - (occlusion / float(SAMPLE_COUNT));
    // ambientOcclusion = 1.0;
}

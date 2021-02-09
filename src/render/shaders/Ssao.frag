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

#ifdef MSAA // defined in enableMSAA.h
uniform sampler2DMS positions; 
uniform sampler2DMS normals;

// average sampling
#define invMSAA 1.0/MSAA
vec4 textureMSavg(sampler2DMS tex, ivec2 P) {
    // return texelFetch(tex, P, 0);
    vec4 val = vec4(0.0);
    for (int i = 0; i < MSAA; i++)
        val += texelFetch(tex, P, i);
    return val * invMSAA;
}
vec4 textureMSsum(sampler2DMS tex, ivec2 P) {
    // return texelFetch(tex, P, 0);
    vec4 val = vec4(0.0);
    for (int i = 0; i < MSAA; i++)
        val += texelFetch(tex, P, i);
    return val;
}
#else
uniform sampler2D positions;
uniform sampler2D normals;
#endif

uniform sampler2D noise;

uniform mat4 projection;
uniform float bias = 0.025;//0.2;
uniform float radius = 0.5;//1.;
uniform vec3 samples[SAMPLE_COUNT];

out float ambientOcclusion;

#ifdef MSAA
vec2 noiseScale = vec2(textureSize(positions))/vec2(textureSize(noise,0));
#else
vec2 noiseScale = vec2(textureSize(positions,0))/vec2(textureSize(noise,0));
#endif

void main()
{
    // NOTE modified and imperfect SSAO implementation
    // splitting the SSAO samples over the MSAA samples to maintain constant number
    #ifdef MSAA
    float occlusion = 0.0;
    // ambientOcclusion = 0.0;
    ivec2 texsize = textureSize(positions); // assume same for normals
    int N_AO = SAMPLE_COUNT / MSAA; // assume divisible
    float invN = 1.0 / float(N_AO);
    for (int m = 0; m < MSAA; m++) {
        ivec2 P = ivec2(textureCoordinate * texsize);
        vec3 position = texelFetch(positions, P, m).xyz;
        vec3 normal = normalize(texelFetch(normals, P, m).xyz);

         // NOTE: using arbitrary -1000 as background depth
        if(position.z < -999.9) {
            continue; // for multisampling reject any sample that happens to land in the background
        }

        // if(position.z == 0) {
        //     ambientOcclusion = 1;
        //     return;
        // }
        vec3 randomVector = normalize(texture(noise, textureCoordinate*noiseScale).xyz);

        /* tangent-space to view-space */
        vec3 tangent = normalize(randomVector - normal*dot(randomVector, normal));
        vec3 bitangent = cross(normal, tangent);
        mat3 TBN = mat3(tangent, bitangent, normal);

        /* iterate over the sample kernel and calculate occlusion factor */
        // float occlusion = 0.0;
        for(int i = 0; i < N_AO; ++i) {
            vec3 randomSample = TBN*samples[m*N_AO + i];
            randomSample = position + randomSample * radius;

            vec4 sampleNDC = projection*vec4(randomSample, 1.);
            sampleNDC.xyz /= sampleNDC.w; /* perspective division */
            vec2 sampleCoords = sampleNDC.xy*0.5 + 0.5;

            float sampleDepth = texelFetch(positions, ivec2(sampleCoords * texsize), m).z;
            // float sampleDepth = texture(positions, sampleCoords).z;

            /* range check & accumulate */
            // float rangeCheck = smoothstep(0.0, 1.0, radius/abs(position.z - sampleDepth));
            float rangeCheck = 1.0; // do this to shadow at any distance

            occlusion += (sampleDepth >= randomSample.z + bias ? 1.0 : 0.0)*rangeCheck;
        }
    }

    ambientOcclusion = 1.0 - (occlusion * invN * invMSAA);

    #else

    vec3 position = texture(positions, textureCoordinate).xyz;
    // if(position.z == 0) {
    //     ambientOcclusion = 1;
    //     return;
    // }
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
    #endif
}

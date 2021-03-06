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

#ifdef MSAA
uniform sampler2DMS positionTexture; 
uniform sampler2DMS albedoTexture;
#define invMSAA 1.0/MSAA
vec4 textureMSavg(sampler2DMS sampler, ivec2 P) {
    vec4 val = vec4(0.0);
    for (int i = 0; i < MSAA; i++)
        val += texelFetch(sampler, P, i);
    return val * invMSAA;
}
#else
uniform sampler2D positionTexture;
uniform sampler2D albedoTexture;
#endif
uniform sampler2D ambientOcclusionTexture;

// uniform sampler2D normalTexture;
// uniform vec3 lightPosition;
// uniform vec3 lightColor;
// uniform float shininess;
// uniform vec3 specularColor;

out vec4 fragmentColor;

uniform int ao_blur_radius = 2;
uniform float ao_blur_feature = 1.0; // exp(-dz^2/2sigmafeature^2) = exp(-dz^2/2 blur_feature^2)
uniform float ao_pow = 4.0;

// bilateral filter:
// If(xi) = 1/Wp sum_j ( I(xj) fr(|I(xi)-I(xj)|) gs(|xi-xj|) )
// Wp = sum_j ( fr(|I(xi)-I(xj)|) gs(|xi-xj) )

// joint bilateral filter: ubstead if using main signal difference in weight, use other signal's dif
// If(xi) = 1/W sum_j ( I(xj) w(xi,xj) )
//          w: wx(i,j) wy(zi,zj)  instead of wy(xi,xj)
//https://bartwronski.com/2019/09/22/local-linear-models-guided-filter/
// w: exp( - (dij)^2/(2 sigd^2) - dz^2/(2 sigr^2) ) 

float bilateral_weight() {
    float w = 0.0;
    return w;
}

// https://stackoverflow.com/a/17897228
vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{

#ifdef MSAA
    float blurredOcclusion = texture(ambientOcclusionTexture, textureCoordinate).x;
    blurredOcclusion = pow(blurredOcclusion, ao_pow);
    #ifndef DRAW_OCCLUSION_FACTOR 
    vec3 albedo = textureMSavg(albedoTexture, ivec2(textureCoordinate*textureSize(albedoTexture))).rgb;
    fragmentColor.rgb = albedo*blurredOcclusion;
    #else
    fragmentColor.rgb = vec3(blurredOcclusion);
    #endif

    return; // MSAA blur not implemented ...
#else


    /* retrieve data from gbuffer */
    vec3 position = texture(positionTexture, textureCoordinate).rgb;

#ifndef DRAW_OCCLUSION_FACTOR
    // vec3 normal = texture(normalTexture, textureCoordinate).rgb;
    vec3 albedo = texture(albedoTexture, textureCoordinate).rgb;
#endif



    // /* average blur occlusion factor */
    vec2 texelSize = 1.0 / vec2(textureSize(ambientOcclusionTexture, 0));
    // float blurredOcclusion = 0.0;
    // for (int x = -ao_blur_radius; x <= ao_blur_radius; ++x) {
    //     for (int y = -ao_blur_radius; y <= ao_blur_radius; ++y) {
    //         vec2 offset = vec2(float(x), float(y)) * texelSize;
    //         blurredOcclusion += texture(ambientOcclusionTexture, textureCoordinate + offset).x;
    //     }
    // }
    // blurredOcclusion /= float((ao_blur_radius*2 + 1)*(ao_blur_radius*2 + 1));

    float blurredOcclusion = 0.0;

    if (ao_blur_radius > 0) {
        float wtotal = 0.0;

        // sigma from radius for gaussian
        /*const*/ float BlurSigma = float(ao_blur_radius) * 0.5;
        /*const*/ float BlurFalloff = 1.0 / (2.0*BlurSigma*BlurSigma);
        /*const*/ float BlurFeatureFalloff = 0.5 * ao_blur_feature*ao_blur_feature;

        for (int i = -ao_blur_radius; i <= ao_blur_radius; ++i) {
            for (int j = -ao_blur_radius; j <= ao_blur_radius; ++j) {
                vec2 uvij = textureCoordinate + vec2(i,j) * texelSize;
                float w = 1.0;
                float exponent = (i*i+j*j)*BlurFalloff;
                if (ao_blur_feature > 0.00001) {
                    // TODO gl_FragCoord.z ?
                    float dz = texture(positionTexture, uvij).z - position.z;
                    exponent -= dz*dz*BlurFeatureFalloff;
                }
                w = exp(exponent);
                // exp( - (dij)^2/(2 sigd^2) - dz^2/(2 sigr^2) ) 
                blurredOcclusion += w * texture(ambientOcclusionTexture, uvij).x;
                wtotal += w;
            }
        }
        blurredOcclusion /= wtotal;
    } else {
        blurredOcclusion = texture(ambientOcclusionTexture, textureCoordinate).x;
    }

    // ssao power to strengthen effect
    blurredOcclusion = pow(blurredOcclusion, ao_pow); 

#ifndef DRAW_OCCLUSION_FACTOR
    fragmentColor.rgb = albedo*blurredOcclusion;
#else
    fragmentColor.rgb = vec3(blurredOcclusion);
#endif

#endif
}


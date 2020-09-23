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

in highp vec3 viewPosition;
in highp vec3 viewNormal;

in ExtraData {
    highp vec2 uv;
} gs_out;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;

uniform vec4 diffuseColor;
uniform sampler2D matcap;
uniform sampler2D tex_cloth;
uniform float tex_scale = 1;



void main() {
    // color = diffuseColor; // TODO SAMPLE TEXTURE based on uv = viewnormal.xy 
    // color = vec4(texture(matcap, viewNormal.xy).rgb, 1.0);
    // vec2 uv = viewNormal.xy * 2 - 1;
    vec2 mat_uv = viewNormal.xy * 0.5 + 0.5;
    color = vec4(texture(matcap, mat_uv).rgb, 1.0) * diffuseColor;
    position = viewPosition;
    normal = viewNormal; //;vec3(0.0,0.0,1.0);

    color.rgb *= texture(tex_cloth, gs_out.uv * tex_scale).rgb;
}
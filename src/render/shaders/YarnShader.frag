in highp vec3 viewPosition;
in highp vec3 viewNormal;

in G2F {
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
    vec2 mat_uv = viewNormal.xy * 0.5 + 0.5;
    color = vec4(texture(matcap, mat_uv).rgb, 1.0) * diffuseColor;
    position = viewPosition;
    normal = viewNormal; //;vec3(0.0,0.0,1.0);

    color.rgb *= texture(tex_cloth, gs_out.uv * tex_scale).rgb;
}
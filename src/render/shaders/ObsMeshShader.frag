in highp vec3 viewPosition;
in highp vec3 viewNormal;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;

uniform sampler2D matcap;

void main() {
    vec2 mat_uv = viewNormal.xy * 0.5 + 0.5; //?normalize normal? flat or bary?
    color = vec4(texture(matcap, mat_uv).rgb, 1.0);
    position = viewPosition;
    normal = viewNormal;
}
in highp vec3 viewPosition;
in highp vec3 viewNormal;
in highp vec3 B;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;


#ifdef LINES
  uniform vec4 diffuseColor;
#else
  uniform sampler2D matcap;
  uniform float wiresize = 0.015;
#endif

void main() {
  position = viewPosition;
  normal = viewNormal;
#ifdef LINES
  uniform vec4 diffuseColor;
  color = diffuseColor;
#else
  bool backfacing = !gl_FrontFacing;//normal.z < 0;
  if (backfacing)
    normal *= -1;
  vec2 mat_uv = normal.xy * 0.5 + 0.5;
  color = vec4(texture(matcap, mat_uv).rgb, 1.0);
  if (backfacing)
    color.rgb = vec3(length(color.rgb)*0.5);
  if (B.x < wiresize || B.y < wiresize || B.z < wiresize)
    color.rgb *= 0;
#endif
}
in highp vec3 viewPosition;
in highp vec3 viewNormal;
in highp vec3 B;
in highp vec2 uv;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;

// #define LINES
#ifdef LINES
  uniform vec4 diffuseColor = vec4(0.3,0.3,1.0,1.0);
#else
  uniform sampler2D matcap;
  uniform sampler2D tex_cloth;
  uniform float uv_scale = 1;
  uniform vec2 uv_offset = vec2(0);
  // uniform float wiresize = 0.015;
  // uniform float wiresize = 0.03;
  uniform float wiresize = 0.05;
#endif

void main() {
  position = viewPosition;
  normal = viewNormal;
#ifdef LINES
  color = diffuseColor;
#else
  bool backfacing = !gl_FrontFacing; // normal.z < 0;
  if (backfacing)
    normal *= -1;
  // texture
  color = vec4(1.0);
  // color = vec4(texture(tex_cloth, uv * uv_scale * vec2(0.5,1)+ uv_offset).rgb, 1.0); // NOTE: stock scale around 190/130? rib around 120 with offset -0.5x
  // color = vec4(texture(tex_cloth, uv).rgb, 1.0);
  // matcap
  vec2 mat_uv = normal.xy * 0.5 + 0.5;
  color *= vec4(texture(matcap, mat_uv).rgb, 1.0);
  if (backfacing)
    color.rgb = vec3(length(color.rgb)*0.5);
  if (B.x < wiresize || B.y < wiresize || B.z < wiresize)
    color.rgb *= 0;
#endif
}
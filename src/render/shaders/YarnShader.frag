in G2F {
  vec2 plycoord; // ply texture coordinate: N a, l/L
  // float a;
  float r;
  vec3 p;
  mat3 Q; // tangential surface frame T1,T2,N
  vec2 uv;
} gs_out;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;

uniform sampler2D matcap;
uniform sampler2D tex_cloth;
// uniform sampler1D normalMap;
uniform sampler2D normalMap;
uniform float uv_scale = 1;
uniform vec2 uv_offset = vec2(0);

uniform float radius;
uniform float plyTwist;
uniform float plyNum;
uniform float plyHeight; // multiple  of radius
uniform float plyLen; // multiple  of radius/Num
uniform float aostrength = 1.0;//0.6; // 0.85

void main() {

  // ply/fiber twisted normal map / ambient occlusion

  vec3 tx = texture(normalMap, gs_out.plycoord).rgb; // nx,ny,ao
  tx.rg = (tx.rg - 0.5) * 2;
  float nz = sqrt(1 - tx.r*tx.r - tx.g*tx.g);
  vec3 nflat = vec3(
    0.31830988*plyHeight*plyNum*tx.r,
    2*plyHeight*plyNum/plyLen * (tx.g - plyTwist*tx.r),
    nz);

  // NOTE not orthonormalizing Q after interpolation, but normalizing n afterwards
  normal = normalize(gs_out.Q * nflat);

  float ao = 1 - aostrength*(1-tx.b);

  // mat cap shading
  vec2 mat_uv = normal.xy * 0.5 + 0.5;
  color = vec4(texture(matcap, mat_uv).rgb, 1.0);
  position = gs_out.p;

  // cloth-scale texture
  color.rgb *= texture(tex_cloth, gs_out.uv * uv_scale + uv_offset).rgb;

  // ply/fiber ambient occlusion
  color.rgb *= ao;
  // color.rgb = normal;
}
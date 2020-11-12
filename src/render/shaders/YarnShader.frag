in G2F {
  float a;
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

uniform vec4 diffuseColor;
uniform sampler2D matcap;
uniform sampler2D tex_cloth;
uniform sampler1D normalMap;
uniform float tex_scale = 1;

uniform float radius;
uniform float normalTwist;
uniform float normalNum;
uniform float normalHeight = 0.1; // multiple  of radius

void main() {

  float ao;
  if (false) {
  // normal = gs_out.n; //;vec3(0.0,0.0,1.0);
    normal = normalize(gs_out.Q[2]);
    ao = 1;
  }
  else {
    vec3 tx = texture(normalMap, gs_out.a).rgb;
    float nx = (tx.r - 0.5) * 2;
    float H = normalHeight * gs_out.r;
    vec3 nflat = vec3(H*normalNum*nx, -normalTwist * H * normalNum * nx, 6.28318530 * gs_out.r);
    // note not normalizing Q after interpolation, but normalizing n afterwards
    normal = normalize(gs_out.Q * nflat);

    float aostrength = 0.5;
    ao = 1 - aostrength*(1-tx.b);
  }

  // TODO multiply with AO tx.z



  vec2 mat_uv = normal.xy * 0.5 + 0.5;
  color = vec4(texture(matcap, mat_uv).rgb, 1.0) * diffuseColor;
  position = gs_out.p;

  color.rgb *= texture(tex_cloth, gs_out.uv * tex_scale).rgb;


  // float tx = pow(abs(sin(gs_out.a)),0.5);
  // color.rgb *= 0.4 + 0.6*tx;
  // if (mod(abs(gs_out.a),1.0)<0.1) {
  //     color.rgb *= 0.4;//vec3(1,0,0);
  // }

  // sample val = sampler(map, N a - N V l / 2 pi r)
  // color.rgb = tx;

  color.rgb *= ao;
  // color.rgb = normal;
}
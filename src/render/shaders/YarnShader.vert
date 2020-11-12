in highp vec4 position;
in float arc;
in highp vec3 d1;
in highp vec2 uv;
in float local_radius;

uniform mat4 transformation;
uniform mat3 normalMatrix;

out V2G {
  vec3 d1;
  float arc;
  float th;
  vec2 uv;
  float r;
} vs_out;

void main(){
  gl_Position = transformation * vec4(position.rgb,1); // pos.w == twist, shoudlnt use for vec4 trafos
  // gl_Position = transformation * vec4(uv,0,1);
  vs_out.d1 = normalMatrix * d1;
  vs_out.arc = arc;
  vs_out.uv = uv;
  vs_out.r = local_radius;
  vs_out.th = position.a;
  // vs_out.viewPosition_vert = transformation*position;
}
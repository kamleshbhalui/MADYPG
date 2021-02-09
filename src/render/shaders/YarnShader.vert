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
  // float r; // deprecate
  // float pad0;
} vs_out;

void main(){
  // NOTE: pos.w or pos.a is the edge twist.
  gl_Position = transformation * vec4(position.rgb,1);
  vs_out.d1 = normalMatrix * d1;
  vs_out.arc = arc;
  vs_out.uv = uv;
  // vs_out.r = local_radius; // deprecate
  vs_out.th = position.a;
}
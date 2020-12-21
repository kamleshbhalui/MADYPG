in highp vec4 position;

out V2G {
  vec2 uv;
} vs_out;

uniform mat4 transformation;
uniform float scale;
uniform float dY;

void main() {
  gl_Position = vec4(position.x*scale,position.y+dY,position.z*scale,1.0);
  vs_out.uv = vec2(gl_Position.x, -gl_Position.z);
  gl_Position = transformation*gl_Position;
}
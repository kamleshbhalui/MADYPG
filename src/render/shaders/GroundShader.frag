in highp vec3 viewPosition;
in highp vec3 viewNormal;
in highp vec2 uv;

layout(location = 0)
out highp vec4 color;
layout(location = 1)
out highp vec3 position;
layout(location = 2)
out highp vec3 normal;

uniform float axw = 0.005;
uniform float gcw = 0.1; // grid cell width
uniform float glw = 0.01; // grid line width
uniform sampler2D tex;

void main() {
  position = viewPosition;
  normal = viewNormal;
  color = vec4(1.0);

  float grid = texture(tex, uv / gcw).r;
  float a = 0.7;
  color.rgb *= 1-a*(1-grid);

  // if (abs(fract(abs(U.x / gcw - 0.5))-0.5) < glw)
  //   color = vec4(0.3,0.3,0.3,1.0);
  // if (abs(fract(abs(U.y / gcw - 0.5))-0.5) < glw)
  //   color = vec4(0.3,0.3,0.3,1.0);

  if (abs(uv.x) < axw)
    color = vec4(0.4,1.0,0.6,1.0);
  if (abs(uv.y) < axw)
    color = vec4(1.0,0.5,0.5,1.0);
}
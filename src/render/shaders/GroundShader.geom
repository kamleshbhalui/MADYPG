layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

out highp vec3 viewPosition;
out highp vec2 uv;
out highp vec3 viewNormal;

in V2G {
  vec2 uv;
} gs_in[];

uniform mat4 projection;

void main() { 
  viewNormal = normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz));
  for(int i=0; i<3; i++) { // triangles
    uv = gs_in[i].uv;
    viewPosition = gl_in[i].gl_Position.xyz;
    gl_Position = projection*vec4(viewPosition, 1.0);
    EmitVertex();
  }
  EndPrimitive();
}

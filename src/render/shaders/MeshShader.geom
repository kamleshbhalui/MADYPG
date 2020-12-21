layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;
// #define Nplus1 4
// layout (line_strip, max_vertices = Nplus1) out;

out highp vec3 viewPosition;
out highp vec3 viewNormal;
out highp vec3 B;

uniform mat4 projection;

void main() { 
  viewNormal = normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz));
  // for(int i=0; i<Nplus1; i++) { // lines
    // viewPosition = gl_in[i % (Nplus1-1)].gl_Position.xyz;
  for(int i=0; i<3; i++) { // triangles
    viewPosition = gl_in[i].gl_Position.xyz;
    gl_Position = projection*vec4(viewPosition, 1.0);
    B = vec3(0);
    B[i] = 1;
    EmitVertex();
  }
  EndPrimitive();
}

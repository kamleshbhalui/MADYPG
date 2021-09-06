layout (triangles) in;
// #define LINES
#ifdef LINES
#define Nplus1 4
layout (line_strip, max_vertices = Nplus1) out;
#else
layout (triangle_strip, max_vertices = 3) out;
#endif

in int gl_PrimitiveIDIn;

// assume same order of verts in world-face and mat-face
// then: face index -> uv = U[face[tri], i]
layout(std430, binding = 0) readonly buffer _meshFmsBuffer
{
  uint[3] data[];
} buf_meshFms;
layout(std430, binding = 1) readonly buffer _meshUBuffer
{
  vec2 data[];
} buf_meshU;

out highp vec3 viewPosition;
out highp vec3 viewNormal;
out highp vec3 B;
out highp vec2 uv;

uniform mat4 projection;

void main() { 
  viewNormal = normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz));
  #ifdef LINES
  for(int i=0; i<Nplus1; i++) { // lines
    viewPosition = gl_in[i % (Nplus1-1)].gl_Position.xyz;
  #else
  uint tri = gl_PrimitiveIDIn;
  uint[3] ms_ixs = buf_meshFms.data[tri];
  for(int i=0; i<3; i++) { // triangles
    viewPosition = gl_in[i].gl_Position.xyz;
    uv = buf_meshU.data[ms_ixs[i]];
  #endif
    gl_Position = projection*vec4(viewPosition, 1.0);
    B = vec3(0);
    B[i] = 1;
    EmitVertex();
  }
  EndPrimitive();
}

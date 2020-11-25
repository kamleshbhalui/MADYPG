#define NVERTICES 16 // NOTE: for cylinder segment, including both ends

#define Vprv gl_in[0]
#define V0 gl_in[1]
#define V1 gl_in[2]
#define Vnxt gl_in[3]
layout (lines_adjacency) in;
// layout (lines) in;
layout (triangle_strip, max_vertices = NVERTICES) out;

in V2G {
  vec3 d1;
  float arc;
  float th;
  vec2 uv;
  float r; // this rlocal will be deprecated once we compute rVolPreserve in here from arc! 
  // float pad0;
} gs_in[];

out G2F {
  float a; // radial texture coordinate: actually Na - NVl/2 pi r
  float r;
  vec3 p;
  // vec3 n;
  mat3 Q; // tangential surface frame T1,T2,N, where T1=t x n, T2=t, N=n. for yarn curve tangent n and normal n. 
  vec2 uv;
} gs_out;

uniform float radius = 1.0;
uniform float normalTwist = 1000;
uniform float normalNum = 4;
uniform mat4 projection;


void main() { 
  // note: variable naming: vertex values indexed as 0 1 2 3, edge values A B C.
  // [v0] --eA-- [v1] --eB-- [v2] --eC-- [v3]
  // shader outputs the geometry for segment eB
  
  // edge rest lengths and inverse for weighted averaging
  float rlA = gs_in[1].arc - gs_in[0].arc;
  float rlB = gs_in[2].arc - gs_in[1].arc;
  float rlC = gs_in[3].arc - gs_in[2].arc;
  float invrlAB = 1.0 / (rlA + rlB);
  float invrlBC = 1.0 / (rlB + rlC);

  #define volpreserve
  #ifdef volpreserve
  const float minR = 0.25, maxR = 2.0;
  // r0^2 l0 pi == rt^2 lt pi --> rt = r0 sqrt(l0/lt)
  vec3 dA =  gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
  vec3 dB =  gl_in[2].gl_Position.xyz - gl_in[1].gl_Position.xyz;
  vec3 dC =  gl_in[3].gl_Position.xyz - gl_in[2].gl_Position.xyz;
  float rA = min(max(minR, sqrt(rlA) * pow(dot(dA,dA),-0.25)), maxR);
  float rB = min(max(minR, sqrt(rlB) * pow(dot(dB,dB),-0.25)), maxR);
  float rC = min(max(minR, sqrt(rlC) * pow(dot(dC,dC),-0.25)), maxR);
  float r1 = (rlA * rA + rlB * rB) * invrlAB * radius;
  float r2 = (rlB * rB + rlC * rC) * invrlBC * radius;
  #else
  float r1 = radius; // gs_in[1].r *
  float r2 = radius; // gs_in[2].r *
  #endif

  // segment tangent, and weighted avg vertex tangents
  vec3 t = V1.gl_Position.xyz - V0.gl_Position.xyz;
  vec3 tv0 = normalize(rlB * t + rlA * (V0.gl_Position.xyz - Vprv.gl_Position.xyz));
  vec3 tv1 = normalize(rlB * t + rlC * (Vnxt.gl_Position.xyz - V1.gl_Position.xyz));

  vec3 nv0 = (rlB * gs_in[1].d1 + rlA * gs_in[0].d1);
  vec3 nv1 = (rlC * gs_in[2].d1 + rlB * gs_in[1].d1);
  // NOTE: when using twists (and should!) probably should compute material frame BEFORE averaging since (m1+m2) != R(thetaavg) (d1+d2)
  // vec3 d1_0 = cos(gs_in[0].th) * gs_in[0].d1 + sin(gs_in[0].th) * cross(normalize(V0.gl_Position.xyz - Vprv.gl_Position.xyz),gs_in[0].d1);
  // vec3 d1_1 = cos(gs_in[1].th) * gs_in[1].d1 + sin(gs_in[1].th) * cross(normalize(t),gs_in[1].d1);
  // vec3 d1_2 = cos(gs_in[2].th) * gs_in[2].d1 + sin(gs_in[2].th) * cross(normalize(Vnxt.gl_Position.xyz - V1.gl_Position.xyz),gs_in[2].d1);
  // vec3 nv0 = (d1_1+d1_0);
  // vec3 nv1 = (d1_2+d1_1);

  // normals and binormals
  nv0 = normalize(nv0 - tv0 * dot(tv0,nv0));///dot(tv0,tv0));
  vec3 bv0 = normalize(cross( tv0, nv0 ));
  nv1 = normalize(nv1 - tv1 * dot(tv1,nv1));///dot(tv1,tv1));
  vec3 bv1 = normalize(cross( tv1, nv1 ));

  #define usetheta
  #ifdef usetheta
  // NOTE actually instead of computing ca1 ca2 sa1 sa2 for each circle vertex,
  // could instead just apply the twist beforehand to each edge segment, and
  // then average the material frames. this might cut down on a bunch of
  // sin/cos, and should work as long as consecutive material frames are
  // less than +180deg apart
  float th1 = (rlB * gs_in[1].th + rlA * gs_in[0].th) * invrlAB;
  float th2 = (rlC * gs_in[2].th + rlB * gs_in[1].th) * invrlBC;
  #endif

  int segs = NVERTICES/2;
  vec3 n;
  for(int i=0; i<segs; i++) {
    float alpha = i/float(segs-1);
    float a = alpha * 2.0 * 3.14159;
    #ifdef usetheta
    float ca1 = cos(a + th1); float sa1 = sin(a + th1);
    float ca2 = cos(a + th2); float sa2 = sin(a + th2); // NOTE: could speed up by ignoring twist theta
    #else
    float ca1 = cos(a); float sa1 = sin(a);
    float ca2 = ca1; float sa2 = sa1; // NOTE: could speed up by ignoring twist theta
    #endif
    // float ca = cos(a + gs_in[1].th); float sa = sin(a + gs_in[1].th);

    // viewNormal = vec3( ca*b.x + sa*n.x,
    //                 ca*b.y + sa*n.y,
    //                 ca*b.z + sa*n.z );
    
    gs_out.r = radius;
    // gs_out.r = r2; // doesnt work well with texture
    // NOTE: -a == (1-a) because of expecting lookup orientation in normal map
    gs_out.a = normalNum*(-alpha - 0.1591549 * normalTwist * gs_in[2].arc / gs_out.r);
    gs_out.uv = gs_in[2].uv;
    // gs_out.uv = vec2(twistSpeed * gs_in[2].arc,0);
    n = ca2*nv1 + sa2 * bv1;
    gs_out.Q = mat3(cross(tv1, n), tv1, n);
    gs_out.p = V1.gl_Position.xyz + r2 * n;
    gl_Position = projection*vec4(gs_out.p, 1.0);
    EmitVertex();   
    
    gs_out.r = radius;
    // gs_out.r = r1; // doesnt work well with texture
    gs_out.a = normalNum*(-alpha - 0.1591549 * normalTwist * gs_in[1].arc / gs_out.r);
    gs_out.uv = gs_in[1].uv;
    // gs_out.uv = vec2(twistSpeed * gs_in[1].arc,0);
    n = ca1*nv0 + sa1 * bv0;
    gs_out.Q = mat3(cross(tv0, n), tv0, n);
    gs_out.p = V0.gl_Position.xyz + r1 * n;  
    gl_Position = projection*vec4(gs_out.p, 1.0);
    EmitVertex();  
  }
  EndPrimitive();  
}

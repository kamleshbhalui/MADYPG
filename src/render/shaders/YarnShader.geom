
#define randomize // randomize cylinder vertices to fake fuzz
#define volpreserve // scale radius to preserve volume
#define usetheta // use twist kinematics to twist the cylinders and textures
// #define DRAW_LINES // if this is active, draws simple lines instead of cylinders, incorrect lighting due to constant normal
#define NVERTICES 16 // NOTE: for cylinder segment, including both ends, ie circular crosssection will have nverts/2 vertices
#define LEVEL_OF_DETAIL // reduce min. number of cylinder vertices for far away geometry


#define Vprv gl_in[0]
#define V0 gl_in[1]
#define V1 gl_in[2]
#define Vnxt gl_in[3]
layout (lines_adjacency) in;
// layout (lines) in;

#ifdef DRAW_LINES
layout (line_strip, max_vertices = 2) out;
#else
layout (triangle_strip, max_vertices = NVERTICES) out;
#endif

in V2G {
  vec3 d1;
  float arc;
  float th;
  vec2 uv;
  float r; // this rlocal will be deprecated once we compute rVolPreserve in here from arc! 
  // float pad0;
} gs_in[];

out G2F {
  vec2 plycoord; // ply texture coordinate: N a, l/L
  // float a; // radial texture coordinate: actually Na - NVl/2 pi r
  float r;
  vec3 p;
  // vec3 n;
  mat3 Q; // tangential surface frame T1,T2,N, where T1=t x n, T2=t, N=n. for yarn curve tangent n and normal n. 
  vec2 uv;
} gs_out;

uniform float radius = 1.0;
uniform float plyTwist = 1000;
uniform float plyLen = 1; // multiple  of radius
uniform float plyNum = 4;
uniform mat4 projection;

float rand(int x, int y){
    return fract(sin(x*12.9898 + y*78.233) * 43758.5453);
}

// hard coded level of detail, from 1 at 10cm distance to 0 at 1m distance
const float lod_z_near = -0.1;
const float lod_z_far = -1.0;
const float lod_z_inv = 1.0/(lod_z_far-lod_z_near);
const int NSEGSCLOSE=NVERTICES/2;
const int NSEGSSFAR=4;
float LOD(float z){
  return 1 - clamp((z-lod_z_near)*lod_z_inv,0,1);
}

void main() { 
#ifdef DRAW_LINES
  gs_out.plycoord = vec2(0,0);
  gs_out.r = radius;
  vec3 t = normalize(V1.gl_Position.xyz - V0.gl_Position.xyz);
  vec3 n = vec3(0,0,1) - t * t.z; // n - n.t t
  float n2 = dot(n,n);
  if (n2 < 0.0001)
    n = normalize(vec3(1,0,0) - t * t.x);
  else
    n /= sqrt(n2);
  gs_out.Q = mat3(cross(t, n), t, n);
  // gs_out.Q = mat3(vec3(1,0,0),vec3(0,1,0),vec3(0,0,1));
  gs_out.uv = gs_in[1].uv;
  gs_out.p = gl_in[1].gl_Position.xyz;
  gl_Position = projection*vec4(gs_out.p, 1.0);
  EmitVertex();
  gs_out.uv = gs_in[2].uv;
  gs_out.p = gl_in[2].gl_Position.xyz;
  gl_Position = projection*vec4(gs_out.p, 1.0);
  EmitVertex();
  EndPrimitive();
  return;
#else


  // note: variable naming: vertex values indexed as 0 1 2 3, edge values A B C.
  // [v0] --eA-- [v1] --eB-- [v2] --eC-- [v3]
  // shader outputs the geometry for segment eB
  
  // edge rest lengths and inverse for weighted averaging
  float rlA = gs_in[1].arc - gs_in[0].arc;
  float rlB = gs_in[2].arc - gs_in[1].arc;
  float rlC = gs_in[3].arc - gs_in[2].arc;
  float invrlAB = 1.0 / (rlA + rlB);
  float invrlBC = 1.0 / (rlB + rlC);

  // rlA=1;
  // rlB=1;
  // rlC=1;
  // invrlAB=0.5;
  // invrlBC=0.5;

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

  // normals and binormals
  nv0 = normalize(nv0 - tv0 * dot(tv0,nv0));///dot(tv0,tv0));
  vec3 bv0 = normalize(cross( tv0, nv0 ));
  nv1 = normalize(nv1 - tv1 * dot(tv1,nv1));///dot(tv1,tv1));
  vec3 bv1 = normalize(cross( tv1, nv1 ));
  

  #ifdef usetheta
  // NOTE actually instead of computing ca1 ca2 sa1 sa2 for each circle vertex,
  // could instead just apply the twist beforehand to each edge segment, and
  // then average the material frames. this might cut down on a bunch of
  // sin/cos, and should work as long as consecutive material frames are
  // less than +180deg apart
  float th1 = (rlB * gs_in[1].th + rlA * gs_in[0].th) * invrlAB;
  float th2 = (rlC * gs_in[2].th + rlB * gs_in[1].th) * invrlBC;
  #endif

  #ifdef LEVEL_OF_DETAIL
  float lod = LOD(0.5*(V1.gl_Position.z+V0.gl_Position.z));
  int segs = int(lod*NSEGSCLOSE+(1-lod)*NSEGSSFAR);
  #else
  int segs = NVERTICES/2;
  #endif
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
    // gs_out.a = plyNum*(-alpha - 0.1591549 * plyTwist * gs_in[2].arc / gs_out.r);
    gs_out.plycoord = vec2(plyNum*-alpha - gs_in[2].arc/(plyLen*radius)*plyTwist*plyNum , gs_in[2].arc / (plyLen*radius)); // TODO precomp invradius and l/L
    // gs_out.plycoord = vec2(plyNum*-alpha - gs_in[2].arc/(plyLen*radius)*plyTwist , gs_in[2].arc / (plyLen*radius)); // TODO precomp invradius and l/L
    gs_out.uv = gs_in[2].uv;
    // gs_out.uv = vec2(segs/8.0,0);
    // gs_out.uv = vec2(twistSpeed * gs_in[2].arc,0);
    n = ca2*nv1 + sa2 * bv1;
    gs_out.Q = mat3(cross(tv1, n), tv1, n);
    // noise? xyz + (r+rand(vertex,i)) * n;
    #ifdef randomize
    // gs_out.p = V1.gl_Position.xyz +  (0.5+1.5*rand(gl_PrimitiveIDIn+1,i)) * r2 * n;
    gs_out.p = V1.gl_Position.xyz +  (0.5+1.5*rand(gl_PrimitiveIDIn,2*i)) * r2 * n;
    #else
    gs_out.p = V1.gl_Position.xyz + r2 * n;
    #endif
    gl_Position = projection*vec4(gs_out.p, 1.0);
    EmitVertex();   
    
    gs_out.r = radius;
    // gs_out.r = r1; // doesnt work well with texture
    // gs_out.a = plyNum*(-alpha - 0.1591549 * plyTwist * gs_in[1].arc / gs_out.r);
    gs_out.plycoord = vec2(plyNum*-alpha - gs_in[1].arc/(plyLen*radius)*plyTwist*plyNum , gs_in[1].arc / (plyLen*radius)); // TODO precomp invradius and l/L
    // gs_out.plycoord = vec2(plyNum*-alpha - gs_in[1].arc/(plyLen*radius)*plyTwist , gs_in[1].arc / (plyLen*radius)); // TODO precomp invradius and l/L
    gs_out.uv = gs_in[1].uv;
    // gs_out.uv = vec2(segs/8.0,0);
    // gs_out.uv = vec2(twistSpeed * gs_in[1].arc,0);
    n = ca1*nv0 + sa1 * bv0;
    gs_out.Q = mat3(cross(tv0, n), tv0, n);
    #ifdef randomize
    // gs_out.p = V0.gl_Position.xyz +  (0.5+1.5*rand(gl_PrimitiveIDIn+0,i)) * r1 * n;
    gs_out.p = V0.gl_Position.xyz +  (0.5+1.5*rand(gl_PrimitiveIDIn,2*i+1)) * r1 * n;
    #else
    gs_out.p = V0.gl_Position.xyz + r1 * n;
    #endif
    gl_Position = projection*vec4(gs_out.p, 1.0);
    EmitVertex();  
  }
  EndPrimitive();  

#endif // DRAW_LINES

}

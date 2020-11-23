#define NVERTICES 16

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

// https://github.com/torbjoern/polydraw_scripts/blob/master/geometry/drawcone_geoshader.pss
vec3 createPerp(vec3 invec)
{
  vec3 ret = cross( invec, vec3(0.0, 0.0, 1.0) );
  if ( length(ret) == 0.0 )
  {
    ret = cross( invec, vec3(0.0, 1.0, 0.0) );
  }
  return ret;
}


// void main() { 

//    float r1 = 1.0 * radius;
//    float r2 = 1.0 * radius;

//    vec3 axis = V1.gl_Position.xyz - V0.gl_Position.xyz;

//    vec3 perpx = normalize(createPerp( axis ));
//    vec3 perpy = normalize(cross( axis, perpx ));
//    int segs = NVERTICES/2;
//    for(int i=0; i<segs; i++) {
//       float a = i/float(segs-1) * 2.0 * 3.14159;
//       float ca = cos(a); float sa = sin(a);

//       viewNormal = vec3( ca*perpx.x + sa*perpy.x,
//                      ca*perpx.y + sa*perpy.y,
//                      ca*perpx.z + sa*perpy.z );


//       vec3 p1 = V0.gl_Position.xyz + r1*viewNormal;
//       vec3 p2 = V1.gl_Position.xyz + r2*viewNormal;

//       viewPosition = p2;
//       gl_Position = projection*vec4(p2, 1.0); EmitVertex();    
//       viewPosition = p1;
//       gl_Position = projection*vec4(p1, 1.0); EmitVertex();   
//    }
//    EndPrimitive();  
// }



void main() { 
  // edge rest lengths and inverse for weighted averaging
  float rlA = gs_in[1].arc - gs_in[0].arc;
  float rlB = gs_in[2].arc - gs_in[1].arc;
  float rlC = gs_in[3].arc - gs_in[2].arc;
  float invrlAB = 1.0 / (rlA + rlB);
  float invrlBC = 1.0 / (rlB + rlC);


  // current dists
  vec3 dA =  gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
  vec3 dB =  gl_in[2].gl_Position.xyz - gl_in[1].gl_Position.xyz;
  vec3 dC =  gl_in[3].gl_Position.xyz - gl_in[2].gl_Position.xyz;

  #define volpreserve
  #ifdef volpreserve
  const float minR = 0.25, maxR = 2.0;
  // r0^2 l0 pi == rt^2 lt pi --> rt = r0 sqrt(l0/lt)
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
  // TODO weigh by arclength!
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
  nv0 = normalize(nv0 - tv0 * dot(tv0,nv0)/dot(tv0,tv0));
  // vec3 nv0 = normalize(createPerp( tv0 ));
  vec3 bv0 = normalize(cross( tv0, nv0 ));
  nv1 = normalize(nv1 - tv1 * dot(tv1,nv1)/dot(tv1,tv1));
  // vec3 nv1 = normalize(createPerp( tv1 ));
  vec3 bv1 = normalize(cross( tv1, nv1 ));

  // vec3 n = normalize(createPerp( t ));
  // vec3 b = normalize(cross( t, n ));

  #define usetheta
  #ifdef usetheta
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

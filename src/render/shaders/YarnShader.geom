#define NVERTICES 16

#define Vprv gl_in[0]
#define V0 gl_in[1]
#define V1 gl_in[2]
#define Vnxt gl_in[3]
layout (lines_adjacency) in;
// layout (lines) in;
layout (triangle_strip, max_vertices = NVERTICES) out;

// out VS_OUT {
//     highp vec4 viewPosition_vert;
// } gs_in[];

out highp vec3 viewPosition;
out highp vec3 viewNormal;

uniform float radius = 1.0;
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

  float r1 = 1.0 * radius;
  float r2 = 1.0 * radius;

  // segment tangent, and weighted avg vertex tangents
  vec3 t = V1.gl_Position.xyz - V0.gl_Position.xyz;
  vec3 tv0 = (t + (V0.gl_Position.xyz - Vprv.gl_Position.xyz));
  vec3 tv1 = (t + (Vnxt.gl_Position.xyz - V1.gl_Position.xyz));

  // normals and binormals
  vec3 nv0 = normalize(createPerp( tv0 ));
  vec3 bv0 = normalize(cross( tv0, nv0 ));
  vec3 nv1 = normalize(createPerp( tv1 ));
  vec3 bv1 = normalize(cross( tv1, nv1 ));

  // vec3 n = normalize(createPerp( t ));
  // vec3 b = normalize(cross( t, n ));

  int segs = NVERTICES/2;
  for(int i=0; i<segs; i++) {
    float a = i/float(segs-1) * 2.0 * 3.14159;
    float ca = cos(a); float sa = sin(a);

    // viewNormal = vec3( ca*b.x + sa*n.x,
    //                 ca*b.y + sa*n.y,
    //                 ca*b.z + sa*n.z );

    viewNormal = vec3( ca*nv1.x + sa*bv1.x,
                ca*nv1.y + sa*bv1.y,
                ca*nv1.z + sa*bv1.z );
    viewPosition = V1.gl_Position.xyz + r2*viewNormal;
    gl_Position = projection*vec4(viewPosition, 1.0); EmitVertex();   
    
    viewNormal = vec3( ca*nv0.x + sa*bv0.x,
                ca*nv0.y + sa*bv0.y,
                ca*nv0.z + sa*bv0.z );
    viewPosition = V0.gl_Position.xyz + r1*viewNormal;  
    gl_Position = projection*vec4(viewPosition, 1.0); EmitVertex();  

  }
  EndPrimitive();  
}

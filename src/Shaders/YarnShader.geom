#define NVERTICES 16
layout (lines) in;
layout (triangle_strip, max_vertices = NVERTICES) out;

// out VS_OUT {
//     highp vec4 viewPosition_vert;
// } gs_in[];

out highp vec3 viewPosition;
out highp vec3 viewNormal;

uniform float radius = 1.0;
uniform mat4 projection;




// void main() {   
//     vec3 e0 = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
//     vec3 e1 = (gl_in[0].gl_Position.xyz + gl_in[1].gl_Position.xyz) * 0.5 + vec3(0.0, 50.0, 0.0) - gl_in[0].gl_Position.xyz;

//     // viewNormal = vec3(0.0,0.0,1.0); 
//     viewNormal = normalize(cross(e0, e1));//vec3(0.0,0.0,1.0); 

//     gl_Position = gl_in[0].gl_Position ;
//     EmitVertex();
//     gl_Position = gl_in[1].gl_Position ;//+ vec4( 10, 0.0, 0.0, 0.0);
//     EmitVertex();
//     gl_Position = (gl_in[0].gl_Position + gl_in[1].gl_Position) * 0.5 + vec4(0.0, 10.0, 0.0, 0.0); 
//     EmitVertex();
    
//     EndPrimitive();

// }    

// void main() { 

//     // input: gl_Position ... viewspace coordinates
//     vec3 e0 = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
//     vec3 e1 = (gl_in[0].gl_Position.xyz + gl_in[1].gl_Position.xyz) * 0.5 + vec3(0.0, 50.0, 0.0) - gl_in[0].gl_Position.xyz;
//     viewNormal = normalize(cross(e0, e1));

//     // viewNormal = vec3(0.0,0.0,1.0); 

//     viewPosition = gl_in[0].gl_Position.xyz;
//     gl_Position = projection * vec4(viewPosition,1.0);
//     EmitVertex();
//     viewPosition = gl_in[1].gl_Position.xyz;
//     gl_Position = projection * vec4(viewPosition,1.0);
//     EmitVertex();
//     viewPosition = (gl_in[0].gl_Position.xyz + gl_in[1].gl_Position.xyz) * 0.5 + vec3(0.0, 10.0, 0.0);
//     gl_Position = projection * vec4(viewPosition,1.0);
//     EmitVertex();
    
//     EndPrimitive();
// }


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

void main() { 

    // // input: gl_Position ... viewspace coordinates
    // vec3 e0 = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
    // vec3 e1 = (gl_in[0].gl_Position.xyz + gl_in[1].gl_Position.xyz) * 0.5 + vec3(0.0, 50.0, 0.0) - gl_in[0].gl_Position.xyz;
    // viewNormal = normalize(cross(e0, e1));

    // // viewNormal = vec3(0.0,0.0,1.0); 

    // viewPosition = gl_in[0].gl_Position.xyz;
    // gl_Position = projection * vec4(viewPosition,1.0);
    // EmitVertex();
    // viewPosition = gl_in[1].gl_Position.xyz;
    // gl_Position = projection * vec4(viewPosition,1.0);
    // EmitVertex();
    // viewPosition = (gl_in[0].gl_Position.xyz + gl_in[1].gl_Position.xyz) * 0.5 + vec3(0.0, 10.0, 0.0);
    // gl_Position = projection * vec4(viewPosition,1.0);
    // EmitVertex();
    
    // EndPrimitive();


   float r1 = 1.0 * radius;
   float r2 = 1.0 * radius;

   vec3 axis = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;

   vec3 perpx = normalize(createPerp( axis ));
   vec3 perpy = normalize(cross( axis, perpx ));
   int segs = NVERTICES/2;
   for(int i=0; i<segs; i++) {
      float a = i/float(segs-1) * 2.0 * 3.14159;
      float ca = cos(a); float sa = sin(a);

      viewNormal = vec3( ca*perpx.x + sa*perpy.x,
                     ca*perpx.y + sa*perpy.y,
                     ca*perpx.z + sa*perpy.z );


      vec3 p1 = gl_in[0].gl_Position.xyz + r1*viewNormal;
      vec3 p2 = gl_in[1].gl_Position.xyz + r2*viewNormal;
      
      viewPosition = p2;
      gl_Position = projection*vec4(p2, 1.0); EmitVertex();    
      viewPosition = p1;
      gl_Position = projection*vec4(p1, 1.0); EmitVertex();   
   }
   EndPrimitive();  
}

// in highp vec4 position;

// uniform mat4 transformation;
// uniform mat4 projection;

// out highp vec3 viewPosition;

// void main(){
//    vec4 position4 = transformation*position;
//    viewPosition = position4.xyz;
//    gl_Position = projection*position4;
// }
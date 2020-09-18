in highp vec4 position;
uniform mat4 transformation;
uniform float dz = 0;

void main(){
    gl_Position = transformation * (position + vec4(0,0,dz,0));
}
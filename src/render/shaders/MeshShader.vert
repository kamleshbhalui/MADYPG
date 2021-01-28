in highp vec4 position;
uniform mat4 transformation;

void main(){
    gl_Position = transformation * position;
    // gl_Position = transformation * (position+vec4(0,-0.005,0,0));
    // gl_Position = transformation * (position+vec4(0,0,-0.005,0));
}
in highp vec4 position;
uniform mat4 transformation;

void main(){
    gl_Position = transformation * position;
}
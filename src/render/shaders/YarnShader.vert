in highp vec4 position;
in highp vec2 uv;
in float local_radius;

uniform mat4 transformation;

out V2G {
    highp vec2 uv;
    float r;
} vs_out;

void main(){
    gl_Position = transformation * position;
    // gl_Position = transformation * vec4(uv,0,1);
    vs_out.uv = uv;
    vs_out.r = local_radius;
    // vs_out.viewPosition_vert = transformation*position;
}
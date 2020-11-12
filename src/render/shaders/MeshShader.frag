in highp vec3 viewPosition;
in highp vec3 viewNormal;

layout(location = 0)
out highp vec4 color;
// layout(location = 1)
// out highp vec3 position;
// layout(location = 2)
// out highp vec3 normal;

uniform vec4 diffuseColor;
// uniform sampler2D matcap;

void main() {
    color = diffuseColor;// vec4(0.6,0.8,0.95,1.0);
    // position = viewPosition;
    // normal = viewNormal;
}
#define DEFAULT_WRKGRPSIZE 8
// simple test on rib sweater & yarn model
// 4: 70fps
// 8: 74fps (best)
// 16: 66fps
// 32: 50fps
// 64: 48fps
// haven't tested using different values for the shaders
// also note that this might strongly depend on the hardware used
#define SHELLMAPSHADER_WRKGRPSIZE DEFAULT_WRKGRPSIZE
#define DEFORMSHADER_WRKGRPSIZE DEFAULT_WRKGRPSIZE
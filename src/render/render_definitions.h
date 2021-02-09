#ifndef __RENDER_DEFINITIONS__H__
#define __RENDER_DEFINITIONS__H__

#define SSAO_SAMPLES 64 // assumed divisible by MSAA

// #define SUPERSAMPLING 1
#define MSAA 4

// #define SUPERSAMPLING 2 
// #define MSAA 2

// #define SUPERSAMPLING 1
// #define MSAA 16

// NOTE: most figures in the paper have been generated with SSAAx2 MSAAx1
// and most videos with SSAAx2 MSAAx2 or SSAAx1 MSAAx16

// NOTE:
// currently supersampling is at best wonky, might need to more carefully
// implement it, or do some better/complicated AA technique that work at extreme
// scales

#endif // __RENDER_DEFINITIONS__H__
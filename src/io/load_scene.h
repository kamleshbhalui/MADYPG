#ifndef __LOAD_SCENE__H__
#define __LOAD_SCENE__H__

#include <string>

#include "../arcball/ArcBall.h"
#include "../yarns/YarnMapper.h"

void load_scene(int select_scene, std::vector<YarnMapper::Settings>& ymSettings,
                bool& rotate_scene, bool& render_ground, float& ground_height,
                float& ground_scale, std::string& cloth_texture_file,
                Magnum::ArcBall& arcball);
// arcball

#endif  // __LOAD_SCENE__H__

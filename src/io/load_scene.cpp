#include "load_scene.h"

// #include <tinyxml2/tinyxml2.h>
#include "../utils/debug_logging.h"

using namespace Magnum;

void load_scene(int select_scene, std::vector<YarnMapper::Settings>& ymSettings,
                bool& rotate_scene, bool& render_ground, float& ground_height,
                float& ground_scale, std::string& cloth_texture_file,
                Magnum::ArcBall& arcball) {
  ::Debug::logf("Selecting scene setup #%d\n", select_scene);
  switch (select_scene) {
    default:
    case 0:
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/sxsy_const.bin";
      ymSettings[0].objseq_settings.folder   = "data/objseqs/sxsy_const";
      ymSettings[0].objseq_settings.constant_material_space = true;
      ymSettings[0].modelfolder = "data/yarnmodels/model_stock_9";

      rotate_scene       = false;
      render_ground      = false;
      ground_height      = 0;
      ground_scale       = 1;
      cloth_texture_file = "data/textures/colorgridy.jpg";
      break;
    case 1:  // stockinette stretch X 30x30
      rotate_scene = true;
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath =
          "data/binseqs/HYLC_30x30/stock_stretchX.bin";
      ymSettings[0].modelfolder = "data/yarnmodels/model_stock_9";
      cloth_texture_file        = "data/textures/clr_blue.jpg";
      break;
    case 2:  // rib stretch X 30x30
      rotate_scene = true;
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath =
          "data/binseqs/HYLC_30x30/rib_stretchX.bin";
      ymSettings[0].modelfolder = "data/yarnmodels/model_rib_9";
      cloth_texture_file        = "data/textures/clr_red.jpg";
      break;
    case 3:  // sock PBD

      arcball.setViewParameters(2 * Vector3(-0.4f, 0.3f, 0.5f),
                                Vector3(0), Vector3(0, 1, 0));
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::PBD;
      ymSettings[0].modelfolder   = "data/yarnmodels/model_rib_9";
      cloth_texture_file          = "data/textures/clr_red.jpg";
      break;
    case 4:  // cylinder
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/cyldef.bin";
      ymSettings[0].binseq_settings.scale    = 0.125;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_stock_9";
      cloth_texture_file                     = "data/textures/clr_blue.jpg";
      break;
    case 5:  // cylinder
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/cyldef.bin";
      ymSettings[0].binseq_settings.scale    = 0.125;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_rib_9";
      cloth_texture_file                     = "data/textures/clr_red.jpg";
      break;
    case 6:  // cylinder
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/cyldef.bin";
      ymSettings[0].binseq_settings.scale    = 0.125;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_honey_9";
      cloth_texture_file                     = "data/textures/clr_yellow.jpg";
      break;
    case 7:  // cylinder
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/cyldef.bin";
      ymSettings[0].binseq_settings.scale    = 0.125;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_basket_9";
      cloth_texture_file                     = "data/textures/clr_green.jpg";
      break;
    case 8:  // cylinder
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/cyldef.bin";
      ymSettings[0].binseq_settings.scale    = 0.125;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_satin_9";
      cloth_texture_file                     = "data/textures/clr_satin.jpg";
      break;
    case 9:  // shearing
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/shear_grid_x.bin";
      ymSettings[0].modelfolder              = "data/yarnmodels/model_basket_9";
      cloth_texture_file                     = "data/textures/clr_green.jpg";
      break;
    case 10:  // sleeve
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath = "data/binseqs/sleeve.bin";
      ymSettings[0].binseq_settings.scale    = 0.5;
      ymSettings[0].modelfolder              = "data/yarnmodels/model_stock_9";
      cloth_texture_file                     = "data/textures/clr_blue.jpg";
      break;
    case 11:  // sweater stockinette
      arcball.setViewParameters(8 * Vector3(+0.4f, 0.3f, 0.5f), Vector3(0),
                                 Vector3(0, 1, 0));
      render_ground = true;
      ground_height = -0.69;
      ground_scale  = 3;
      rotate_scene  = true;
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath =
          "data/binseqs/HYLC_sweaters/shirt_stock.bin";
      ymSettings[0].modelfolder = "data/yarnmodels/model_stock_9";
      cloth_texture_file        = "data/textures/clr_blue.jpg";
      break;
    case 12:  // sweater rib
      arcball.setViewParameters(8 * Vector3(+0.4f, 0.3f, 0.5f), Vector3(0),
                                 Vector3(0, 1, 0));
      render_ground = true;
      ground_height = -0.69;
      ground_scale  = 3;
      rotate_scene  = true;
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath =
          "data/binseqs/HYLC_sweaters/shirt_rib.bin";
      ymSettings[0].modelfolder = "data/yarnmodels/model_rib_9";
      cloth_texture_file        = "data/textures/clr_red.jpg";
      break;
    case 13:  // sweater tiny stockinette
      arcball.setViewParameters(8 * Vector3(+0.4f, 0.3f, 0.5f), Vector3(0),
                                 Vector3(0, 1, 0));
      render_ground = true;
      ground_height = -0.69;
      ground_scale  = 3;
      rotate_scene  = true;
      ymSettings.emplace_back();
      ymSettings[0].provider_type = YarnMapper::Settings::Provider::BinSeq;
      ymSettings[0].binseq_settings.filepath =
          "data/binseqs/HYLC_sweaters/shirt_stock.bin";
      ymSettings[0].modelfolder = "data/yarnmodels/model_stocktiniest_9";
      cloth_texture_file        = "data/textures/clr_blue.jpg";
      break;
  }
}
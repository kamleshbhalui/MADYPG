#ifndef __MAINAPPLICATION__H__
#define __MAINAPPLICATION__H__

#include <Corrade/Containers/Optional.h>
#include <Corrade/PluginManager/Manager.h>
#include <Corrade/Utility/Arguments.h>
#include <Corrade/Utility/Directory.h>
#include <Corrade/Utility/Resource.h>
#include <Magnum/DebugTools/FrameProfiler.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Framebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/MultisampleTexture.h>
#include <Magnum/GL/Renderbuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Texture.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Primitives/Axis.h>
#include <Magnum/Trade/AbstractImporter.h>
#include <Magnum/Trade/ImageData.h>
#include <Magnum/Trade/MeshData.h>
#include <MagnumPlugins/AnyImageImporter/AnyImageImporter.h>
#include <MagnumPlugins/JpegImporter/JpegImporter.h>
#include <MagnumPlugins/PngImporter/PngImporter.h>
#include <imgui.h>

#include <Magnum/ImGuiIntegration/Context.hpp>
#include <memory>
#include <random>

#include "arcball/ArcBall.h"
#include "arcball/ArcBallCamera.h"
#include "imfilebrowser/imfilebrowser.h"
#include "render/MeshDrawable.h"
#include "render/ObsMeshDrawable.h"
#include "render/YarnDrawable.h"
#include "render/enableMSAA.h"
#include "render/shaders/GroundShader.h"
#include "render/shaders/MeshShader.h"
#include "render/shaders/ObsMeshShader.h"
#include "render/shaders/SsaoApplyShader.h"
#include "render/shaders/SsaoShader.h"
#include "render/shaders/YarnShader.h"
#include "utils/debug_includes.h"
#include "yarns/YarnMapper.h"

#define SUPERSAMPLING 2 // EXTREMELY SLOW AND UNOPTIMIZED EVEN WITH MSAA x1

namespace Magnum {

using namespace Math::Literals;

class MainApplication : public Platform::Application {
 public:
  explicit MainApplication(const Arguments &arguments);

  struct ApplicationSettings {
  } _settings;

 private:
  void drawEvent() override;

  void viewportEvent(ViewportEvent &event) override;

  void keyPressEvent(KeyEvent &event) override;
  void keyReleaseEvent(KeyEvent &event) override;

  void mousePressEvent(MouseEvent &event) override;
  void mouseReleaseEvent(MouseEvent &event) override;
  void mouseMoveEvent(MouseMoveEvent &event) override;
  void mouseScrollEvent(MouseScrollEvent &event) override;
  void textInputEvent(TextInputEvent &event) override;

  void setupFramebuffer(const Vector2i &size);
  void drawSettings();

  void reset_simulation();

  std::vector<YarnDrawable> _yarnDrawable;
  std::unique_ptr<MeshDrawable<MeshShader>> _meshdrawable;
  std::vector<ObsMeshDrawable<ObsMeshShader>> _obsmeshdrawables;
  std::unique_ptr<YarnMapper> _yarnMapper;
  YarnMapper::Settings _yarnMapperSettings;
  std::unique_ptr<ImGui::FileBrowser> _fileDialog;
  std::unique_ptr<ImGui::FileBrowser> _folderDialog;
  GL::Texture2D _matcapObs{NoCreate};
  GL::Texture2D _matcapMesh{NoCreate};
  GL::Texture2D _matcap{NoCreate};
  GL::Texture2D _clothTexture{NoCreate};
  GL::Texture2D _gridTexture{NoCreate};
  // GL::Texture1D _normalMap{NoCreate};
  GL::Texture2D _normalMap{NoCreate};
  bool _gui                      = true;
  bool _paused                   = false;
  int _pauseAt                   = -1;
  int _frame                     = 0;
  int _min_loop_ms               = 16;
  bool _single_step              = false;
  bool _render_mesh              = false;
  bool _render_yarns             = true;
  bool _render_obstacles         = true;
  bool _render_ground            = true;
  bool _rotate_scene             = false;
  std::string _matcap_file       = "data/textures/matcaps/glossyc.jpg";
  std::string _matcapObs_file    = "data/textures/matcaps/lighting1.jpg";
  std::string _matcapMesh_file   = "data/textures/matcaps/blue.jpg";
  std::string _clothtexture_file = "data/textures/colorgridy.jpg";
  std::string _gridtexture_file  = "data/textures/groundgrid.jpg";
  // std::string _normalMap_file = "data/textures/normalMap.jpg";
  std::string _normalMap_file = "data/textures/PLY.jpg";
  float _render_radius_mult   = 1.0f;
  float _render_nmtwist       = 1.0f;
  float _render_nmnum         = 4.0f;
  float _render_nmheight      = 0.4f;  // relative to radius
  float _render_nmlen         = 6.0f;  // relative to radius and 1/num
  float _mesh_dz              = 0.0f;
  float _clothUV_scale   = 1.0f;
  Color4 _bgColor             = Color4(1.0f);  // Color4(Color3(0.2f), 1.0f);

  // std::unique_ptr<ArcBallCamera> _arcballCamera;

  GroundShader _ground{NoCreate};
  MeshShader _meshShader{NoCreate};
  ObsMeshShader _obsMeshShader{NoCreate};
  YarnShader _yarnGeometryShader{NoCreate};
  SsaoApplyShader _ssaoApplyShader{NoCreate};
  SsaoShader _ssaoShader{NoCreate};

  Containers::Optional<ArcBall> _arcball;
  Matrix4 _projection;
  Deg _proj_fov    = 45.0_degf;
  float _proj_near = 0.01f;  // 0.0001f;
  float _proj_far  = 10.0f;  // TODO reduce further?

  GL::Mesh _screenAlignedTriangle{NoCreate};

  GL::Framebuffer _fbo_gbuffer{NoCreate};
  GL::Framebuffer _fbo_ssao{NoCreate};
  GL::Framebuffer _fbo_supersample{NoCreate};

#ifdef MSAA
  GL::MultisampleTexture2D _albedo{NoCreate};
  GL::MultisampleTexture2D _positions{NoCreate};
  GL::MultisampleTexture2D _normals{NoCreate};
  GL::MultisampleTexture2D _depth{NoCreate};
#else
  GL::Texture2D _albedo{NoCreate};
  GL::Texture2D _positions{NoCreate};
  GL::Texture2D _normals{NoCreate};
  GL::Texture2D _depth{NoCreate};
#endif
  GL::Texture2D _occlusion{NoCreate};
  GL::Texture2D _noise{NoCreate};

  /* Profiling */
  DebugTools::GLFrameProfiler _profiler;

  // TODO make a rendersettings struct with default init
  // which can the also be used to reset the settings!
  Color4 _specularColor{0.3};  // TODO remove
#ifdef MSAA
  Float _ao_radius = 0.015f;  // m
  Float _ao_bias   = 0.001f;  // m
  Float _ao_pow    = 1.8f;
#else
  Float _ao_radius = 0.004f;   // m
  Float _ao_bias   = 0.0003f;  // m
  Float _ao_pow    = 2.0f;
#endif
  int _ao_blur_radius    = 0;      // pixels
  float _ao_blur_feature = 25.0f;  // 1/m

  SsaoApplyShader::Flag _ssaoApplyFlag = {};

  ImGuiIntegration::Context _imgui{NoCreate};
};
}  // namespace Magnum

#endif  // __MAINAPPLICATION__H__

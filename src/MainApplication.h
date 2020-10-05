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
#include "render/shaders/MeshShader.h"
#include "render/shaders/ObsMeshShader.h"
#include "render/shaders/SsaoApplyShader.h"
#include "render/shaders/SsaoShader.h"
#include "render/shaders/YarnShader.h"
#include "utils/debug_includes.h"
#include "yarns/YarnMapper.h"

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
  GL::Texture2D _matcap{NoCreate};
  GL::Texture2D _clothTexture{NoCreate};
  bool _paused                   = false;
  int _min_loop_ms               = 16;
  bool _single_step              = false;
  bool _render_mesh              = false;
  bool _render_yarns             = true;
  bool _render_obstacles         = true;
  bool _rotate_scene             = false;
  std::string _matcap_file       = "matcaps/lighting1.jpg";
  std::string _matcapObs_file    = "matcaps/lighting1.jpg";
  std::string _clothtexture_file = "textures/colorgridy.jpg";
  float _render_radius_mult      = 1.0f;
  float _mesh_dz                 = 0.0f;
  float _clothTexture_scale      = 1.0f;
  Color4 _bgColor                = Color4(Color3(0.2f), 1.0f);

  // std::unique_ptr<ArcBallCamera> _arcballCamera;

  MeshShader _meshShader{NoCreate};
  ObsMeshShader _obsMeshShader{NoCreate};
  YarnShader _yarnGeometryShader{NoCreate};
  SsaoApplyShader _ssaoApplyShader{NoCreate};
  SsaoShader _ssaoShader{NoCreate};

  Containers::Optional<ArcBall> _arcball;
  Matrix4 _projection;
  Deg _proj_fov    = 45.0_degf;
  float _proj_near = 0.0001f;
  float _proj_far  = 100.0f;  // TODO reduce further?

  GL::Mesh _screenAlignedTriangle{NoCreate};

  GL::Framebuffer _framebuffer{NoCreate};

  GL::Texture2D _albedo{NoCreate};
  GL::Texture2D _positions{NoCreate};
  GL::Texture2D _normals{NoCreate};
  GL::Texture2D _occlusion{NoCreate};
  GL::Texture2D _noise{NoCreate};

  GL::Texture2D _depth{NoCreate};

  /* Profiling */
  DebugTools::GLFrameProfiler _profiler;

  // TODO hotkeys for various camera views and/or distances (with modifiers: no
  // mod both, ctrl view, alt dist)!
  // TODO make a rendersettings struct with default init
  // which can the also be used to reset the settings!
  Color4 _specularColor{0.3};        // TODO remove
  Float _ao_radius       = 0.004f;   // m
  Float _ao_bias         = 0.0003f;  // m
  int _ao_blur_radius    = 0;        // pixels
  float _ao_blur_feature = 25.0f;    // 1/m
  Float _ao_pow          = 2.0f;

  SsaoApplyShader::Flag _ssaoApplyFlag = {};

  ImGuiIntegration::Context _imgui{NoCreate};
};
}  // namespace Magnum

#endif  // __MAINAPPLICATION__H__

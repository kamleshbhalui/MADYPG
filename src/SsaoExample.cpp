#define DECLARE_UNUSED(x) ((void)x);

#define HELPSCALE 1.0f

#include <Magnum/Primitives/Axis.h>

#include "arcball/ArcBall.h"
#include "arcball/ArcBallCamera.h"

// #include "configure.h"
#include <Corrade/Containers/Optional.h>

#include "render/MeshDrawable.h"
#include "render/YarnDrawable.h"
#include "render/shaders/MeshShader.h"
#include "render/shaders/SsaoApplyShader.h"
#include "render/shaders/SsaoShader.h"
#include "render/shaders/YarnShader.h"
#include "yarns/YarnMapper.h"
// #include <Corrade/Containers/String.h>
// #include <Corrade/Containers/StringStl.h> /* until Corrade/Directory can
// handle Containers::String */
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
#include <Magnum/Shaders/Phong.h>
// #include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Trade/AbstractImporter.h>
#include <Magnum/Trade/ImageData.h>
#include <Magnum/Trade/MeshData.h>
#include <MagnumPlugins/AnyImageImporter/AnyImageImporter.h>
#include <MagnumPlugins/JpegImporter/JpegImporter.h>
#include <MagnumPlugins/PngImporter/PngImporter.h>

// #include <Magnum/SceneGraph/Camera.h>
// #include <Magnum/SceneGraph/Drawable.h>
// #include <Magnum/SceneGraph/MatrixTransformation3D.h>
// #include <Magnum/SceneGraph/Object.h>
// #include <Magnum/SceneGraph/Scene.h>

#include <imgui.h>

#include <Magnum/ImGuiIntegration/Context.hpp>
#include <memory>
#include <random>

namespace Magnum {
void setupTexture(GL::Texture2D &texture, Vector2i const &size,
                  GL::TextureFormat format);
void setupTexture(GL::Texture2D &texture, Vector2i const &size,
                  GL::TextureFormat format) {
  texture = GL::Texture2D{};
  texture.setMagnificationFilter(GL::SamplerFilter::Linear)
      .setMinificationFilter(GL::SamplerFilter::Linear)
      .setWrapping(GL::SamplerWrapping::ClampToEdge)
      .setStorage(1, format, size);
}

using namespace Math::Literals;

class SsaoExample : public Platform::Application {
 public:
  explicit SsaoExample(const Arguments &arguments);

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

  std::vector<YarnDrawable<YarnShader>> _yarnDrawable;
  std::unique_ptr<MeshDrawable<MeshShader>> _meshdrawable;
  std::unique_ptr<YarnMapper> _yarnMapper;
  GL::Texture2D _matcap{NoCreate};
  bool _paused              = false;
  bool _render_mesh         = true;
  bool _render_yarns        = true;
  float _render_radius_mult = 1.0f;
  Color4 _bgColor           = Color4(Color3(0.2f), 1.0f);

  // std::unique_ptr<ArcBallCamera> _arcballCamera;

  MeshShader _meshShader{NoCreate};
  YarnShader _yarnGeometryShader{NoCreate};
  SsaoApplyShader _ssaoApplyShader{NoCreate};
  SsaoShader _ssaoShader{NoCreate};

  Magnum::Shaders::Phong _phong{NoCreate};
  // Magnum::Shaders::VertexColor3D _vcShader{NoCreate};
  // GL::Mesh _axes = MeshTools::compile(Primitives::axis3D());

  Containers::Optional<ArcBall> _arcball;
  Matrix4 _projection;
  Deg _proj_fov = 45.0_degf;
  // float _proj_near = 0.01f;//0.001f *HELPSCALE;
  // float _proj_far = 10000.0f;//50.0f*HELPSCALE; // TODO reduce far/near to
  // proper scale
  float _proj_near = 0.0001f * HELPSCALE;  // 0.001f *HELPSCALE;
  float _proj_far =
      100.0f *
      HELPSCALE;  // 50.0f*HELPSCALE; // TODO reduce far/near to proper scale

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
  Color4 _specularColor{0.3};
  Float _ao_radius       = 0.004f * HELPSCALE;  // 1.5f;
  Float _ao_bias         = 0.0003f;             // 0.5f;
  int _ao_blur_radius    = 0;
  float _ao_blur_feature = 25.0f / HELPSCALE;
  Float _ao_pow          = 2.0f;

  SsaoApplyShader::Flag _ssaoApplyFlag = {};

  ImGuiIntegration::Context _imgui{NoCreate};
};

SsaoExample::SsaoExample(const Arguments &arguments)
    : Platform::Application{arguments, NoCreate} {
  std::string meshPath;
  {
    // if(Utility::Directory::exists(SSAO_EXAMPLE_DIR))
    //     meshPath = SSAO_EXAMPLE_DIR;
    // else if(Utility::Directory::exists(SSAO_EXAMPLE_INSTALL_DIR))
    //     meshPath = SSAO_EXAMPLE_INSTALL_DIR;
    // else
    meshPath =
        Utility::Directory::path(Utility::Directory::executableLocation());
    meshPath = Utility::Directory::join(meshPath, "../../src/Armadillo.ply");
    // meshPath = Utility::Directory::join(meshPath, "Armadillo.ply");

    /* Finally, provide a way for the user to override the model directory */
    Utility::Arguments args;
    args.addFinalOptionalArgument("mesh", meshPath)
        .setHelp("mesh", "Path to the mesh you want to import")
        .addSkippedPrefix("magnum", "engine-specific options")
        .setGlobalHelp(
            "Press P to toggle between phong shading and phong shading + ssao\n"
            "Press H in debug mode to hot reload the shaders\n"
            "Press C in debug mode to use a compute shader in the ssao pass\n"
            "Press R to reset the camera\n"
            "Press L to toggle lagging for the camera controls\n")
        .parse(arguments.argc, arguments.argv);
    /* relative paths are brittle, so prepend CWD to them if needed */
    meshPath = Utility::Directory::join(Utility::Directory::current(),
                                        args.value("mesh"));
  }

  /* Setup window */
  {
    const Vector2 dpiScaling = this->dpiScaling({});
    Configuration conf;
    conf.setTitle("Magnum SSAO Example")
        .setSize(conf.size(), dpiScaling)
        .setWindowFlags(Configuration::WindowFlag::Resizable);
    GLConfiguration glConf;
    glConf.setSampleCount(0);
    if (!tryCreate(conf, glConf)) {
      create(conf, glConf.setSampleCount(0));
    }
  }

  GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
  GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

  /* setup imgui context and proper blending */
  {
    _imgui = ImGuiIntegration::Context(Vector2{windowSize()} / dpiScaling(),
                                       windowSize(), framebufferSize());
    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
                                   GL::Renderer::BlendEquation::Add);
    GL::Renderer::setBlendFunction(
        GL::Renderer::BlendFunction::SourceAlpha,
        GL::Renderer::BlendFunction::OneMinusSourceAlpha);
  }

  /* Load Mesh, setup textures and framebuffer*/
  {
    _screenAlignedTriangle = GL::Mesh{};
    _screenAlignedTriangle.setCount(3);

    Containers::Array<Vector4> noise(Containers::NoInit, 16);
    std::random_device device;
    std::default_random_engine engine(device());
    // uniform sampling: square
    // std::uniform_real_distribution<float> distr(-1, 1);
    // for (Vector4 &n : noise)
    //   n = Vector4{distr(engine), distr(engine), 0, 0};
    // uniform sampling: circle
    std::normal_distribution<float> ndistr(0.0f, 0.4f);
    for (Vector4 &n : noise) {
      n = Vector4{ndistr(engine), ndistr(engine), 0, 0}.normalized();
    }

    _noise = GL::Texture2D{};
    ImageView2D view{PixelFormat::RGBA32F, {4, 4}, noise};
    _noise.setMagnificationFilter(GL::SamplerFilter::Linear)
        .setMinificationFilter(GL::SamplerFilter::Linear)
        .setWrapping(GL::SamplerWrapping::Repeat)
        .setStorage(1, GL::TextureFormat::RGBA32F, {4, 4})
        .setSubImage(0, {}, view);

    const Range2Di viewport = GL::defaultFramebuffer.viewport();
    const Vector2i vpSize   = viewport.size();

    _framebuffer = GL::Framebuffer{viewport};
    setupFramebuffer(vpSize);
    // GL::Renderer::setClearColor({});
    // GL::Renderer::setClearColor(_bgColor);
    // GL::defaultFramebuffer.clearColor(_bgColor);

    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _ssaoShader         = SsaoShader{};
    _ssaoApplyShader    = SsaoApplyShader{};
    _phong              = Magnum::Shaders::Phong{};

    {
      _yarnMapper = std::make_unique<YarnMapper>("models/model_rib/");
      _yarnDrawable.emplace_back(_yarnGeometryShader);
      _yarnDrawable.back().setIndices(_yarnMapper->getIndices());
      _yarnDrawable.back().setVertices(_yarnMapper->getVertexData());
      _yarnDrawable.back().m_radius =
          _yarnMapper->getRadius() * _render_radius_mult;

      _meshdrawable = std::make_unique<MeshDrawable<MeshShader>>(_meshShader);
      const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
      _meshdrawable->setIndices(mesh.F);
      _meshdrawable->setVertices(mesh.X);
    }

    PluginManager::Manager<Trade::AbstractImporter> manager;
    Trade::AnyImageImporter importer = Trade::AnyImageImporter(manager);
    if (!importer.openFile("matcaps/mcb.jpg"))
      // if (!importer.openFile("matcaps/hughsk/00036.png"))
      // if (!importer.openFile("matcaps/mua/test_gold.jpg"))
      std::exit(2);
    Containers::Optional<Trade::ImageData2D> image = importer.image2D(0);
    CORRADE_INTERNAL_ASSERT(image);
    _matcap = GL::Texture2D();
    _matcap.setWrapping(GL::SamplerWrapping::ClampToEdge)
        .setMagnificationFilter(GL::SamplerFilter::Linear)
        .setMinificationFilter(GL::SamplerFilter::Linear,
                               GL::SamplerMipmap::Linear)
        .setStorage(Math::log2(image->size().min()) + 1,
                    GL::textureFormat(image->format()), image->size())
        .setSubImage(0, {}, *image)
        .generateMipmap();
  }

  /* Set up the arcball and projection */
  {
    const Vector3 eye = Vector3::zAxis(1.0f) * HELPSCALE;
    const Vector3 center{};
    const Vector3 up = Vector3::yAxis();
    _arcball.emplace(eye, center, up, 45.0_degf, windowSize());
    _arcball->setLagging(0.85f);

    _projection = Matrix4::perspectiveProjection(
        _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near,
        _proj_far);
  }

  _profiler = DebugTools::GLFrameProfiler{
      DebugTools::GLFrameProfiler::Value::FrameTime |
          DebugTools::GLFrameProfiler::Value::GpuDuration |
          DebugTools::GLFrameProfiler::Value::CpuDuration,
      60};

  /* Loop at 60 Hz max */
  setSwapInterval(1);
  setMinimalLoopPeriod(16);
}

void SsaoExample::drawEvent() {
  _profiler.beginFrame();

  if (!_paused) {  // SIM
    _yarnMapper->step();
    // assume no update to yarn indices
    _yarnDrawable.back().setVertices(_yarnMapper->getVertexData());

    const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
    if (_yarnMapper->getMeshSimulation()->meshIndicesDirty())
      _meshdrawable->setIndices(mesh.F);
    // always assume changes to vertices
    _meshdrawable->setVertices(mesh.X);
  }

  const bool camChanged = _arcball->updateTransformation();
  DECLARE_UNUSED(camChanged);

  GL::defaultFramebuffer.clear(GL::FramebufferClear::Color |
                               GL::FramebufferClear::Depth);
  bool require_redraw = !_paused || camChanged;
  if (require_redraw) {
    const Matrix4 tf = _arcball->viewMatrix();

    /* render the scene into g-buffer */
    _framebuffer
        .mapForDraw(
            {{YarnShader::AlbedoOutput, GL::Framebuffer::ColorAttachment{0}},
             {YarnShader::PositionsOutput, GL::Framebuffer::ColorAttachment{1}},
             {YarnShader::NormalsOutput, GL::Framebuffer::ColorAttachment{2}}})
        // .clear(GL::FramebufferClear::Depth | GL::FramebufferClear::Color)
        .clearColor(0, _bgColor)
        .clearColor(1, Color4(0.0, 0.0, 0.0, 1.0))
        .clearColor(2, Color4(0.0, 0.0, 1.0, 1.0))
        .clearDepth(1.0)
        .bind();

    if (_render_yarns) {
      _yarnGeometryShader.bindTexture(_matcap);
      _yarnGeometryShader.setProjection(_projection);
      _yarnDrawable.back().m_radius =
          _yarnMapper->getRadius() * _render_radius_mult;
      for (auto &line : _yarnDrawable) line.draw(tf);
    }

    if (_render_mesh) {
      _meshShader.setProjection(_projection);
      _meshdrawable->draw(tf);
    }

    _ssaoShader.bindNormalTexture(_normals)
        .bindNoiseTexture(_noise)
        .bindPositionTexture(_positions)
        .setProjectionMatrix(_projection)
        .setSampleRadius(_ao_radius)
        .setBias(_ao_bias);

    _framebuffer
        .mapForDraw({{SsaoShader::AmbientOcclusionOutput,
                      GL::Framebuffer::ColorAttachment{3}}})
        .clear(GL::FramebufferClear::Color);
    _ssaoShader.draw(_screenAlignedTriangle);
  }
  GL::defaultFramebuffer.bind();
  _ssaoApplyShader.bindAlbedoTexture(_albedo)
      .bindOcclusionTexture(_occlusion)
      .bindNormalTexture(_normals)
      .bindPositionTexture(_positions)
      .setLightPosition({5.0f, 5.0f, 7.0f})
      .setLightColor(Color3{1.f})
      .setShininess(80)
      .setSpecularColor(_specularColor.rgb())
      .setAOBlurRadius(_ao_blur_radius)
      .setAOBlurFeature(_ao_blur_feature)
      .setAOPow(_ao_pow)
      .draw(_screenAlignedTriangle);

  _profiler.endFrame();

  _imgui.newFrame();
  if (ImGui::GetIO().WantTextInput && !isTextInputActive())
    startTextInput();
  else if (!ImGui::GetIO().WantTextInput && isTextInputActive())
    stopTextInput();

  drawSettings();
  _imgui.updateApplicationCursor(*this);

  GL::Renderer::enable(GL::Renderer::Feature::Blending);
  GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
  GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
  GL::Renderer::disable(GL::Renderer::Feature::DepthTest);

  _imgui.drawFrame();

  GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
  GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
  GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
  GL::Renderer::disable(GL::Renderer::Feature::Blending);

  swapBuffers();

  redraw();
}

void SsaoExample::viewportEvent(ViewportEvent &event) {
  const Vector2i fbSize = event.framebufferSize();
  const Vector2i wSize  = event.windowSize();

  GL::defaultFramebuffer.setViewport({{}, fbSize});
  _framebuffer.setViewport({{}, fbSize});

  _arcball->reshape(wSize);
  _projection = Matrix4::perspectiveProjection(
      _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near,
      _proj_far);

  setupFramebuffer(fbSize);

  _imgui.relayout(Vector2{wSize} / event.dpiScaling(), event.windowSize(),
                  fbSize);
}

void SsaoExample::setupFramebuffer(const Vector2i &size) {
  setupTexture(_albedo, size, GL::TextureFormat::RGBA32F);
  setupTexture(_positions, size, GL::TextureFormat::RGBA32F);
  setupTexture(_normals, size, GL::TextureFormat::RGBA32F);
  setupTexture(_occlusion, size, GL::TextureFormat::R32F);
  setupTexture(_depth, size, GL::TextureFormat::DepthComponent32F);

  _framebuffer
      .attachTexture(GL::Framebuffer::BufferAttachment::Depth, _depth, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _albedo, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{3}, _occlusion, 0);

  CORRADE_INTERNAL_ASSERT(
      _framebuffer.checkStatus(GL::FramebufferTarget::Draw) ==
      GL::Framebuffer::Status::Complete);
}

void SsaoExample::drawSettings() {
  ImGui::Begin("Render Options");
  const float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));

  if (ImGui::CollapsingHeader("SSAO", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();
    ImGui::PushItemWidth(100.0f);
    {
      float val = _ao_radius * 100.0f;
      if (ImGui::SliderFloat("radius (cm)", &val, 0.001f, 5.0f))
        _ao_radius = val * 0.01f;
    }
    {
      float val = _ao_bias * 100.0f;
      if (ImGui::SliderFloat("bias (cm)", &val, 0.001f, 0.5f))
        _ao_bias = val * 0.01f;
    }
    ImGui::SliderInt("blur radius", &_ao_blur_radius, 0, 5);
    //  (ImGui::DragFloat("blur feature (cm)", &_ao_blur_feature, 0.01f,
    //  0.0f,100.0f,"%.2e"));
    {
      float val = _ao_blur_feature * 0.01f;
      if (ImGui::SliderFloat("blur feature (1/cm)", &val, 0.0f, 50.0f))
        _ao_blur_feature = val * 100.0f;
    }
    ImGui::SliderFloat("strength", &_ao_pow, 0.0f, 10.0f);
    ImGui::PopItemWidth();
    ImGui::Unindent();
  }

  ImGui::Checkbox("Yarns", &_render_yarns);
  ImGui::SameLine();
  ImGui::Checkbox("Mesh", &_render_mesh);

  {
    ImGui::PushItemWidth(100.0f);
    ImGui::DragFloat("yarn radius mult", &_render_radius_mult, 0.01f, 0.0f,
                     2.0f);
    ImGui::PopItemWidth();
  }

  ImGui::Checkbox("Pause", &_paused);

  ImGui::SameLine();

  static bool drawOcclusion = false;
  if (ImGui::Checkbox("Show Occlusion Factor", &drawOcclusion)) {
    if (drawOcclusion)
      _ssaoApplyFlag = SsaoApplyShader::Flag::DrawAmbientOcclusion;
    else
      _ssaoApplyFlag = {};
    _ssaoApplyShader = SsaoApplyShader{_ssaoApplyFlag};
  }

  if (ImGui::Button("Hot Reload Shader")) {
    Utility::Resource::overrideGroup("ssao-data",
                                     "src/render/shaders/resources.conf");
    _ssaoShader      = SsaoShader{};
    _ssaoApplyShader = SsaoApplyShader{_ssaoApplyFlag};
  }

  if (ImGui::Button("Reset Camera"))
    _arcball->reset();

  ImGui::SameLine();
  static bool lagging = true;
  if (ImGui::Checkbox("Camera Lagging", &lagging)) {
    _arcball->setLagging(float(lagging) * 0.75f);
  }

  {
    static float col[3] = {_bgColor.r(), _bgColor.g(), _bgColor.b()};
    if (ImGui::ColorEdit3("BG Color", col)) {
      _bgColor = Color4(col[0], col[1], col[2], 1.0f);
      // GL::Renderer::setClearColor(_bgColor);
    }
  }

  std::string stats = _profiler.statistics();
  // ImGui::TextUnformatted(stats.begin(), stats.end());
  ImGui::TextUnformatted(stats.c_str());

  ImGui::PopStyleVar();
  ImGui::End();
}

void SsaoExample::keyPressEvent(KeyEvent &event) {
  if (_imgui.handleKeyPressEvent(event))
    return;

  switch (event.key()) {
    case KeyEvent::Key::Esc:
      this->exit();
      break;
    case KeyEvent::Key::Space:
      _paused = !_paused;
      break;
    default:
      break;
  }
}

void SsaoExample::keyReleaseEvent(KeyEvent &event) {
  if (_imgui.handleKeyReleaseEvent(event))
    return;
}

void SsaoExample::textInputEvent(TextInputEvent &event) {
  if (_imgui.handleTextInputEvent(event))
    return;
}

void SsaoExample::mousePressEvent(MouseEvent &event) {
  if (_imgui.handleMousePressEvent(event))
    return;
  /* Enable mouse capture so the mouse can drag outside of the window */
  /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
  SDL_CaptureMouse(SDL_TRUE);
  _arcball->initTransformation(event.position());
  event.setAccepted();
  redraw();
}

void SsaoExample::mouseReleaseEvent(MouseEvent &event) {
  if (_imgui.handleMouseReleaseEvent(event))
    return;
  /* Disable mouse capture again */
  /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
  SDL_CaptureMouse(SDL_FALSE);
}

void SsaoExample::mouseMoveEvent(MouseMoveEvent &event) {
  if (_imgui.handleMouseMoveEvent(event))
    return;
  if (!event.buttons())
    return;

  if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
    _arcball->translate(event.position());
  else
    _arcball->rotate(event.position());

  event.setAccepted();
  redraw();
}

void SsaoExample::mouseScrollEvent(MouseScrollEvent &event) {
  if (_imgui.handleMouseScrollEvent(event)) {
    /* Prevent scrolling the page */
    event.setAccepted();
    return;
  }

  const Float delta = event.offset().y();
  if (Math::abs(delta) < 1.0e-2f)
    return;

  _arcball->zoom(delta * 0.1f * HELPSCALE);

  event.setAccepted();
  redraw();
}
}  // namespace Magnum

MAGNUM_APPLICATION_MAIN(Magnum::SsaoExample)

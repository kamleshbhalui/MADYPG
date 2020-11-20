#include "MainApplication.h"

#include <imgui_stdlib.h>

// #ifdef MSAA
// #define SSAO_SAMPLES 32
// #else
#define SSAO_SAMPLES 64
// #endif

using namespace Magnum;
using namespace Math::Literals;

void setupTexture(GL::Texture2D &texture, Magnum::Vector2i const &size,
                  GL::TextureFormat format);
void setupTexture(GL::Texture2D &texture, Magnum::Vector2i const &size,
                  GL::TextureFormat format) {
  texture = GL::Texture2D{};
  texture.setMagnificationFilter(GL::SamplerFilter::Linear)
      .setMinificationFilter(GL::SamplerFilter::Linear)
      .setWrapping(GL::SamplerWrapping::ClampToEdge)
      .setStorage(1, format, size);
}
#ifdef MSAA
void setupTexture(GL::MultisampleTexture2D &texture,
                  Magnum::Vector2i const &size, GL::TextureFormat format);
void setupTexture(GL::MultisampleTexture2D &texture,
                  Magnum::Vector2i const &size, GL::TextureFormat format) {
  texture = GL::MultisampleTexture2D{};
  texture.setStorage(MSAA, format, size);
}
#endif

void makeDefaultTexture(GL::Texture2D &tex2D);
void makeDefaultTexture(GL::Texture2D &tex2D) {
  Containers::Array<Vector4> data(Containers::NoInit, 1);
  data[0] = Vector4(1.0f);
  ImageView2D view{PixelFormat::RGBA32F, {1, 1}, data};

  tex2D = GL::Texture2D{};
  tex2D.setMagnificationFilter(GL::SamplerFilter::Linear)
      .setMinificationFilter(GL::SamplerFilter::Linear)
      .setWrapping(GL::SamplerWrapping::Repeat)
      .setStorage(1, GL::TextureFormat::RGBA32F, view.size())
      .setSubImage(0, {}, view);
}

#include <functional>
namespace ImGui {
bool TextBrowser(
    ImGui::FileBrowser &browser, std::string &txt,
    std::function<void()> onChange = []() {});
bool TextBrowser(ImGui::FileBrowser &browser, std::string &txt,
                 std::function<void()> onChange) {
  ImGui::PushID(&txt);
  ImGui::PushItemWidth(150);
  if (ImGui::InputText("##txt", &txt, ImGuiInputTextFlags_EnterReturnsTrue)) {
    onChange();
    ImGui::PopItemWidth();
    ImGui::SameLine();
    ImGui::Button("Browse");
    ImGui::PopID();
    return true;
  }
  ImGui::PopItemWidth();
  ImGui::SameLine();
  if (ImGui::Button("Browse")) {
    browser.SetPwd(std::filesystem::path(txt).parent_path());
    browser.Open();
  }
  browser.Display();
  if (browser.HasSelected()) {
    txt = std::filesystem::relative(browser.GetSelected().string());
    browser.ClearSelected();
    onChange();
    ImGui::PopID();
    return true;
  }
  ImGui::PopID();
  return false;
}
}  // namespace ImGui

bool loadTexture(const std::string &file, GL::Texture2D &tex2D,
                 GL::SamplerWrapping wrapping = GL::SamplerWrapping::Repeat);
bool loadTexture1D(const std::string &file, GL::Texture1D &tex1D,
                 GL::SamplerWrapping wrapping = GL::SamplerWrapping::Repeat);
bool loadTexture(const std::string &file, GL::Texture2D &tex2D,
                 GL::SamplerWrapping wrapping) {
  PluginManager::Manager<Trade::AbstractImporter> manager;
  Trade::AnyImageImporter importer = Trade::AnyImageImporter(manager);
  if (!importer.openFile(file)) {
    ::Debug::error("Couldn't load texture:", file);
    makeDefaultTexture(tex2D);
    return false;
  }
  Containers::Optional<Trade::ImageData2D> image = importer.image2D(0);
  CORRADE_INTERNAL_ASSERT(image);
  if (!image) {
    ::Debug::error("Couldn't load texture image:", file);
    makeDefaultTexture(tex2D);
    return false;
  }
  tex2D = GL::Texture2D();
  tex2D.setWrapping(wrapping)
      .setMagnificationFilter(GL::SamplerFilter::Linear)
      .setMinificationFilter(GL::SamplerFilter::Linear,
                             GL::SamplerMipmap::Linear)
      .setStorage(Math::log2(image->size().min()) + 1,
                  GL::textureFormat(image->format()), image->size())
      .setSubImage(0, {}, *image)
      .generateMipmap();
  return true;
}
bool loadTexture1D(const std::string &file, GL::Texture1D &tex1D,
                 GL::SamplerWrapping wrapping) {
  PluginManager::Manager<Trade::AbstractImporter> manager;
  Trade::AnyImageImporter importer = Trade::AnyImageImporter(manager);
  if (!importer.openFile(file)) {
    ::Debug::error("Couldn't load texture:", file);
    return false;
  }
  Containers::Optional<Trade::ImageData2D> image = importer.image2D(0);
  CORRADE_INTERNAL_ASSERT(image);
  if (!image) {
    ::Debug::error("Couldn't load texture image:", file);
    return false;
  }
  // 1D view on the first row
  ImageView1D image1D{image->format(), image->size().x(), image->data()};
  tex1D = GL::Texture1D();
  tex1D.setWrapping(wrapping)
      .setMagnificationFilter(GL::SamplerFilter::Linear)
      .setMinificationFilter(GL::SamplerFilter::Linear,
                             GL::SamplerMipmap::Linear)
      .setStorage(Math::log2(image1D.size()[0]) + 1,
                  GL::textureFormat(image1D.format()), image1D.size())
      .setSubImage(0, {}, image1D)
      .generateMipmap();
  return true;
}



void MainApplication::reset_simulation() {
  if (!_yarnMapper)
    _yarnMapper = std::make_unique<YarnMapper>();
  else
    _yarnMapper.reset(new YarnMapper());

  _yarnMapper->m_settings = _yarnMapperSettings;

  _yarnMapper->initialize();
  _yarnDrawable.clear();
  _yarnDrawable.emplace_back(_yarnGeometryShader,
                             _yarnMapper->getVertexBuffer(),
                             _yarnMapper->getIndexBuffer());

  // _yarnDrawable.back().setIndices(_yarnMapper->getIndices());
  // _yarnDrawable.back().setVertices(_yarnMapper->getVertexData());
  _yarnDrawable.back().m_radius =
      _yarnMapper->getRadius() * _render_radius_mult;
  _yarnDrawable.back().m_nmtwist = _render_nmtwist;
  _yarnDrawable.back().m_nmnum = _render_nmnum;
  _yarnDrawable.back().m_nmheight = _render_nmheight;

  auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
  _meshdrawable.release();
  _meshdrawable =
      std::make_unique<MeshDrawable<MeshShader>>(_meshShader, mesh.F, mesh.X);
  // _meshdrawable->setIndices(mesh.F);
  // _meshdrawable->setVertices(mesh.X);

  _obsmeshdrawables.clear();
  for (const auto &obs : _yarnMapper->getMeshSimulation()->getObstacles()) {
    _obsmeshdrawables.emplace_back(_obsMeshShader);
    _obsmeshdrawables.back().setIndices(obs.mesh.F);
    _obsmeshdrawables.back().setVertices(obs.mesh.X);
  }
}

MainApplication::MainApplication(const Arguments &arguments)
    : Platform::Application{arguments, NoCreate} {
#ifdef MSAA
  ::Debug::logf("Using %dx MSAA\n", MSAA);
  GL::Renderer::enable(GL::Renderer::Feature::Multisampling);
#endif

  /* Setup window */
  {
    const Vector2 dpiScaling = this->dpiScaling({});
    Configuration conf;
    conf.setSize({1200, 800}, dpiScaling);
    conf.setTitle("Mesh 2 Yarns")
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

    _fbo_gbuffer = GL::Framebuffer{viewport};
#ifdef MSAA
    _fbo_ssao = GL::Framebuffer{viewport};
#endif
    setupFramebuffer(vpSize);
    // GL::Renderer::setClearColor({});
    // GL::Renderer::setClearColor(_bgColor);
    // GL::defaultFramebuffer.clearColor(_bgColor);

    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _obsMeshShader      = ObsMeshShader{};
    _ssaoShader         = SsaoShader{SSAO_SAMPLES};
    _ssaoApplyShader    = SsaoApplyShader{};

    {  // default settings
      _yarnMapperSettings.modelfolder = "models/model_stock";
      // _yarnMapperSettings.modelfolder = "models/V0/model_rib";
      _yarnMapperSettings.provider_type =
          YarnMapper::Settings::Provider::BinSeq;
      // _yarnMapperSettings.binseq_settings.filepath = "binseqs/sock.bin";
      // #ifdef NDEBUG
      _yarnMapperSettings.binseq_settings.filepath =
          "binseqs/sxsy_const.bin";
          // "binseqs/HYLC_30x30/basket_drapeX.bin";
      // #else
      // _yarnMapperSettings.binseq_settings.filepath =
      // "binseqs/sxsy_const.bin"; #endif
      _yarnMapperSettings.objseq_settings.folder = "objseqs/sxsy_const";
      _yarnMapperSettings.objseq_settings.constant_material_space = true;
    }
    _folderDialog = std::make_unique<ImGui::FileBrowser>(
        ImGuiFileBrowserFlags_SelectDirectory |
        ImGuiFileBrowserFlags_CloseOnEsc);
    _folderDialog->SetTitle("Select folder");
    _fileDialog =
        std::make_unique<ImGui::FileBrowser>(ImGuiFileBrowserFlags_CloseOnEsc);
    _fileDialog->SetTitle("Select file");

    reset_simulation();

    // PluginManager::Manager<Trade::AbstractImporter> manager;
    // Trade::AnyImageImporter importer = Trade::AnyImageImporter(manager);
    loadTexture(_matcap_file, _matcap, GL::SamplerWrapping::ClampToEdge);
    loadTexture(_matcapObs_file, _matcapObs, GL::SamplerWrapping::ClampToEdge);
    loadTexture(_clothtexture_file, _clothTexture, GL::SamplerWrapping::Repeat);
    loadTexture1D(_normalMap_file, _normalMap, GL::SamplerWrapping::Repeat);
  }

  /* Set up the arcball and projection */
  {
    const Vector3 eye =
        2.0f * Vector3(+0.4f, 0.3f, 0.5f);  // Vector3::yAxis(1.0f);
    const Vector3 center{};
    const Vector3 up = Vector3::yAxis();
    ;
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
  setSwapInterval(0);  // 0 no vsync, 1 vsync
  setMinimalLoopPeriod(_min_loop_ms);
}

void MainApplication::drawEvent() {
  _profiler.beginFrame();

  bool simChanged = false;
  if (!_paused || _single_step) {  // SIM
    if (_yarnMapper->isInitialized()) {
      // force update some settings
      _yarnMapper->m_settings = _yarnMapperSettings;
      _yarnMapper->step();

      const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
      // if (_yarnMapper->getMeshSimulation()->meshIndicesDirty())
      _meshdrawable->updateIndexCount(mesh.F.getGPUSize() * 3);
    }
    _single_step = false;
    simChanged   = true;
  }

  const bool camChanged = _arcball->updateTransformation();
  DECLARE_UNUSED(camChanged);

  GL::defaultFramebuffer.clear(GL::FramebufferClear::Color |
                               GL::FramebufferClear::Depth);
  bool require_redraw = simChanged || camChanged;
  if (require_redraw) {
    Matrix4 tf = _arcball->viewMatrix();
    if (_rotate_scene) {
      tf = tf * Matrix4::rotation(Rad(-1.57079632679f), Vector3::xAxis());
    }

    /* render the scene into g-buffer */
    _fbo_gbuffer
        .mapForDraw(
            {{YarnShader::AlbedoOutput, GL::Framebuffer::ColorAttachment{0}},
             {YarnShader::PositionsOutput, GL::Framebuffer::ColorAttachment{1}},
             {YarnShader::NormalsOutput, GL::Framebuffer::ColorAttachment{2}}})
        // .clear(GL::FramebufferClear::Depth | GL::FramebufferClear::Color)
        .clearColor(0, _bgColor)
        .clearColor(1, Color4(0.0, 0.0, -1000.0, 1.0)) // NOTE: using arbitrary -1000 as background depth
        .clearColor(2, Color4(0.0, 0.0, 1.0, 1.0))
        .clearDepth(1.0)
        .bind();

    if (_render_yarns && _yarnMapper->isInitialized()) {
      _yarnGeometryShader.bindMatCap(_matcap);
      _yarnGeometryShader.bindClothTexture(_clothTexture);
      _yarnGeometryShader.bindNormalMap(_normalMap);
      _yarnGeometryShader.setProjection(_projection);
      _yarnGeometryShader.setTextureScale(_clothTexture_scale);
      _yarnDrawable.back().m_radius =
          _yarnMapper->getRadius() * _render_radius_mult;
      _yarnDrawable.back().m_nmtwist = _render_nmtwist;
      _yarnDrawable.back().m_nmnum = _render_nmnum;
      _yarnDrawable.back().m_nmheight = _render_nmheight;
      for (auto &line : _yarnDrawable) line.draw(tf);
    }

    if (_render_mesh) {
      _meshShader.setProjection(_projection);
      _meshShader.setDZ(_mesh_dz);
      _meshdrawable->draw(tf);
    }

    if (_render_obstacles) {
      _obsMeshShader.bindMatCap(_matcapObs);
      _obsMeshShader.setProjection(_projection);

      const auto &obs = _yarnMapper->getMeshSimulation()->getObstacles();
      ::Debug::msgassert(
          "Different number of obstacles and obstacle drawables!",
          obs.size() == _obsmeshdrawables.size());
      for (size_t i = 0; i < obs.size(); i++) {
        _obsmeshdrawables[i].draw(tf, obs[i].transformation);
      }
    }


    _ssaoShader.bindNormalTexture(_normals)
        .bindNoiseTexture(_noise)
        .bindPositionTexture(_positions)
        .setProjectionMatrix(_projection)
        .setSampleRadius(_ao_radius)
        .setBias(_ao_bias);

    #ifdef MSAA
    _fbo_ssao
        .mapForDraw({{SsaoShader::AmbientOcclusionOutput,
                      GL::Framebuffer::ColorAttachment{0}}})
        .clear(GL::FramebufferClear::Color)
        // .clearColor(3, Color4(0.0, 0.0, 1.0, 1.0))
        // .clear(GL::FramebufferClear::Color) // TODO? WHAT IS COLORATTACHMENT, replace with clearColor(3, ...)
        .bind();
        
  // ::Debug::error("TEST",_fbo_ssao.checkStatus(GL::FramebufferTarget::Draw) ==
  //                     GL::Framebuffer::Status::Complete);
    _ssaoShader.draw(_screenAlignedTriangle);
    #else
    _fbo_gbuffer
        .mapForDraw({{SsaoShader::AmbientOcclusionOutput,
                      GL::Framebuffer::ColorAttachment{3}}})
        .clear(GL::FramebufferClear::Color);
    _ssaoShader.draw(_screenAlignedTriangle);
    #endif
  }

  GL::defaultFramebuffer.bind();
  _ssaoApplyShader.bindAlbedoTexture(_albedo)
      .bindOcclusionTexture(_occlusion)
      // .bindNormalTexture(_normals)
      .bindPositionTexture(_positions)
      .setLightPosition({5.0f, 5.0f, 7.0f})
      .setLightColor(Color3{1.f})
      .setShininess(80)
      .setSpecularColor(_specularColor.rgb())
      .setAOBlurRadius(_ao_blur_radius)
      .setAOBlurFeature(_ao_blur_feature)
      .setAOPow(_ao_pow)
      .draw(_screenAlignedTriangle);

    
    // GL::AbstractFramebuffer::blit(_fbo_gbuffer, GL::defaultFramebuffer, {Vector2i(0,0),_positions.imageSize()}, {Vector2i(0,0),_positions.imageSize()},GL::FramebufferBlit::Color,GL::FramebufferBlitFilter::Nearest); // TODO CONTINUE HERE. THIS SEEMS TO WORK, SO GBUFFER FBO AT LEAST HAS COLOR CORRECT (and probably the rest.) so at what stage is it failing then. is the ssao shader not binding gbuffer correct or reading its values? try to blit the output of ssao   
    // it also doesnt seem to work to blit gbuffer into ssao and then into default. so maybe fbossao is not even set up correctly (like to have a color thing? or bc blit color assumes that it wants to use colorattachment 0 which doesnt exist for ssaofbo?)
    // i can map any of gbuffer stuff for read into the default framebuffer with working MSAA
    // however going over ssaofbo does not work, maybe again i dont know to which attachment to draw 

    // _fbo_gbuffer.mapForRead(GL::Framebuffer::ColorAttachment{2});
    // GL::AbstractFramebuffer::blit(_fbo_gbuffer, GL::defaultFramebuffer, {Vector2i(0,0),_positions.imageSize()}, GL::FramebufferBlit::Color); 
    
    // maybe cant blit bc different format texture RGB->R?
    // _fbo_gbuffer.mapForRead(GL::Framebuffer::ColorAttachment{2});
    // _fbo_ssao.mapForDraw({{SsaoShader::AmbientOcclusionOutput,
    //                   GL::Framebuffer::ColorAttachment{3}}});
    // GL::AbstractFramebuffer::blit(_fbo_gbuffer, _fbo_ssao, {Vector2i(0,0),_positions.imageSize()}, GL::FramebufferBlit::Color); 


    // _fbo_ssao.mapForRead(GL::Framebuffer::ColorAttachment{3});    GL::defaultFramebuffer.mapForDraw(GL::DefaultFramebuffer::DrawAttachment::Back);
    // GL::AbstractFramebuffer::blit(_fbo_ssao, GL::defaultFramebuffer, {Vector2i(0,0),_positions.imageSize()},GL::FramebufferBlit::Color);
/*
0th attempt: all below but use texelfetch of fixed sample 0


1st attempt

--enable opengl multisample
--make albedo normal and position multisample
--keep aotexture single sample

--make second framebuffer for ssaoshading

ssao shader: use texelfetch to get avg position, avg normal, for neighborlookup use some random single sample position (eg based on the random normal vector using tangent t or so... biased?)
  ~~~~ position/normal

ssaoapply shader: use texelfetch for avg color, (if blur use texelfetch for random or average position)
  try just getting and showing occlusion first using occludraw
  ~~~~ albedo/position

use second fbo

define MSAA somewhere else


multisampled depth??
*/


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

  _profiler.endFrame();

  redraw();
}

void MainApplication::viewportEvent(ViewportEvent &event) {
  const Vector2i fbSize = event.framebufferSize();
  const Vector2i wSize  = event.windowSize();

  GL::defaultFramebuffer.setViewport({{}, fbSize});
  _fbo_gbuffer.setViewport({{}, fbSize});
  #ifdef MSAA
  _fbo_ssao.setViewport({{},fbSize});
  #endif

  _arcball->reshape(wSize);
  _projection = Matrix4::perspectiveProjection(
      _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near,
      _proj_far);

  setupFramebuffer(fbSize);

  _imgui.relayout(Vector2{wSize} / event.dpiScaling(), event.windowSize(),
                  fbSize);
}

void MainApplication::setupFramebuffer(const Vector2i &size) {
  setupTexture(_albedo, size, GL::TextureFormat::RGBA32F);
  setupTexture(_positions, size, GL::TextureFormat::RGBA32F);
  setupTexture(_normals, size, GL::TextureFormat::RGBA32F);
  setupTexture(_occlusion, size, GL::TextureFormat::R32F);
  setupTexture(_depth, size, GL::TextureFormat::DepthComponent32F);

// _depth.setStorage(128,GL::TextureFormat::DepthComponent32F,size); TODO WHY DOESNT THIS THROW AN ERROR, MAKING VARZING
#ifdef MSAA
  _fbo_gbuffer
      .attachTexture(GL::Framebuffer::BufferAttachment::Depth, _depth)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _albedo)
      .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions)
      .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals);
#else
  _fbo_gbuffer
      .attachTexture(GL::Framebuffer::BufferAttachment::Depth, _depth, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _albedo, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{3}, _occlusion,
                             0);
#endif
  CORRADE_INTERNAL_ASSERT(
      _fbo_gbuffer.checkStatus(GL::FramebufferTarget::Draw) ==
      GL::Framebuffer::Status::Complete);

#ifdef MSAA
// TODO maybe only attach the textures that need to be drawn?
  _fbo_ssao
      // .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions)
      // .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals)
      .attachTexture(GL::Framebuffer::ColorAttachment{0},  _occlusion, 0); // TEST ATTACHMENT {0}
  CORRADE_INTERNAL_ASSERT(_fbo_gbuffer.checkStatus(GL::FramebufferTarget::Read) ==
                          GL::Framebuffer::Status::Complete);
  CORRADE_INTERNAL_ASSERT(_fbo_ssao.checkStatus(GL::FramebufferTarget::Draw) ==
                          GL::Framebuffer::Status::Complete);
  CORRADE_INTERNAL_ASSERT(_fbo_ssao.checkStatus(GL::FramebufferTarget::Read) ==
                          GL::Framebuffer::Status::Complete);
  // ::Debug::error("TEST",_fbo_ssao.checkStatus(GL::FramebufferTarget::Draw) ==
  //                     GL::Framebuffer::Status::Complete);
  // ::Debug::error("TEST",_fbo_gbuffer.checkStatus(GL::FramebufferTarget::Draw) ==
  //                     GL::Framebuffer::Status::Complete);
#endif
}

void MainApplication::drawSettings() {
  ImGui::Begin("Render Options");
  const float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);

  if (ImGui::Button("Hot Reload Shaders")) {
    Utility::Resource::overrideGroup("ssao-data",
                                     "src/render/shaders/resources.conf");
    _ssaoShader         = SsaoShader{SSAO_SAMPLES};
    _ssaoApplyShader    = SsaoApplyShader{_ssaoApplyFlag};
    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _obsMeshShader      = ObsMeshShader{};
  }

  if (ImGui::CollapsingHeader("SSAO" /*, ImGuiTreeNodeFlags_DefaultOpen*/)) {
    ImGui::Indent();
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

    static bool drawOcclusion = false;
    if (ImGui::Checkbox("Show Occlusion Factor", &drawOcclusion)) {
      if (drawOcclusion)
        _ssaoApplyFlag = SsaoApplyShader::Flag::DrawAmbientOcclusion;
      else
        _ssaoApplyFlag = {};
      _ssaoApplyShader = SsaoApplyShader{_ssaoApplyFlag};
    }

    ImGui::Unindent();
  }

  if (ImGui::CollapsingHeader("Other", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();

    ImGui::TextBrowser(*_fileDialog.get(), _matcap_file, [&]() {
      loadTexture(_matcap_file, _matcap, GL::SamplerWrapping::ClampToEdge);
      _yarnGeometryShader.bindMatCap(_matcap);
    });
    static bool _rep = true;
    ImGui::TextBrowser(*_fileDialog.get(), _clothtexture_file, [&]() {
      loadTexture(_clothtexture_file, _clothTexture,
                  _rep ? GL::SamplerWrapping::Repeat : GL::SamplerWrapping::ClampToEdge);
      _yarnGeometryShader.bindClothTexture(_clothTexture);
    });
    ImGui::SameLine();
    ImGui::Checkbox("Rep", &_rep);

    ImGui::Checkbox("Yarns", &_render_yarns);
    ImGui::SameLine();
    ImGui::Checkbox("Mesh", &_render_mesh);
    ImGui::SameLine();
    ImGui::Checkbox("Obs", &_render_obstacles);
    ImGui::Checkbox("Rotate Scene", &_rotate_scene);

    {
      ImGui::PushItemWidth(100.0f);
      ImGui::DragFloat("yarn radius mult", &_render_radius_mult, 0.01f, 0.0f,
                       2.0f);
      static float _ts = _render_nmtwist  * 0.01f;
      if (ImGui::DragFloat("NM twist(1/cm)", &_ts, 0.1f, 0.0f,
                       100.0f))
                       _render_nmtwist = _ts * 100;
      ImGui::DragFloat("NM num", &_render_nmnum, 1.0f, 0.0f,
                       10.0f);
      ImGui::DragFloat("NM height", &_render_nmheight, 1.0f, 0.0f,
                       300.0f);
      // ImGui::DragFloat("mesh offset", &_mesh_dz, 0.001f, -1.0f, 1.0f);
      ImGui::DragFloat("tex scale", &_clothTexture_scale, 0.1f, 0.0f, 100.0f);
      ImGui::PopItemWidth();
    }

    ImGui::TextBrowser(*_fileDialog.get(), _matcapObs_file, [&]() {
      loadTexture(_matcapObs_file, _matcapObs,
                  GL::SamplerWrapping::ClampToEdge);
      _obsMeshShader.bindMatCap(_matcapObs);
    });

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
    ImGui::Unindent();
  }

  ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();

  ///

  ImGui::Begin("Sim Options");
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);
  ImGui::Checkbox("Pause", &_paused);
  ImGui::SameLine();
  if (ImGui::Button("[S]tep")) {
    _single_step = true;
    _paused      = true;
  }
  if (ImGui::DragInt("Min. Loop", &_min_loop_ms, 1, 4, 100, "%d (ms)"))
    setMinimalLoopPeriod(_min_loop_ms);
  ImGui::SameLine();
  if (ImGui::Button("(16)##loop16")) {
    _min_loop_ms = 16;
    setMinimalLoopPeriod(_min_loop_ms);
  }
  ImGui::Separator();
  ImGui::Checkbox("Repeat Mesh Frame", &_yarnMapperSettings.repeat_frame);
  ImGui::Checkbox("Flat Normals", &_yarnMapperSettings.flat_normals);
  ImGui::Checkbox("Flat Strains", &_yarnMapperSettings.flat_strains);
  ImGui::Checkbox("Def.SameTri", &_yarnMapperSettings.default_same_tri);
  // if (_yarnMapperSettings.flat_normals) {
  //   ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
  // }
  ImGui::Checkbox("Shepard Weights", &_yarnMapperSettings.shepard_weights);
  // if (_yarnMapperSettings.flat_normals)
  //   ImGui::PopStyleVar();
  ImGui::SliderFloat("Deform Ref.", &_yarnMapperSettings.deform_reference, 0.0f,
                     1.0f);
  ImGui::Checkbox("Shell Map", &_yarnMapperSettings.shell_map);
  ImGui::DragFloat("Phong", &_yarnMapperSettings.phong_deformation, 0.01f, 0.0f, 1.0f);
  ImGui::DragFloat("SVDClamp", &_yarnMapperSettings.svdclamp, 0.01f, 0.0f, 1.0f);
  ImGui::Checkbox("GPU", &_yarnMapperSettings.gpu_compute);
  // ImGui::Checkbox("DBG", &_yarnMapperSettings.debug_toggle);

  ImGui::TextBrowser(*_folderDialog.get(), _yarnMapperSettings.modelfolder);

  static int elem                = _yarnMapperSettings.provider_type;
  static constexpr auto count    = YarnMapper::Settings::Provider::COUNT;
  const char *elems_names[count] = {"ObjSeq", "BinSeq", "XPBD"};
  const char *elem_name =
      (elem >= 0 && elem < count) ? elems_names[elem] : "Unknown";
  if (ImGui::SliderInt("slider enum", &elem, 0, count - 1, elem_name))
    _yarnMapperSettings.provider_type =
        static_cast<YarnMapper::Settings::Provider>(elem);

  if (_yarnMapperSettings.provider_type ==
      YarnMapper::Settings::Provider::ObjSeq) {
    ImGui::TextBrowser(*_folderDialog.get(),
                       _yarnMapperSettings.objseq_settings.folder);
  } else if (_yarnMapperSettings.provider_type ==
             YarnMapper::Settings::Provider::BinSeq) {
    ImGui::TextBrowser(*_fileDialog.get(),
                       _yarnMapperSettings.binseq_settings.filepath);
  } else {
    ImGui::TextUnformatted("- not implemented -");
  }

  ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();

  ///

  ImGui::Begin("Stats");
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  // ImGui::PushItemWidth(100.0f);

  if (_profiler.isMeasurementAvailable(0)) {
    double fps = 1e9 / _profiler.frameTimeMean();
    ImGui::Text("FPS: %.1f", fps);
    ImGui::Separator();
    std::string stats = _profiler.statistics();
    ImGui::TextUnformatted(stats.c_str());
    ImGui::Separator();
  }
  {
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 4));
    ImGui::Columns(2, NULL, true);
    ImGui::SetColumnWidth(1, 100);
    auto averages          = _yarnMapper->m_timer.getAverageList();
    constexpr double scale = 0.001;
    for (const auto &avg : averages) ImGui::Text("%s:", avg.first.c_str());
    ImGui::NextColumn();
    for (const auto &avg : averages)
      ImGui::Text("%8.2f ms", scale * avg.second);
    ImGui::PopStyleVar();
    ImGui::Columns(1);
  }

  // ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();


  ImGui::Begin("DBG");
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);
  ImGui::Checkbox("D", &_yarnMapper->m_dbg.toggle);
  {
    ImGui::TextUnformatted("S"); ImGui::SameLine();
    ImGui::PushItemWidth(20.0f);
    ImGui::DragFloat("##S0", &_yarnMapper->m_dbg.strain_toggle[0], 0.01f, 0.0f, 1.0f); ImGui::SameLine();
    ImGui::DragFloat("##S1", &_yarnMapper->m_dbg.strain_toggle[1], 0.01f, 0.0f, 1.0f); ImGui::SameLine();
    ImGui::DragFloat("sx sa sy##S2", &_yarnMapper->m_dbg.strain_toggle[2], 0.01f, 0.0f, 1.0f);
    // ImGui::DragFloat("##S3", &_yarnMapper->m_dbg.strain_toggle[3], 0.01f, 0.0f, 1.0f); ImGui::SameLine();
    // ImGui::DragFloat("##S4", &_yarnMapper->m_dbg.strain_toggle[4], 0.01f, 0.0f, 1.0f); ImGui::SameLine();
    // ImGui::DragFloat("IIxx IIxy IIyy##S5", &_yarnMapper->m_dbg.strain_toggle[5], 0.01f, 0.0f, 1.0f);
    ImGui::DragFloat("IIxx IIxy IIyy##S5", &_yarnMapper->m_dbg.strain_toggle[3], 0.01f, 0.0f, 1.0f);
    ImGui::PopItemWidth();
  }
  #ifdef DO_DEBUG_STATS
  {
    // Treeview with imgui histograms
    if (ImGui::TreeNode("Histograms")) {
      static float hist_scale = 0.25f;
      ImGui::DragFloat("scale", &hist_scale, 0.05f, 0.1f, 1.0f);
      ImGui::Text("%.2f -- %.2f",double(_yarnMapper->m_dbg.hist_min), double(_yarnMapper->m_dbg.hist_max));
      // ImGui::PlotHistogram("sx", _yarnMapper->m_dbg.hist_counts[0].data(),  _yarnMapper->m_dbg.hist_nbins, 0, NULL, 0.0f, 1.0f, ImVec2(160.0f,80.0f));
      ImGui::PlotHistogram("##sx", _yarnMapper->m_dbg.hist_counts[0].data(),  _yarnMapper->m_dbg.hist_nbins, 0, "SX", 0.0f, float(_yarnMapper->m_dbg.hist_stepcount)*hist_scale, ImVec2(256.0f,64.0f)); 
      ImGui::PlotHistogram("##sa", _yarnMapper->m_dbg.hist_counts[1].data(),  _yarnMapper->m_dbg.hist_nbins, 0, "SA", 0.0f, float(_yarnMapper->m_dbg.hist_stepcount)*hist_scale, ImVec2(256.0f,64.0f));
      ImGui::PlotHistogram("##sy", _yarnMapper->m_dbg.hist_counts[2].data(),  _yarnMapper->m_dbg.hist_nbins, 0, "SY", 0.0f, float(_yarnMapper->m_dbg.hist_stepcount)*hist_scale, ImVec2(256.0f,64.0f));

      if (ImGui::Button("Reset")) {
        for (size_t i = 0; i < _yarnMapper->m_dbg.hist_counts.size(); i++)
        {
          _yarnMapper->m_dbg.hist_counts[i].clear();
          _yarnMapper->m_dbg.hist_counts[i].resize(_yarnMapper->m_dbg.hist_nbins, 0);
          _yarnMapper->m_dbg.hist_stepcount = 0;
        }
      }
    }
  }
  #endif

  ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();
}

void MainApplication::keyPressEvent(KeyEvent &event) {
  if (_imgui.handleKeyPressEvent(event))
    return;

  static bool MS = true;
  bool success;

  // key press events, that are only allowed while not editing a text field
  if (!isTextInputActive()) {
    float cam_d = 2.0f;
    if ((event.modifiers() & KeyEvent::Modifier::Shift))
      cam_d *= 2.0f;
    if ((event.modifiers() & KeyEvent::Modifier::Ctrl))
      cam_d *= 2.0f;
    if ((event.modifiers() & KeyEvent::Modifier::Alt))
      cam_d *= 0.5f;
    switch (event.key()) {
      case KeyEvent::Key::Esc:
        this->exit();
        break;
      case KeyEvent::Key::E:
        ::Debug::log("Exporting FBX...");
        success = _yarnMapper->export2fbx("test.fbx"); // TODO filename either by date or incremented number
        ::Debug::log("Result:", success);
        break;
      case KeyEvent::Key::Space:
        _paused = !_paused;
        break;
      case KeyEvent::Key::R:
        reset_simulation();
        break;
      case KeyEvent::Key::S:
        if ((event.modifiers() & KeyEvent::Modifier::Shift)) {
          _paused                          = false;
          _yarnMapperSettings.repeat_frame = !_yarnMapperSettings.repeat_frame;
        } else {
          _single_step = true;
          _paused      = true;
        }
        break;
      case KeyEvent::Key::B:
      if (MS)
        GL::Renderer::disable(GL::Renderer::Feature::Multisampling);
      else
        GL::Renderer::enable(GL::Renderer::Feature::Multisampling);
        MS = !MS;
        break;
      case KeyEvent::Key::M:
        _render_mesh = !_render_mesh;
        break;
      case KeyEvent::Key::O:
        _render_obstacles = !_render_obstacles;
        break;
      case KeyEvent::Key::Y:
        _render_yarns = !_render_yarns;
        break;
      case KeyEvent::Key::F:
        _rotate_scene = !_rotate_scene;
        break;
      // case KeyEvent::Key::T:
      //   _yarnMapperSettings.debug_toggle = !_yarnMapperSettings.debug_toggle;
      //   break;
      case KeyEvent::Key::G:
        _yarnMapperSettings.gpu_compute = !_yarnMapperSettings.gpu_compute;
        break;
      case KeyEvent::Key::NumOne:
      case KeyEvent::Key::One:
        _arcball->setViewParameters(cam_d * Vector3(+0.4f, 0.3f, 0.5f),
                                    Vector3(0), Vector3(0, 1, 0));
        break;
      case KeyEvent::Key::NumZero:
      case KeyEvent::Key::Zero:
      case KeyEvent::Key::NumTwo:
      case KeyEvent::Key::Two:
        _arcball->setViewParameters(cam_d * Vector3(0, 1, 0), Vector3(0),
                                    Vector3(0, 0, -1));
        break;
      case KeyEvent::Key::NumThree:
      case KeyEvent::Key::Three:
        _arcball->setViewParameters(cam_d * Vector3(0, 0, 1), Vector3(0),
                                    Vector3(0, 1, 0));
        break;
      default:
        break;
    }
  }
}

void MainApplication::keyReleaseEvent(KeyEvent &event) {
  if (_imgui.handleKeyReleaseEvent(event))
    return;
}

void MainApplication::textInputEvent(TextInputEvent &event) {
  if (_imgui.handleTextInputEvent(event))
    return;
}

void MainApplication::mousePressEvent(MouseEvent &event) {
  if (_imgui.handleMousePressEvent(event))
    return;
  /* Enable mouse capture so the mouse can drag outside of the window */
  /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
  SDL_CaptureMouse(SDL_TRUE);
  _arcball->initTransformation(event.position());
  event.setAccepted();
  redraw();
}

void MainApplication::mouseReleaseEvent(MouseEvent &event) {
  if (_imgui.handleMouseReleaseEvent(event))
    return;
  /* Disable mouse capture again */
  /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
  SDL_CaptureMouse(SDL_FALSE);
}

void MainApplication::mouseMoveEvent(MouseMoveEvent &event) {
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

void MainApplication::mouseScrollEvent(MouseScrollEvent &event) {
  if (_imgui.handleMouseScrollEvent(event)) {
    /* Prevent scrolling the page */
    event.setAccepted();
    return;
  }

  Float delta = event.offset().y();

  if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
    delta *= 0.5f;
  if (Math::abs(delta) < 1.0e-2f)
    return;

  _arcball->zoom(delta * 0.1f);

  event.setAccepted();
  redraw();
}

MAGNUM_APPLICATION_MAIN(Magnum::MainApplication)

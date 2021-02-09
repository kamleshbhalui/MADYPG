#include "MainApplication.h"

#include <imgui_stdlib.h>

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
    const std::string &label, ImGui::FileBrowser &browser, std::string &txt,
    std::function<void()> onChange = []() {});
bool TextBrowser(const std::string &label, ImGui::FileBrowser &browser,
                 std::string &txt, std::function<void()> onChange) {
  ImGui::PushID(&txt);
  ImGui::PushItemWidth(150);
  if (ImGui::InputText(label.c_str(), &txt,
                       ImGuiInputTextFlags_EnterReturnsTrue)) {
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
  // create new YarnMapper instance
  if (!_yarnMapper)
    _yarnMapper = std::make_unique<YarnMapper>();
  else
    _yarnMapper.reset(new YarnMapper());

  // copy settings
  _yarnMapper->m_settings = _yarnMapperSettings;

  // init
  _yarnMapper->initialize();
  
  // set up drawable objects for rendering of yarns, mesh, obstacles
  _yarnDrawable.clear();
  _yarnDrawable.emplace_back(_yarnGeometryShader,
                             _yarnMapper->getVertexBuffer(),
                             _yarnMapper->getIndexBuffer());

  auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
  _meshdrawable.release();
  _meshdrawable =
      std::make_unique<MeshDrawable<MeshShader>>(_meshShader, mesh.F, mesh.X);

  _obsmeshdrawables.clear();
  for (const auto &obs : _yarnMapper->getMeshSimulation()->getObstacles()) {
    _obsmeshdrawables.emplace_back(_obsMeshShader);
    _obsmeshdrawables.back().setIndices(obs.mesh.F);
    _obsmeshdrawables.back().setVertices(obs.mesh.X);
  }

  _frame = 0;
}

MainApplication::MainApplication(const Arguments &arguments)
    : Platform::Application{arguments, NoCreate} {
#ifdef MSAA
  ::Debug::logf("Using %dx MSAA\n", MSAA);
  GL::Renderer::enable(GL::Renderer::Feature::Multisampling);
#endif
#ifdef SUPERSAMPLING
  ::Debug::logf("Using %dx SSAA\n", SUPERSAMPLING);
#endif

  /* Setup window */
  {
    const Vector2 dpiScaling = this->dpiScaling({});
    Configuration conf;
    // conf.setSize({1200, 800}, dpiScaling);
    // conf.setSize({1920, 1080}, dpiScaling);
    conf.setSize({1600, 1000}, dpiScaling);
    conf.setTitle(
            "Mechanics-Aware Deformation of Yarn Pattern Geometry [Sperl et "
            "al. 2021]")
        .setSize(conf.size(), dpiScaling)
        .setWindowFlags(Configuration::WindowFlag::Maximized |
                        // Configuration::WindowFlag::Borderless |
                        Configuration::WindowFlag::Resizable);
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

  /* SSAO and other setup */
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

    // #ifdef SUPERSAMPLING
    // GL::defaultFramebuffer.setViewport({{}, vpSize});
    //     _fbo_gbuffer.setViewport({{}, vpSize * SUPERSAMPLING});
    //     _fbo_ssao.setViewport({{}, vpSize * SUPERSAMPLING});
    // #ifdef MSAA
    //     _fbo_ssao = GL::Framebuffer{viewport};
    // #endif
    // setupFramebuffer(vpSize * SUPERSAMPLING);
    // #else
    setupFramebuffer(vpSize);
    // #endif

    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _ground             = GroundShader{};
    _obsMeshShader      = ObsMeshShader{};
    _ssaoShader         = SsaoShader{SSAO_SAMPLES};
    _ssaoApplyShader    = SsaoApplyShader{};

    {  // initial yarnmapper settings overrides
      _yarnMapperSettings.modelfolder = "data/yarnmodels/model_stock_9";
      _yarnMapperSettings.provider_type =
          YarnMapper::Settings::Provider::BinSeq;
      _yarnMapperSettings.binseq_settings.filepath =
          "data/binseqs/sxsy_const.bin";
      _yarnMapperSettings.objseq_settings.folder = "data/objseqs/sxsy_const";
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
    loadTexture(_matcapMesh_file, _matcapMesh,
                GL::SamplerWrapping::ClampToEdge);
    loadTexture(_clothtexture_file, _clothTexture, GL::SamplerWrapping::Repeat);
    loadTexture(_gridtexture_file, _gridTexture, GL::SamplerWrapping::Repeat);
    // loadTexture1D(_normalMap_file, _normalMap, GL::SamplerWrapping::Repeat);
    loadTexture(_normalMap_file, _normalMap, GL::SamplerWrapping::Repeat);
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

  setSwapInterval(0);  // 0 no vsync, 1 vsync
  setMinimalLoopPeriod(_min_loop_ms); // limit fps (default: ~60)
}

void MainApplication::drawEvent() {
  _profiler.beginFrame();

  // update simulation
  bool simChanged = false;
  if (!_paused || _single_step) {  // SIM
    if (_yarnMapper->isInitialized()) {
      // force update some settings
      _yarnMapper->m_settings = _yarnMapperSettings;
      
      // actual update
      _yarnMapper->step();

      const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
      // if (_yarnMapper->getMeshSimulation()->meshIndicesDirty())
      _meshdrawable->updateIndexCount(mesh.F.getGPUSize() * 3);
    }
    _single_step = false;
    simChanged   = true;

    if (_pauseAt == _frame)
      _paused = true;

    if (!_yarnMapperSettings.repeat_frame)
      ++_frame;
  }

  const bool camChanged = _arcball->updateTransformation();

  GL::defaultFramebuffer.clear(GL::FramebufferClear::Color |
                               GL::FramebufferClear::Depth);
  bool require_redraw = simChanged || camChanged;
  if (require_redraw) {
    Matrix4 tf = _arcball->viewMatrix();

    /* render the scene into g-buffer */
    // ie. render colors, positions, normals as textures to use for SSAO
    _fbo_gbuffer
        .mapForDraw(
            {{YarnShader::AlbedoOutput, GL::Framebuffer::ColorAttachment{0}},
             {YarnShader::PositionsOutput, GL::Framebuffer::ColorAttachment{1}},
             {YarnShader::NormalsOutput, GL::Framebuffer::ColorAttachment{2}}})
        // .clear(GL::FramebufferClear::Depth | GL::FramebufferClear::Color)
        .clearColor(0, _bgColor)
        .clearColor(
            1, Color4(0.0, 0.0, -1000.0,
                      1.0))  // NOTE: using arbitrary -1000 as background depth
        .clearColor(2, Color4(0.0, 0.0, 1.0, 1.0))
        .clearDepth(1.0)
        .bind();

    if (_render_ground) {  // NOTE: before rotating rest of scene
      _ground.bindTexture(_gridTexture);
      _ground.renderQuad(tf, _projection);
    }

    if (_rotate_scene) {
      tf = tf * Matrix4::rotation(Rad(-1.57079632679f), Vector3::xAxis());
    }

    if (_render_yarns && _yarnMapper->isInitialized()) {
      _yarnGeometryShader.bindMatCap(_matcap);
      _yarnGeometryShader.bindClothTexture(_clothTexture);
      _yarnGeometryShader.bindNormalMap(_normalMap);
      _yarnGeometryShader.setProjection(_projection);
      _yarnGeometryShader.setTextureScale(_clothUV_scale);
      _yarnGeometryShader.setTextureOffset(_clothUV_offset);
      _yarnDrawable.back().m_radius =
          _yarnMapper->getRadius() * _render_radius_mult;
      _yarnDrawable.back().m_nmtwist  = _render_nmtwist;
      _yarnDrawable.back().m_nmnum    = _render_nmnum;
      _yarnDrawable.back().m_nmheight = _render_nmheight;
      _yarnDrawable.back().m_nmlen    = _render_nmlen;
      for (auto &line : _yarnDrawable) line.draw(tf);
    }

    if (_render_mesh) {
      // for debug rendering only: overwrite mesh shader data with material
      // space
      auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
      if (!_yarnMapperSettings.shell_map && _render_mesh) {
        // m_mesh.F.
        // draw instead: Fms & U+0
        std::vector<Mesh::WSVertex> tmpX;
        tmpX.reserve(mesh.U.getCPUSize());
        for (const auto &u : mesh.U.cpu()) tmpX.push_back({u.u, u.v, 0});
        std::vector<Mesh::Face> tmpF;
        tmpF.reserve(mesh.Fms.getCPUSize());
        for (const auto &f : mesh.Fms.cpu()) tmpF.push_back(f);
        mesh.X.gpu().setData({tmpX.data(), uint32_t(tmpX.size())});
        mesh.F.gpu().setData({tmpF.data(), uint32_t(tmpF.size())});

        _meshdrawable->updateIndexCount(mesh.Fms.getGPUSize() * 3);
      }

      _meshShader.bindMatCap(_matcapMesh);
      _meshShader.setProjection(_projection);
      GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
      _meshdrawable->draw(tf);
      GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

      // undo debug render
      if (!_yarnMapperSettings.shell_map && _render_mesh) {
        mesh.F.bufferData();
      }
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

    // Compute SSAO occlusion from positions and normals
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
        .bind();

    _ssaoShader.draw(_screenAlignedTriangle);
#else
    _fbo_gbuffer
        .mapForDraw({{SsaoShader::AmbientOcclusionOutput,
                      GL::Framebuffer::ColorAttachment{3}}})
        .clear(GL::FramebufferClear::Color);
    _ssaoShader.draw(_screenAlignedTriangle);
#endif
  }

  // combine color and occlusion for final image
  GL::defaultFramebuffer.bind();
  _ssaoApplyShader.bindAlbedoTexture(_albedo)
      .bindOcclusionTexture(_occlusion)
      .bindPositionTexture(_positions)
      .setAOBlurRadius(_ao_blur_radius)
      .setAOBlurFeature(_ao_blur_feature)
      .setAOPow(_ao_pow)
      .draw(_screenAlignedTriangle);

  // draw GUI on top
  _imgui.newFrame();
  if (ImGui::GetIO().WantTextInput && !isTextInputActive())
    startTextInput();
  else if (!ImGui::GetIO().WantTextInput && isTextInputActive())
    stopTextInput();

  if (_simple_gui) {
    drawGUISimple();
  } else {
    drawGUIStats();
    drawGUISettings();
    drawGUIRender();
  }
  drawGUISliders();

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

#ifdef SUPERSAMPLING
  _fbo_gbuffer.setViewport({{}, fbSize * SUPERSAMPLING});
#ifdef MSAA
  _fbo_ssao.setViewport({{}, fbSize * SUPERSAMPLING});
#endif
#else
  _fbo_gbuffer.setViewport({{}, fbSize});
#ifdef MSAA
  _fbo_ssao.setViewport({{}, fbSize});
#endif
#endif

  _arcball->reshape(wSize);
  _projection = Matrix4::perspectiveProjection(
      _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near,
      _proj_far);

#ifdef SUPERSAMPLING
  setupFramebuffer(fbSize * SUPERSAMPLING);
#else
  setupFramebuffer(fbSize);
#endif

  _imgui.relayout(Vector2{wSize} / event.dpiScaling(), event.windowSize(),
                  fbSize);
}

void MainApplication::setupFramebuffer(const Vector2i &size) {
  setupTexture(_albedo, size, GL::TextureFormat::RGBA32F);
  setupTexture(_positions, size, GL::TextureFormat::RGBA32F);
  setupTexture(_normals, size, GL::TextureFormat::RGBA32F);
  setupTexture(_occlusion, size, GL::TextureFormat::R32F);
  setupTexture(_depth, size, GL::TextureFormat::DepthComponent32F);

#ifdef MSAA
  _fbo_gbuffer.attachTexture(GL::Framebuffer::BufferAttachment::Depth, _depth)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _albedo)
      .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions)
      .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals);
#else
  _fbo_gbuffer
      .attachTexture(GL::Framebuffer::BufferAttachment::Depth, _depth, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _albedo, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals, 0)
      .attachTexture(GL::Framebuffer::ColorAttachment{3}, _occlusion, 0);
#endif
  CORRADE_INTERNAL_ASSERT(
      _fbo_gbuffer.checkStatus(GL::FramebufferTarget::Draw) ==
      GL::Framebuffer::Status::Complete);

#ifdef MSAA
  _fbo_ssao
      // .attachTexture(GL::Framebuffer::ColorAttachment{1}, _positions)
      // .attachTexture(GL::Framebuffer::ColorAttachment{2}, _normals)
      .attachTexture(GL::Framebuffer::ColorAttachment{0}, _occlusion,
                     0); 
  CORRADE_INTERNAL_ASSERT(
      _fbo_gbuffer.checkStatus(GL::FramebufferTarget::Read) ==
      GL::Framebuffer::Status::Complete);
  CORRADE_INTERNAL_ASSERT(_fbo_ssao.checkStatus(GL::FramebufferTarget::Draw) ==
                          GL::Framebuffer::Status::Complete);
  CORRADE_INTERNAL_ASSERT(_fbo_ssao.checkStatus(GL::FramebufferTarget::Read) ==
                          GL::Framebuffer::Status::Complete);
#endif
}

void MainApplication::drawGUIStats() {
  ImGui::Begin("Stats/Playback");
  constexpr float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  // ImGui::PushItemWidth(100.0f);

  if (_profiler.isMeasurementAvailable(0)) {
    double fps = 1e9 / _profiler.frameTimeMean();
    ImGui::Text("FPS:    %5.0f", fps);
    // std::string stats = _profiler.statistics();
    // ImGui::TextUnformatted(stats.c_str());
    // ImGui::Separator();
  }
  ImGui::TextUnformatted("Min. Frame Time:");
  ImGui::SameLine();
  ImGui::PushItemWidth(60.0f);
  if (ImGui::DragInt("##Min. Frame", &_min_loop_ms, 1, 4, 100, "%d ms"))
    setMinimalLoopPeriod(_min_loop_ms);
  ImGui::PopItemWidth();

  ImGui::Text(
      "# vertices: %s",
      ::Debug::format_locale(_yarnMapper->getNumVertices(), "en_US.UTF-8")
          .c_str());
  {
    static const std::vector<std::string> labels{
        "mesh: update", "mesh: strains", "yarns: deform", "yarns: map"};
    const auto &timer = _yarnMapper->m_timer;

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 4));
    ImGui::Columns(2, NULL, true);
    ImGui::SetColumnWidth(1, 100);
    for (const auto &label : labels) {
      ImGui::TextUnformatted(label.c_str());
    }
    ImGui::NextColumn();
    for (const auto &label : labels) {
      ImGui::Text("%8.2f ms", 0.001 * timer.getAverage(label));
    }
    ImGui::PopStyleVar();
    ImGui::Columns(1);
  }
  if (ImGui::TreeNode("All Timers")) {
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 4));
    ImGui::Columns(2, NULL, true);
    ImGui::SetColumnWidth(1, 100);
    auto averages          = _yarnMapper->m_timer.getAverageList();
    constexpr double scale = 0.001;
    for (const auto &avg : averages) ImGui::Text("%s", avg.first.c_str());
    ImGui::NextColumn();
    for (const auto &avg : averages)
      ImGui::Text("%8.2f ms", scale * avg.second);
    ImGui::PopStyleVar();
    ImGui::Columns(1);
    ImGui::TreePop();
  }

  ImGui::Separator();  // Playback

  ImGui::Text("Frame: %d", _frame);
  ImGui::Checkbox("[P]ause", &_paused);
  ImGui::SameLine();
  if (ImGui::Button("[S]tep")) {
    _single_step = true;
    _paused      = true;
  }
  ImGui::SameLine();
  if (ImGui::Button("[R]eset")) {
    reset_simulation();
  }
  // ImGui::SameLine();
  // ImGui::PushItemWidth(30.0f);
  // ImGui::DragInt("pause@", &_pauseAt, 1, -1, 1000);
  // ImGui::PopItemWidth();
  ImGui::Checkbox("Repeat Frame [Space]", &_yarnMapperSettings.repeat_frame);

  // ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();
}

void MainApplication::drawGUISettings() {
  ImGui::Begin("Animation Settings");
  constexpr float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);

  // ImGui::Checkbox("Flat Normals", &_yarnMapperSettings.flat_normals);
  // ImGui::Checkbox("Flat Strains", &_yarnMapperSettings.flat_strains);
  // ImGui::Checkbox("Def.SameTri", &_yarnMapperSettings.default_same_tri);
  // if (_yarnMapperSettings.flat_normals) {
  //   ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
  // }
  // ImGui::Checkbox("Shepard Weights", &_yarnMapperSettings.shepard_weights);
  // if (_yarnMapperSettings.flat_normals)
  //   ImGui::PopStyleVar();
  ImGui::DragFloat("Min. length", &_yarnMapperSettings.min_yarn_length_per_r,
                   0.1f, 0.0f, 100.0f, "%.1f R");
  ImGui::SliderFloat("Deform Ref.", &_yarnMapperSettings.deform_reference, 0.0f,
                     1.0f);
  ImGui::SliderFloat("Lin. Bending", &_yarnMapperSettings.linearized_bending,
                     0.0f, 1.0f);
  ImGui::Checkbox("Shell Map", &_yarnMapperSettings.shell_map);
  ImGui::DragFloat("Phong", &_yarnMapperSettings.phong_deformation, 0.01f, 0.0f,
                   1.0f);
  ImGui::DragFloat("SVDClamp", &_yarnMapperSettings.svdclamp, 0.01f, 0.0f,
                   1.0f);
  ImGui::Checkbox("GPU", &_yarnMapperSettings.gpu_compute);
  // ImGui::Checkbox("DBG", &_yarnMapperSettings.debug_toggle);

  ImGui::TextBrowser("##txt", *_folderDialog.get(),
                     _yarnMapperSettings.modelfolder);

  static int elem                = _yarnMapperSettings.provider_type;
  static constexpr auto count    = YarnMapper::Settings::Provider::COUNT;
  const char *elems_names[count] = {"ObjSeq", "BinSeq", "PBD"};
  const char *elem_name =
      (elem >= 0 && elem < count) ? elems_names[elem] : "Unknown";
  if (ImGui::SliderInt("Input Method", &elem, 0, count - 1, elem_name))
    _yarnMapperSettings.provider_type =
        static_cast<YarnMapper::Settings::Provider>(elem);

  ImGui::Indent();
  if (_yarnMapperSettings.provider_type ==
      YarnMapper::Settings::Provider::ObjSeq) {
    ImGui::TextBrowser("##txt", *_folderDialog.get(),
                       _yarnMapperSettings.objseq_settings.folder);
    ImGui::Checkbox(
        "Const UV",
        &_yarnMapperSettings.objseq_settings.constant_material_space);
    ImGui::DragFloat("scale", &_yarnMapperSettings.objseq_settings.scale, 0.01f,
                     0.01f, 10.0f);
  } else if (_yarnMapperSettings.provider_type ==
             YarnMapper::Settings::Provider::BinSeq) {
    ImGui::TextBrowser("##txt", *_fileDialog.get(),
                       _yarnMapperSettings.binseq_settings.filepath);
    ImGui::DragFloat("scale", &_yarnMapperSettings.binseq_settings.scale, 0.01f,
                     0.01f, 10.0f);
  } else {
    ImGui::SliderInt("MethodS",
                     &_yarnMapperSettings.pbd_settings.simulationMethod, 0, 4);
    ImGui::SliderInt("MethodB", &_yarnMapperSettings.pbd_settings.bendingMethod,
                     0, 2);
    ImGui::DragFloat3("Stiffness", _yarnMapperSettings.pbd_settings.stiffness,
                      0.001f, 0.0f, 10.0f);
    ImGui::DragFloat("Bending",
                     &_yarnMapperSettings.pbd_settings.bending_stiffness,
                     0.001f, 0.0f, 10.0f);
    ImGui::DragInt("Iterations", &_yarnMapperSettings.pbd_settings.iterations,
                   1, 1, 100);
    ImGui::SliderInt("Substeps", &_yarnMapperSettings.pbd_settings.substeps, 1,
                     20);
    ImGui::DragFloat("dt", &_yarnMapperSettings.pbd_settings.timestep, 0.001f,
                     0.0001f, 1.0f);
    ImGui::DragFloat("Density", &_yarnMapperSettings.pbd_settings.density,
                     0.01f, 0.0f, 10.0f);
    ImGui::DragFloat3("Gravity", _yarnMapperSettings.pbd_settings.gravity,
                      0.01f, -10.0f, 10.0f);
  }
  ImGui::Unindent();

  ImGui::PopItemWidth();
  ImGui::PopStyleVar();
  ImGui::End();
}

void MainApplication::drawGUIRender() {
  ImGui::Begin("Render Options");
  constexpr float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);

  ImGui::TextUnformatted("Simple GUI: [H]");

  if (ImGui::Button("Hot Reload Shaders")) {
    Utility::Resource::overrideGroup("ssao-data",
                                     "src/render/shaders/resources.conf");
    _ssaoShader         = SsaoShader{SSAO_SAMPLES};
    _ssaoApplyShader    = SsaoApplyShader{_ssaoApplyFlag};
    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _ground             = GroundShader{};
    _obsMeshShader      = ObsMeshShader{};

    Utility::Resource::overrideGroup("compute-shaders",
                                     "src/yarns/shaders/resources.conf");
    _yarnMapper->reloadShaders();
  }

  if (ImGui::CollapsingHeader("Camera")) {
    ImGui::Indent();
    float params[2] = {_proj_near, _proj_far};
    if (ImGui::DragFloat2("near/far", params, 0.01f, 0.001f, 20.0f)) {
      _proj_near  = params[0];
      _proj_far   = params[1];
      _projection = Matrix4::perspectiveProjection(
          _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near,
          _proj_far);
    }

    if (ImGui::Button("Reset Camera"))
      _arcball->reset();
    ImGui::SameLine();
    static bool lagging = true;
    if (ImGui::Checkbox("Camera Lagging", &lagging)) {
      _arcball->setLagging(float(lagging) * 0.75f);
    }

    ImGui::TextUnformatted("Views: [1] [2] [3] (+ modifiers)");

    ImGui::Unindent();
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
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    ImGui::SliderInt("blur radius", &_ao_blur_radius, 0, 5);
    //  (ImGui::DragFloat("blur feature (cm)", &_ao_blur_feature, 0.01f,
    //  0.0f,100.0f,"%.2e"));
    {
      float val = _ao_blur_feature * 0.01f;
      if (ImGui::SliderFloat("blur feature (1/cm)", &val, 0.0f, 50.0f))
        _ao_blur_feature = val * 100.0f;
    }
    ImGui::PopStyleVar();
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
  if (ImGui::CollapsingHeader("Yarns", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();
    ImGui::Checkbox("Render [Y]arns", &_render_yarns);

    ImGui::DragFloat("yarn radius mult", &_render_radius_mult, 0.01f, 0.0f,
                     2.0f);

    ImGui::TextBrowser("MatCap ##Y", *_fileDialog.get(), _matcap_file, [&]() {
      loadTexture(_matcap_file, _matcap, GL::SamplerWrapping::ClampToEdge);
      _yarnGeometryShader.bindMatCap(_matcap);
    });
    static bool _rep = true;
    ImGui::TextBrowser("Texture##Y", *_fileDialog.get(), _clothtexture_file,
                       [&]() {
                         loadTexture(_clothtexture_file, _clothTexture,
                                     _rep ? GL::SamplerWrapping::Repeat
                                          : GL::SamplerWrapping::ClampToEdge);
                         _yarnGeometryShader.bindClothTexture(_clothTexture);
                       });
    ImGui::SameLine();
    ImGui::Checkbox("Rep.", &_rep);
    ImGui::DragFloat("uv scale", &_clothUV_scale, 0.1f, 0.0f, 100.0f);
    ImGui::DragFloat2("uv offset", _clothUV_offset.data(), 0.1f, -10.0f, 10.0f);

    if (ImGui::TreeNode("Ply/Fiber Texture")) {
      ImGui::DragFloat("twist(*r)", &_render_nmtwist, 0.1f, 0.0f, 10.0f);
      ImGui::DragFloat("# ply", &_render_nmnum, 1.0f, 0.0f, 10.0f);
      ImGui::DragFloat("height (*R)", &_render_nmheight, 0.01f, 0.0f, 2.0f);
      ImGui::DragFloat("length", &_render_nmlen, 0.1f, 0.0f, 20.0f);
      ImGui::TreePop();
    }

    ImGui::Unindent();
  }

  if (ImGui::CollapsingHeader("Obstacle", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();
    ImGui::Checkbox("Render [O]bstacle", &_render_obstacles);

    ImGui::TextBrowser("MatCap##O", *_fileDialog.get(), _matcapObs_file, [&]() {
      loadTexture(_matcapObs_file, _matcapObs,
                  GL::SamplerWrapping::ClampToEdge);
      _obsMeshShader.bindMatCap(_matcapObs);
    });

    ImGui::Unindent();
  }
  if (ImGui::CollapsingHeader("Cloth Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();
    ImGui::Checkbox("Render Cloth [M]esh", &_render_mesh);

    ImGui::TextBrowser("MatCap##M", *_fileDialog.get(), _matcapMesh_file,
                       [&]() {
                         loadTexture(_matcapMesh_file, _matcapMesh,
                                     GL::SamplerWrapping::ClampToEdge);
                         _meshShader.bindMatCap(_matcapMesh);
                       });

    ImGui::Unindent();
  }

  if (ImGui::CollapsingHeader("Other", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();

    {
      ImGui::Checkbox("Render [G]round", &_render_ground);
      ImGui::DragFloat("height", &_ground.dY, 0.01f, -10.0f, 10.0f);
      ImGui::DragFloat("scale", &_ground.scale, 0.01f, 0.001f, 1000.0f);
    }

    ImGui::Checkbox("[F]lip Scene", &_rotate_scene);

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
}

void MainApplication::drawGUISimple() {
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 10));
  ImGui::PushStyleColor(ImGuiCol_ResizeGrip, ImVec4(1, 1, 1, 0));
  ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0, 0, 0, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(1, 1, 1, 0.6));
  ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1, 1, 1, 0));

  static const ImGuiWindowFlags window_flags =
      ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_AlwaysAutoResize;
  static bool open_ptr = true;

  ImGui::Begin("##simple", &open_ptr, window_flags);
  ImGui::SetWindowFontScale(1.5);

  ImGui::Text(
      "# vertices: %s",
      ::Debug::format_locale(_yarnMapper->getNumVertices(), "en_US.UTF-8")
          .c_str());

  static const std::vector<std::string> labels{"mesh: strains", "yarns: deform",
                                               "yarns: map"};
  const auto &timer          = _yarnMapper->m_timer;
  double total_ms_meshupdate = 0.001 * timer.getAverage("mesh: update");
  double total_ms_ours       = 0;
  for (const auto &label : labels) {
    total_ms_ours += 0.001 * timer.getAverage(label);
  }

  ImGui::Text("Yarn Animation: %5.2f ms / frame", total_ms_ours);

  if (_yarnMapperSettings.provider_type ==
      YarnMapper::Settings::Provider::PBD) {
    ImGui::Text("PBD:            %5.2f ms", total_ms_meshupdate);
  }
  ImGui::TextUnformatted("Advanced GUI: [H]");
  ImGui::End();

  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();  // grab triangle
  ImGui::PopStyleVar();
}

void MainApplication::drawGUISliders() {
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 10));
  ImGui::PushStyleColor(ImGuiCol_ResizeGrip, ImVec4(1, 1, 1, 0));
  ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0, 0, 0, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(0.9, 0.9, 0.9, 1));
  ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(1, 1, 1, 0.6));
  ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1, 1, 1, 0));

  static const ImGuiWindowFlags window_flags =
      ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_AlwaysAutoResize;
  static bool open_ptr = true;

  static const bool slider_deform = true;
  static const bool slider_phong  = false;
  static const bool slider_clamp  = false;

  ImGui::Begin("##interact_force", &open_ptr, window_flags);
  ImGui::SetWindowFontScale(1.5);
  ImGui::PushItemWidth(500.0f);
  if (_yarnMapperSettings.provider_type ==
      YarnMapper::Settings::Provider::PBD) {
    static float magn = 0;
    ImGui::TextUnformatted("     ");
    ImGui::SameLine();
    ImGui::SliderFloat("##force", &magn, 0.0f, 3.0f, "force: x%.2f");
    if (magn > 0.0001f && !_paused && !_yarnMapperSettings.repeat_frame)
      _yarnMapper->applyForce(magn * 3, magn * 0.5f, magn * 0);
  }
  if (slider_deform) {
    ImGui::TextUnformatted("naive");
    ImGui::SameLine();
    // static bool _use_ours = true;
    if (ImGui::SliderFloat("ours##sep", &_yarnMapperSettings.deform_reference,
                           0.0f, 1.0f, "")) {
      // _use_ours = _yarnMapperSettings.deform_reference > 0.001f;
    }
  }
  if (slider_phong) {
    ImGui::TextUnformatted("     ");
    ImGui::SameLine();
    ImGui::SliderFloat("##phong", &_yarnMapperSettings.phong_deformation, 0.0f,
                       1.0f, "alpha = %.2f");
  }
  if (slider_clamp) {
    ImGui::TextUnformatted("     ");
    ImGui::SameLine();
    ImGui::SliderFloat("##clamp", &_yarnMapperSettings.svdclamp, 0.0f, 1.0f,
                       "lambda_min = %.2f");
  }
  ImGui::PopItemWidth();
  ImGui::End();

  // ImGui::Begin("##nsamples", &open_ptr, window_flags);
  // ImGui::SetWindowFontScale(1.5);
  // static int nsamples = 5;
  // bool changed = false;
  // changed |= ImGui::RadioButton("    5 x 5 x 5", &nsamples, 5);
  // changed |= ImGui::RadioButton("    9 x 9 x 9", &nsamples, 9);
  // changed |= ImGui::RadioButton(" 15 x 15 x 15", &nsamples, 15);
  // changed |= ImGui::RadioButton(" 31 x 31 x 31", &nsamples, 31);
  // if (changed)
  //   _yarnMapperSettings.modelfolder =
  //   "data/yarnmodels/num_samples/model_stock_" + std::to_string(nsamples);
  // ImGui::End();

  // ImGui::Begin("##bend", &open_ptr, window_flags);
  // ImGui::SetWindowFontScale(1.5);
  // int toggle = int(_yarnMapper->m_dbg.toggle);
  // bool changed = false;
  // changed |= ImGui::RadioButton("linearized bending", &toggle, 0);
  // changed |= ImGui::RadioButton("        4D bending", &toggle, 1);
  // if (changed)
  //   _yarnMapper->m_dbg.toggle = (toggle > 0);
  // ImGui::End();

  // ImGui::Begin("##slide", &open_ptr, window_flags);
  // ImGui::SetWindowFontScale(1.5);
  // bool toggle = _yarnMapper->m_dbg.toggle;
  // ImGui::Checkbox("Sliding Constraint", &toggle);
  // ImGui::End();

  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();
  ImGui::PopStyleColor();  // grab triangle
  ImGui::PopStyleVar();
}

void MainApplication::keyPressEvent(KeyEvent &event) {
  if (_imgui.handleKeyPressEvent(event))
    return;

  bool success;

  // key press events, that are only allowed while not editing a text field
  if (!isTextInputActive()) {
    float force_mult = 1.0f;
    if ((event.modifiers() & KeyEvent::Modifier::Shift))
      force_mult *= 10;
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
      case KeyEvent::Key::H:
        _simple_gui = !_simple_gui;
        break;
      // case KeyEvent::Key::D:
      //   _yarnMapper->m_dbg.toggle = !_yarnMapper->m_dbg.toggle;
      //   break;
      case KeyEvent::Key::E:
        ::Debug::log("Exporting FBX...");
        if ((event.modifiers() & KeyEvent::Modifier::Shift))
          success = _yarnMapper->export2fbx_cloth(
              "cloth.fbx");  // TODO filename either by date or incremented
                             // number
        else
          success = _yarnMapper->export2fbx(
              "yarns.fbx");  // TODO filename either by date or incremented
                             // number
        ::Debug::log("Result:", success);
        break;
      case KeyEvent::Key::P:
        _paused = !_paused;
        break;
      case KeyEvent::Key::Space:
        _yarnMapperSettings.repeat_frame = !_yarnMapperSettings.repeat_frame;
        break;
      case KeyEvent::Key::R:
        reset_simulation();
        break;
      case KeyEvent::Key::S:
        _single_step = true;
        _paused      = true;
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
        if ((event.modifiers() & KeyEvent::Modifier::Alt))
          _yarnMapperSettings.gpu_compute = !_yarnMapperSettings.gpu_compute;
        else
          _render_ground = !_render_ground;
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

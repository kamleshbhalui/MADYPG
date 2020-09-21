#include "MainApplication.h"

#include <imgui_stdlib.h>

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
}  // namespace Magnum

using namespace Magnum;
using namespace Math::Literals;

void MainApplication::reset_simulation() {
  if (!_yarnMapper)
    _yarnMapper = std::make_unique<YarnMapper>();
  else
    _yarnMapper.reset(new YarnMapper());

  _yarnMapper->m_settings = _yarnMapperSettings;

  ::Debug::log(_yarnMapper->m_settings.modelfolder);
  _yarnMapper->initialize();
  _yarnDrawable.back().setIndices(_yarnMapper->getIndices());
  _yarnDrawable.back().setVertices(_yarnMapper->getVertexData());
  _yarnDrawable.back().m_radius =
      _yarnMapper->getRadius() * _render_radius_mult;

  const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
  _meshdrawable->setIndices(mesh.F);
  _meshdrawable->setVertices(mesh.X);
}

MainApplication::MainApplication(const Arguments &arguments)
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

    _framebuffer = GL::Framebuffer{viewport};
    setupFramebuffer(vpSize);
    // GL::Renderer::setClearColor({});
    // GL::Renderer::setClearColor(_bgColor);
    // GL::defaultFramebuffer.clearColor(_bgColor);

    _yarnGeometryShader = YarnShader{};
    _meshShader         = MeshShader{};
    _ssaoShader         = SsaoShader{};
    _ssaoApplyShader    = SsaoApplyShader{};

    {  // default settings
      _yarnMapperSettings.modelfolder            = "models/model_rib";
      _yarnMapperSettings.objseq_settings.folder = "objseqs/sxsy";
      _yarnMapperSettings.objseq_settings.constant_material_space = true;
    }
    _fileDialog = std::make_unique<ImGui::FileBrowser>(
        ImGuiFileBrowserFlags_SelectDirectory |
        ImGuiFileBrowserFlags_CloseOnEsc);
    _fileDialog->SetTitle("File Dialog");
    {
      _yarnDrawable.emplace_back(_yarnGeometryShader);
      _meshdrawable = std::make_unique<MeshDrawable<MeshShader>>(_meshShader);
    }

    reset_simulation();

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
    const Vector3 eye = Vector3::zAxis(1.0f);
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
  setSwapInterval(0);  // 0 no vsync, 1 vsync
  setMinimalLoopPeriod(_min_loop_ms);
}

void MainApplication::drawEvent() {
  _profiler.beginFrame();

  bool simChanged = false;
  if (!_paused || _single_step) {  // SIM
    if (_yarnMapper->initialized()) {
      _yarnMapper->step();
      // assume no update to yarn indices
      _yarnDrawable.back().setVertices(_yarnMapper->getVertexData());

      const auto &mesh = _yarnMapper->getMeshSimulation()->getMesh();
      if (_yarnMapper->getMeshSimulation()->meshIndicesDirty())
        _meshdrawable->setIndices(mesh.F);
      // always assume changes to vertices
      _meshdrawable->setVertices(mesh.X);
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

    if (_render_yarns && _yarnMapper->initialized()) {
      _yarnGeometryShader.bindTexture(_matcap);
      _yarnGeometryShader.setProjection(_projection);
      _yarnDrawable.back().m_radius =
          _yarnMapper->getRadius() * _render_radius_mult;
      for (auto &line : _yarnDrawable) line.draw(tf);
    }

    if (_render_mesh) {
      _meshShader.setProjection(_projection);
      _meshShader.setDZ(_mesh_dz);
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

void MainApplication::viewportEvent(ViewportEvent &event) {
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

void MainApplication::setupFramebuffer(const Vector2i &size) {
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

void MainApplication::drawSettings() {
  ImGui::Begin("Render Options");
  const float spacing = 10;
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));
  ImGui::PushItemWidth(100.0f);

  if (ImGui::CollapsingHeader("SSAO", ImGuiTreeNodeFlags_DefaultOpen)) {
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
    ImGui::Unindent();
  }

  if (ImGui::CollapsingHeader("Other", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Indent();

    ImGui::Checkbox("Yarns", &_render_yarns);
    ImGui::SameLine();
    ImGui::Checkbox("Mesh", &_render_mesh);

    {
      ImGui::PushItemWidth(100.0f);
      ImGui::DragFloat("yarn radius mult", &_render_radius_mult, 0.01f, 0.0f,
                       2.0f);
      ImGui::DragFloat("mesh offset", &_mesh_dz, 0.001f, -1.0f, 1.0f);
      ImGui::PopItemWidth();
    }

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
  ImGui::Checkbox("Flat Normals", &_yarnMapperSettings.flat_normals);
  // if (_yarnMapperSettings.flat_normals) {
  //   ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
  // }
  ImGui::Checkbox("Shepard Weights", &_yarnMapperSettings.shepard_weights);
  // if (_yarnMapperSettings.flat_normals)
  //   ImGui::PopStyleVar();

  {
    std::string &txt = _yarnMapperSettings.modelfolder;
    ImGui::PushID(&txt);
    ImGui::PushItemWidth(150);
    ImGui::InputText("##txt", &txt);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Browse")) {
      _fileDialog->SetPwd(std::filesystem::path(txt).parent_path());
      _fileDialog->Open();
    }
    _fileDialog->Display();
    if (_fileDialog->HasSelected()) {
      txt = std::filesystem::relative(_fileDialog->GetSelected().string());
      _fileDialog->ClearSelected();
    }
    ImGui::PopID();
  }
  {
    std::string &txt = _yarnMapperSettings.objseq_settings.folder;
    ImGui::PushID(&txt);
    ImGui::PushItemWidth(150);
    ImGui::InputText("##txt", &txt);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Browse")) {
      _fileDialog->SetPwd(std::filesystem::path(txt).parent_path());
      _fileDialog->Open();
    }
    _fileDialog->Display();
    if (_fileDialog->HasSelected()) {
      txt = std::filesystem::relative(_fileDialog->GetSelected().string());
      _fileDialog->ClearSelected();
    }
    ImGui::PopID();
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
    ImGui::Columns(2, NULL, false);
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
}

void MainApplication::keyPressEvent(KeyEvent &event) {
  if (_imgui.handleKeyPressEvent(event))
    return;

  switch (event.key()) {
    case KeyEvent::Key::Esc:
      this->exit();
      break;
    case KeyEvent::Key::Space:
      _paused = !_paused;
      break;
    case KeyEvent::Key::R:
      reset_simulation();
      break;
    case KeyEvent::Key::S:
      _single_step = true;
      _paused      = true;
      break;
    default:
      break;
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

  const Float delta = event.offset().y();
  if (Math::abs(delta) < 1.0e-2f)
    return;

  _arcball->zoom(delta * 0.1f);

  event.setAccepted();
  redraw();
}
MAGNUM_APPLICATION_MAIN(Magnum::MainApplication)

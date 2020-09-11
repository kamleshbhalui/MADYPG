#define DECLARE_UNUSED(x) ((void)x);

#include <Magnum/Primitives/Axis.h>
#include "arcball/ArcBall.h"
#include "arcball/ArcBallCamera.h"

// #include "configure.h"
#include "Shaders/YarnShader.h"
#include "Shaders/SsaoApplyShader.h"
#include "Shaders/SsaoShader.h"
#include "render/Lines.h"
#include "render/TestYarns.h"

#include <Corrade/Containers/Optional.h>
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
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Trade/ImageData.h>
#include <MagnumPlugins/PngImporter/PngImporter.h>
#include <MagnumPlugins/JpegImporter/JpegImporter.h>
#include <MagnumPlugins/AnyImageImporter/AnyImageImporter.h>

// #include <Magnum/SceneGraph/Camera.h>
// #include <Magnum/SceneGraph/Drawable.h>
// #include <Magnum/SceneGraph/MatrixTransformation3D.h>
// #include <Magnum/SceneGraph/Object.h>
// #include <Magnum/SceneGraph/Scene.h>

#include <memory>

#include <imgui.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <random>

namespace Magnum
{
  namespace Examples
  {

    void setupTexture(GL::Texture2D &texture, Vector2i const &size,
                      GL::TextureFormat format);
    void setupTexture(GL::Texture2D &texture, Vector2i const &size,
                      GL::TextureFormat format)
    {
      texture = GL::Texture2D{};
      texture.setMagnificationFilter(GL::SamplerFilter::Linear)
          .setMinificationFilter(GL::SamplerFilter::Linear)
          .setWrapping(GL::SamplerWrapping::ClampToEdge)
          .setStorage(1, format, size);
    }

    using namespace Math::Literals;

    class SsaoExample : public Platform::Application
    {
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

      std::vector<Lines<YarnShader>> _lines;
      TestYarns _interface;
      GL::Texture2D _matcap{NoCreate};
      bool _paused = false;
      Color4 _bgColor = Color4(Color3(0.2f), 1.0f);

      // std::unique_ptr<ArcBallCamera> _arcballCamera;

      YarnShader _yarnGeometryShader{NoCreate};
      SsaoApplyShader _ssaoApplyShader{NoCreate};
      SsaoShader _ssaoShader{NoCreate};

      Magnum::Shaders::Phong _phong{NoCreate};
      // Magnum::Shaders::VertexColor3D _vcShader{NoCreate};
      // GL::Mesh _axes = MeshTools::compile(Primitives::axis3D());

      Containers::Optional<ArcBall> _arcball;
      Matrix4 _projection;
      Deg _proj_fov = 45.0_degf;
      float _proj_near = 0.01f;
      float _proj_far = 10000.0f; // TODO reduce far/near to proper scale

      GL::Mesh _mesh{NoCreate}, _screenAlignedTriangle{NoCreate};

      GL::Framebuffer _framebuffer{NoCreate};

      GL::Texture2D _albedo{NoCreate};
      GL::Texture2D _positions{NoCreate};
      GL::Texture2D _normals{NoCreate};
      GL::Texture2D _occlusion{NoCreate};
      GL::Texture2D _noise{NoCreate};

      GL::Texture2D _depth{NoCreate};

      /* Profiling */
      DebugTools::GLFrameProfiler _profiler;

      Color4 _specularColor{0.3};
      Float _radius = 1.5f;
      Float _bias = 0.5f;
      int _ao_blur_radius = 0;
      Float _ao_pow = 4.0f;

      SsaoApplyShader::Flag _ssaoApplyFlag = {};

      ImGuiIntegration::Context _imgui{NoCreate};
    };

    SsaoExample::SsaoExample(const Arguments &arguments)
        : Platform::Application{arguments, NoCreate}
    {
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
        if (!tryCreate(conf, glConf))
        {
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

        std::random_device device;
        std::default_random_engine engine(device());
        std::uniform_real_distribution<float> distr(-1, 1);
        Containers::Array<Vector4> noise(Containers::NoInit, 16);
        for (Vector4 &n : noise)
          n = Vector4{distr(engine), distr(engine), 0, 0};

        _noise = GL::Texture2D{};
        ImageView2D view{PixelFormat::RGBA32F, {4, 4}, noise};
        _noise.setMagnificationFilter(GL::SamplerFilter::Linear)
            .setMinificationFilter(GL::SamplerFilter::Linear)
            .setWrapping(GL::SamplerWrapping::Repeat)
            .setStorage(1, GL::TextureFormat::RGBA32F, {4, 4})
            .setSubImage(0, {}, view);

        const Range2Di viewport = GL::defaultFramebuffer.viewport();
        const Vector2i vpSize = viewport.size();

        _framebuffer = GL::Framebuffer{viewport};
        setupFramebuffer(vpSize);
        // GL::Renderer::setClearColor({});
        // GL::Renderer::setClearColor(_bgColor);
        // GL::defaultFramebuffer.clearColor(_bgColor);

        _yarnGeometryShader = YarnShader{};
        _ssaoShader = SsaoShader{};
        _ssaoApplyShader = SsaoApplyShader{};
        _phong = Magnum::Shaders::Phong{};

        {
          _lines.emplace_back(_yarnGeometryShader);
          _lines.back().setIndices(_interface.getIndices());
          _lines.back().setVertices(_interface.getVertexData());
        }

        PluginManager::Manager<Trade::AbstractImporter> manager;
        Trade::AnyImageImporter importer = Trade::AnyImageImporter(manager);
        if (!importer.openFile("matcaps/mcb.jpg"))
        // if (!importer.openFile("matcaps/hughsk/00036.png"))
          std::exit(2);
        Containers::Optional<Trade::ImageData2D> image = importer.image2D(0);
        CORRADE_INTERNAL_ASSERT(image);
        _matcap = GL::Texture2D();
        _matcap.setWrapping(GL::SamplerWrapping::ClampToEdge)
            .setMagnificationFilter(GL::SamplerFilter::Linear)
            .setMinificationFilter(GL::SamplerFilter::Linear, GL::SamplerMipmap::Linear)
            .setStorage(Math::log2(image->size().min()) + 1, GL::textureFormat(image->format()), image->size())
            .setSubImage(0, {}, *image)
            .generateMipmap();
      }

      /* Set up the arcball and projection */
      {
        const Vector3 eye = Vector3::zAxis(500.0f);
        const Vector3 center{};
        const Vector3 up = Vector3::yAxis();
        _arcball.emplace(eye, center, up, 45.0_degf, windowSize());
        _arcball->setLagging(0.85f);

        _projection = Matrix4::perspectiveProjection(
            _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near, _proj_far);
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

    void SsaoExample::drawEvent()
    {
      GL::defaultFramebuffer.clear(GL::FramebufferClear::Color |
                                   GL::FramebufferClear::Depth);
      _profiler.beginFrame();

      if (!_paused)
      { // SIM
        _interface.step();
        _lines.back().setVertices(_interface.getVertexData());
      }

      const bool camChanged = _arcball->updateTransformation();
      DECLARE_UNUSED(camChanged);
      const Matrix4 tf = _arcball->viewMatrix();

      /* render the scene into g-buffer */
      _framebuffer
          .mapForDraw({{YarnShader::AlbedoOutput,
                        GL::Framebuffer::ColorAttachment{0}},
                       {YarnShader::PositionsOutput,
                        GL::Framebuffer::ColorAttachment{1}},
                       {YarnShader::NormalsOutput,
                        GL::Framebuffer::ColorAttachment{2}}})
          // .clear(GL::FramebufferClear::Depth | GL::FramebufferClear::Color)
          .clearColor(0,_bgColor)
          .clearColor(1,Color4(0.0,0.0,0.0,1.0))
          .clearColor(2,Color4(0.0,0.0,1.0,1.0))
          .clearDepth(1.0)
          .bind();

      _yarnGeometryShader.bindTexture(_matcap);
      _yarnGeometryShader.setProjection(_projection);
      for (auto &line : _lines)
        line.draw(tf);

      _ssaoShader.bindNormalTexture(_normals)
          .bindNoiseTexture(_noise)
          .bindPositionTexture(_positions)
          .setProjectionMatrix(_projection)
          .setSampleRadius(_radius)
          .setBias(_bias);

      _framebuffer
          .mapForDraw({{SsaoShader::AmbientOcclusionOutput,
                        GL::Framebuffer::ColorAttachment{3}}})
          .clear(GL::FramebufferClear::Color);
      _ssaoShader.draw(_screenAlignedTriangle);

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

    void SsaoExample::viewportEvent(ViewportEvent &event)
    {
      const Vector2i fbSize = event.framebufferSize();
      const Vector2i wSize = event.windowSize();

      GL::defaultFramebuffer.setViewport({{}, fbSize});
      _framebuffer.setViewport({{}, fbSize});

      _arcball->reshape(wSize);
      _projection = Matrix4::perspectiveProjection(
          _proj_fov, Vector2{framebufferSize()}.aspectRatio(), _proj_near, _proj_far);

      setupFramebuffer(fbSize);

      _imgui.relayout(Vector2{wSize} / event.dpiScaling(), event.windowSize(),
                      fbSize);
    }

    void SsaoExample::setupFramebuffer(const Vector2i &size)
    {
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

    void SsaoExample::drawSettings()
    {
      ImGui::Begin("SSAO Options");
      const float spacing = 10;
      ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(spacing, spacing));

      ImGui::SliderFloat("SSAO radius", &_radius, 0.001f, 5.0f);
      ImGui::SliderFloat("SSAO bias", &_bias, 0.001f, 1.0f);
      ImGui::SliderInt("SSAO blur radius", &_ao_blur_radius, 0, 5);
      ImGui::SliderFloat("SSAO pow", &_ao_pow, 0.0f, 10.0f);

      if (_lines.size() > 0)
      {
        float rad = _lines[0].m_radius;
        if (ImGui::DragFloat("yarn radius", &rad, 0.1f, 0.001f, 10.0f))
        {
          for (auto &line : _lines)
          {
            line.m_radius = rad;
          }
        }
      }

      ImGui::Checkbox("Pause", &_paused);

      ImGui::SameLine();

      static bool drawOcclusion = false;
      if (ImGui::Checkbox("Show Occlusion Factor", &drawOcclusion))
      {
        if (drawOcclusion)
          _ssaoApplyFlag = SsaoApplyShader::Flag::DrawAmbientOcclusion;
        else
          _ssaoApplyFlag = {};
        _ssaoApplyShader = SsaoApplyShader{_ssaoApplyFlag};
      }

      // if (ImGui::Button("Hot Reload Shader"))
      // {
      //   Utility::Resource::overrideGroup("ssao-data", "../Shaders/resources.conf");
      //   _ssaoShader = SsaoShader{};
      //   _ssaoApplyShader = SsaoApplyShader{_ssaoApplyFlag};
      // }

      if (ImGui::Button("Reset Camera"))
        _arcball->reset();

      ImGui::SameLine();
      static bool lagging = true;
      if (ImGui::Checkbox("Camera Lagging", &lagging))
      {
        _arcball->setLagging(float(lagging) * 0.75f);
      }

      {
        static float col[3] = {_bgColor.r(), _bgColor.g(), _bgColor.b()};
        if (ImGui::ColorEdit3("BG Color", col))
        {
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

    void SsaoExample::keyPressEvent(KeyEvent &event)
    {
      if (_imgui.handleKeyPressEvent(event))
        return;

      switch (event.key())
      {
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

    void SsaoExample::keyReleaseEvent(KeyEvent &event)
    {
      if (_imgui.handleKeyReleaseEvent(event))
        return;
    }

    void SsaoExample::textInputEvent(TextInputEvent &event)
    {
      if (_imgui.handleTextInputEvent(event))
        return;
    }

    void SsaoExample::mousePressEvent(MouseEvent &event)
    {
      if (_imgui.handleMousePressEvent(event))
        return;
      /* Enable mouse capture so the mouse can drag outside of the window */
      /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
      SDL_CaptureMouse(SDL_TRUE);
      _arcball->initTransformation(event.position());
      event.setAccepted();
      redraw();
    }

    void SsaoExample::mouseReleaseEvent(MouseEvent &event)
    {
      if (_imgui.handleMouseReleaseEvent(event))
        return;
      /* Disable mouse capture again */
      /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
      SDL_CaptureMouse(SDL_FALSE);
    }

    void SsaoExample::mouseMoveEvent(MouseMoveEvent &event)
    {
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

    void SsaoExample::mouseScrollEvent(MouseScrollEvent &event)
    {
      if (_imgui.handleMouseScrollEvent(event))
      {
        /* Prevent scrolling the page */
        event.setAccepted();
        return;
      }

      const Float delta = event.offset().y();
      if (Math::abs(delta) < 1.0e-2f)
        return;

      _arcball->zoom(delta * 50); /* the Aramdillo is in mm */

      event.setAccepted();
      redraw();
    }
  } // namespace Examples
} // namespace Magnum

MAGNUM_APPLICATION_MAIN(Magnum::Examples::SsaoExample)

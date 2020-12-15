#include "PBDSimulation.h"

// PBD LIBRARY NEEDS THESE DEFINITIONS SOMEWHERE,
// OTHERWISE IT HAS A LOT OF UNDEFINEDS SYMBOLS !
#include <Utils/Timing.h>
INIT_LOGGING
INIT_TIMING

#include "../io/objio.h"
#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"
#include "PBDConstraints.h"
// #include <Common/Common.h>
#include <Simulation/Simulation.h>
#include <Simulation/TimeManager.h>
#include <Simulation/TimeStepController.h>
// #include <GenericParameters/.h>

#include <filesystem>
namespace fs = std::filesystem;

PBDSimulation::PBDSimulation(const Settings& settings) : m_settings(settings) {
  m_model = std::make_shared<PBD::SimulationModel>();
  m_model->init();
  PBD::Simulation::getCurrent()->setModel(m_model.get());
  PBD::TimeManager::getCurrent()->setTimeStepSize(m_settings.timestep);

  std::string filepath = "pbdscene/sock/cloth.obj";
  std::string obsfilepath = "pbdscene/sock/obstacle.obj";

  // std::string filepath = "pbdscene/cloth.obj";
  // std::string obsfilepath = "pbdscene/obstacle.obj";

  // cloth
  // std::string filepath = "pbdscene/cloth.obj";
  if (!loadClothMesh(filepath))
    return;

  // obstacle
  // std::string obsfilepath = "pbdscene/cube.obj";
  // if (loadBoxObstacle(obsfilepath)) {
  //   PBD::Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*m_model.get(), &m_cd);
  //   m_cd.setTolerance(static_cast<Real>(0.005));

  //   // add collider to cloth
  //   auto &tm = m_model->getTriangleModels();
  //   auto &pd = m_model->getParticles();
  //   for (unsigned int i = 0; i < tm.size(); i++)
  //   {
  //     const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
  //     unsigned int offset = tm[i]->getIndexOffset();
      // tm->setRestitutionCoeff(0.2f); // 0.6
      // tm->setFrictionCoeff(0.1f); // 0.2
  //     m_cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
  //   }
  // }

  if (loadSDFObstacle(obsfilepath)) {
    PBD::Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*m_model.get(), &m_cd2);
    m_cd2.setTolerance(static_cast<Real>(0.005));

    // add collider to cloth
    auto &tms = m_model->getTriangleModels();
    auto &pd = m_model->getParticles();
    for (unsigned int i = 0; i < tms.size(); i++)
    {
      auto& tm = tms[i];
      const unsigned int nVert = tm->getParticleMesh().numVertices();
      unsigned int offset = tm->getIndexOffset();
      tm->setFrictionCoeff(static_cast<Real>(0.1));
      m_cd2.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);


      tm->setRestitutionCoeff(0.05f); // 0.6
      tm->setFrictionCoeff(0.1f); // 0.2
      m_model->setContactStiffnessParticleRigidBody(m_settings.contact_stiffness);
    }
  }

  // set cloth stiffnesses
  m_model->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_XX, m_settings.stiffness[0]);
  m_model->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_XY, m_settings.stiffness[1]);
  m_model->setValue(PBD::SimulationModel::CLOTH_STIFFNESS_YY, m_settings.stiffness[2]);
  m_model->setValue(PBD::SimulationModel::CLOTH_BENDING_STIFFNESS, m_settings.bending_stiffness);
  // for distance constraint stiffness just use average
  m_model->setValue(PBD::SimulationModel::CLOTH_STIFFNESS, (m_settings.stiffness[0]+m_settings.stiffness[1]+m_settings.stiffness[2])/3);
  // skipping:
  // FEM: CLOTH_POISSON_RATIO_XY, CLOTH_POISSON_RATIO_YX
  // STRAIN: CLOTH_NORMALIZE_STRETCH, CLOTH_NORMALIZE_SHEAR
  // SOLID: ...

  // initial flags
  m_indicesDirty = true;
}

PBDSimulation::~PBDSimulation() {
  // can't really store simulation into a unique_ptr instead, because the
  // library might/does use getCurrent in other places...
  delete PBD::Simulation::getCurrent();
  // TimeManager::getCurrent will be deleted by simulation destructor
}


bool PBDSimulation::loadClothMesh(const std::string& filepath) {
  auto& X = m_mesh.X.cpu();
  auto& F = m_mesh.F.cpu();
  auto& U = m_mesh.U.cpu();
  auto& Fms = m_mesh.Fms.cpu();
  if (!load_obj(filepath, X, F, U, Fms))
    return false;

  // TODO try skipping this and directly using X data, similar for indices array
  AlignedVector<Vector3r> points(X.size()); 
  // Vector3r points[X.size()]; 
  for (size_t i = 0; i < X.size(); i++) {
      points[i] = X[i].map();
  }

  PBD::TriangleModel::ParticleMesh::UVs uvs;
  uvs.resize(U.size());
  for (size_t i = 0; i < U.size(); i++) {
      uvs[i] = U[i].map();
  }

  std::vector<unsigned int> indices(F.size() * 3);
  // unsigned int indices[F.size() * 3];
  int index = 0;
  for (size_t i = 0; i < F.size(); i++)
  {
    auto f = F[i].map();
    indices[index] = f[0];
    indices[index + 1] = f[1];
    indices[index + 2] = f[2];
    index += 3;
  }

  PBD::TriangleModel::ParticleMesh::UVIndices uvIndices;
  uvIndices.resize(Fms.size() * 3);
  index = 0;
  for (size_t i = 0; i < Fms.size(); i++)
  {
    auto f = Fms[i].map();
    uvIndices[index] = f[0];
    uvIndices[index + 1] = f[1];
    uvIndices[index + 2] = f[2];
    index += 3;
  }

  m_model->addTriangleModel(X.size(), F.size(), &points[0], &indices[0], uvIndices, uvs);

  // {
  // 	auto& pmesh = m_model->getTriangleModels()[0]->getParticleMesh();
  // 	auto& uvics = pmesh.getUVIndices();
  // 	auto& ics = pmesh.getFaces();
  // 	auto& puvs = pmesh.getUVs();
  // 	F.clear();
  // 	F.reserve(ics.size()/3);
  // 	for (int i = 0; i < ics.size()/3; i++) {
  // 		F.push_back({ics[i*3+0],ics[i*3+1],ics[i*3+2]});
  // 	}
  // 	Fms.clear();
  // 	Fms.reserve(uvics.size()/3);
  // 	for (int i = 0; i < uvics.size()/3; i++) {
  // 		Fms.push_back({uvics[i*3+0],uvics[i*3+1],uvics[i*3+2]});
  // 	}
  // 	// X.clear();
  // 	// X.reserve();
  // 	// for (int i = 0; i < nRows*nCols; i++) {
  // 	// 	X.push_back({points[i][0],points[i][1],points[i][2]});
  // 	// }
  // 	U.clear();
  // 	U.reserve(puvs.size());
  // 	for (size_t i = 0; i < puvs.size(); i++)
  // 		U.push_back({puvs[i][0],puvs[i][1]});
  // }

  // calculate vertex mass
  // TODO currently vertex mass doesnt do anything! maybe just remove computation then...
  std::vector<float> voronoi(X.size(),0.0f);
  std::vector<float> avguvheight(X.size(),0.0f);
  std::vector<int> avguvheight_count(X.size(),0);
  for (size_t i = 0; i < Fms.size(); i++) {
    // compute material space area
    auto ixs     = Fms[i].map();
    Vector2s e01 = U[ixs[1]].map() - U[ixs[0]].map();
    Vector2s e02 = U[ixs[2]].map() - U[ixs[0]].map();
    float A3 = std::abs(e01[0]*e02[1]-e01[1]*e02[0]) * 0.5f / 3;
    // add 1/3 of area to each vertex world-space incident to the triangle
    auto wsixs = F[i].map();
    voronoi[wsixs[0]] += A3;
    voronoi[wsixs[1]] += A3;
    voronoi[wsixs[2]] += A3;

    float v = U[ixs[0]].v+U[ixs[1]].v+U[ixs[2]].v;
    for (size_t j = 0; j < 3; ++j) {
      voronoi[wsixs[j]] += A3;
      avguvheight[wsixs[j]] += v;
      avguvheight_count[wsixs[j]] += 3;
    }
  }
  m_selected.reserve(avguvheight.size());
  for (size_t i = 0; i < avguvheight.size(); i++) {
    float avgv = avguvheight[i]/avguvheight_count[i];
    if (avgv > 0.26f)
      m_selected.push_back(uint32_t(i));
  }

  Debug::log("SELECTED:", m_selected.size());
  



  PBD::ParticleData &pd = m_model->getParticles();
  for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
  {
    pd.setMass(i, voronoi[i] * m_settings.density);
    // pd.setMass(i, 1 * m_settings.density);
    // pd.setMass(i, 1);
  }

  // // TODO SET FIXED VERTS (imgui, by height or list of vertex ids, or both, default none)
  // // Set mass of points to zero => make it static
  // for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
  // {
  //   // if(pd.getPosition(i)[1] > 0.2f)
  //   if(pd.getPosition(i)[1] > 0.09f)
  //     pd.setMass(i, 0.0f);
  // }

  // init constraints
  for (unsigned int cm = 0; cm < m_model->getTriangleModels().size(); cm++)
  {
    PBD::TriangleModel* trimodel = m_model->getTriangleModels()[cm];
    if (m_settings.simulationMethod == 1)
    {
      const unsigned int offset = trimodel->getIndexOffset();
      const unsigned int nEdges = trimodel->getParticleMesh().numEdges();
      const Utilities::IndexedFaceMesh::Edge *edges = trimodel->getParticleMesh().getEdges().data();
      for (unsigned int i = 0; i < nEdges; i++)
      {
        const unsigned int v1 = edges[i].m_vert[0] + offset;
        const unsigned int v2 = edges[i].m_vert[1] + offset;

        m_model->addDistanceConstraint(v1, v2);
      }
    }
    else if (m_settings.simulationMethod == 2)
    {
      const unsigned int offset = trimodel->getIndexOffset();
      PBD::TriangleModel::ParticleMesh &mesh = trimodel->getParticleMesh();
      const unsigned int *tris = mesh.getFaces().data();
      const unsigned int nFaces = mesh.numFaces();
      for (unsigned int i = 0; i < nFaces; i++)
      {
        const unsigned int v1 = tris[3 * i] + offset;
        const unsigned int v2 = tris[3 * i + 1] + offset;
        const unsigned int v3 = tris[3 * i + 2] + offset;
        m_model->addFEMTriangleConstraint(v1, v2, v3);
      }
    }
    else if (m_settings.simulationMethod == 3)
    {
      const unsigned int offset = trimodel->getIndexOffset();
      PBD::TriangleModel::ParticleMesh &mesh = trimodel->getParticleMesh();
      const unsigned int *tris = mesh.getFaces().data();
      const unsigned int nFaces = mesh.numFaces();
      for (unsigned int i = 0; i < nFaces; i++)
      {
        const unsigned int v1 = tris[3 * i] + offset;
        const unsigned int v2 = tris[3 * i + 1] + offset;
        const unsigned int v3 = tris[3 * i + 2] + offset;
        m_model->addStrainTriangleConstraint(v1, v2, v3);
      }
    }
    else if (m_settings.simulationMethod == 4)
    {
      const unsigned int offset = trimodel->getIndexOffset();
      PBD::TriangleModel::ParticleMesh &mesh = trimodel->getParticleMesh();
      auto& ics = mesh.getFaces();
      auto& uvics = mesh.getUVIndices();
      const unsigned int nFacesUV = uvics.size()/3;
      auto& puvs = mesh.getUVs();
      // if (nFacesUV != Fms.size() || Fms.size() != F.size())
      //   Debug::error("BROKE ASSYMPTION");
      for (unsigned int i = 0; i < nFacesUV; i++)
      {
        const Vector2r& u0 = puvs[uvics[i*3+0]];
        const Vector2r& u1 = puvs[uvics[i*3+1]];
        const Vector2r& u2 = puvs[uvics[i*3+2]];
        // const Vector2r u0 = U[Fms[i].v0].map();
        // const Vector2r u1 = U[Fms[i].v1].map();
        // const Vector2r u2 = U[Fms[i].v2].map();
        // PBD::UVFEMTriangleConstraint *c = new PBD::UVFEMTriangleConstraint();
        PBD::UVStrainTriangleConstraint *c = new PBD::UVStrainTriangleConstraint();
        const bool res = c->initConstraint(*m_model.get(), u0, u1, u2);
        c->m_bodies[0] = ics[i*3+0]; // TODO push into initConstraint maybe
        c->m_bodies[1] = ics[i*3+1];
        c->m_bodies[2] = ics[i*3+2];
        // c->m_bodies[0] = F[i].v0; // TODO push into initConstraint maybe
        // c->m_bodies[1] = F[i].v1;
        // c->m_bodies[2] = F[i].v2;
        if (res) {
          m_model->getConstraints().push_back(c);
          m_model->m_groupsInitialized = false;
        }
      }
    }

    if (m_settings.bendingMethod != 0)
    {
      const unsigned int offset = trimodel->getIndexOffset();
      PBD::TriangleModel::ParticleMesh &mesh = trimodel->getParticleMesh();
      unsigned int nEdges = mesh.numEdges();
      const PBD::TriangleModel::ParticleMesh::Edge *edges = mesh.getEdges().data();
      const unsigned int *tris = mesh.getFaces().data();
      for (unsigned int i = 0; i < nEdges; i++)
      {
        const int tri1 = edges[i].m_face[0];
        const int tri2 = edges[i].m_face[1];
        if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff))
        {
          // Find the triangle points which do not lie on the axis
          const int axisPoint1 = edges[i].m_vert[0];
          const int axisPoint2 = edges[i].m_vert[1];
          int point1 = -1;
          int point2 = -1;
          for (int j = 0; j < 3; j++)
          {
            if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2))
            {
              point1 = tris[3 * tri1 + j];
              break;
            }
          }
          for (int j = 0; j < 3; j++)
          {
            if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2))
            {
              point2 = tris[3 * tri2 + j];
              break;
            }
          }
          if ((point1 != -1) && (point2 != -1))
          {
            const unsigned int vertex1 = point1 + offset;
            const unsigned int vertex2 = point2 + offset;
            const unsigned int vertex3 = edges[i].m_vert[0] + offset;
            const unsigned int vertex4 = edges[i].m_vert[1] + offset;
            if (m_settings.bendingMethod == 1)
              m_model->addDihedralConstraint(vertex1, vertex2, vertex3, vertex4);
            else if (m_settings.bendingMethod == 2)
              m_model->addIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
          }
        }
      }
    }
  }

  return true;
}

bool PBDSimulation::loadBoxObstacle(const std::string& filepath) {
  m_obstacles.resize(m_obstacles.size()+1);
  auto& oX = m_obstacles.back().mesh.X.cpu();
  auto& oF = m_obstacles.back().mesh.F.cpu();
  if (!load_obj(filepath, oX, oF)) {
    m_obstacles.resize(m_obstacles.size()-1);
    return false;
  }

  { // obstacle mesh to pbd rigidbody mesh
    Utilities::IndexedFaceMesh mesh;
    PBD::VertexData vd;

    const unsigned int nPoints = static_cast<unsigned int>(oX.size());
    const unsigned int nFaces = static_cast<unsigned int>(oF.size());
    mesh.initMesh(nPoints, nFaces * 2, nFaces);
    vd.reserve(nPoints);
    for (unsigned int i = 0; i < nPoints; i++)
      vd.addVertex(oX[i].map());
    for (unsigned int i = 0; i < nFaces; i++)
    {
      int posIndices[3];
      auto f = oF[i].map();
      for (int j = 0; j < 3; j++)
        posIndices[j] = f[j];
      mesh.addFace(&posIndices[0]);
    }
    mesh.buildNeighbors();
    mesh.updateNormals(vd, 0);
    mesh.updateVertexNormals(vd);


    PBD::SimulationModel::RigidBodyVector &rb = m_model->getRigidBodies();
    rb.resize(1);
    rb[0] = new PBD::RigidBody();
    rb[0]->initBody(1.0,
      // Vector3r(0.0, 0.0, 0.0),
      Vector3r(0.0, -0.095792, -0.0477),
      Quaternionr(1.0, 0.0, 0.0, 0.0),
      vd, mesh,
      // Vector3r(1.0, 1.0, 1.0)); // scale
      Vector3r(0.07,0.07,0.07)); // scale, doesnt seem to matter for box CD
    rb[0]->setMass(0.0); // static

    // my rendering trafo
    m_obstacles[0].transformation.scale = 0.07;
    m_obstacles[0].transformation.translation = Vector3r(0.0, -0.095792, -0.0477);

      
    const std::vector<Vector3r> *vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
    const unsigned int nVert1 = static_cast<unsigned int>(vertices1->size());
    m_cd.addCollisionBox(0, // rigid body index
      PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices1)[0], nVert1, 
      Vector3r(0.07,0.07,0.07)); 
  }

  return true;
}

bool PBDSimulation::loadSDFObstacle(const std::string& filepath) {
  m_obstacles.resize(m_obstacles.size()+1);
  auto& oX = m_obstacles.back().mesh.X.cpu();
  auto& oF = m_obstacles.back().mesh.F.cpu();
  if (!load_obj(filepath, oX, oF)) {
    m_obstacles.resize(m_obstacles.size()-1);
    return false;
  }

  // NOTE: there might be some redundancy/useless copies of mesh data between rigidbody and sdf etc. that could be avoided

  Utilities::IndexedFaceMesh mesh;
  PBD::VertexData vd;

  const unsigned int nPoints = static_cast<unsigned int>(oX.size());
  const unsigned int nFaces = static_cast<unsigned int>(oF.size());
  mesh.initMesh(nPoints, nFaces * 2, nFaces);
  vd.reserve(nPoints);
  for (unsigned int i = 0; i < nPoints; i++)
    vd.addVertex(oX[i].map());
  for (unsigned int i = 0; i < nFaces; i++)
  {
    int posIndices[3];
    auto f = oF[i].map();
    for (int j = 0; j < 3; j++)
      posIndices[j] = f[j];
    mesh.addFace(&posIndices[0]);
  }
  mesh.buildNeighbors();
  mesh.updateNormals(vd, 0);
  mesh.updateVertexNormals(vd);


  PBD::SimulationModel::RigidBodyVector &rbs = m_model->getRigidBodies();
  rbs.push_back(new PBD::RigidBody());
  auto rb = rbs.back();
  rb->initBody(1.0,
    Vector3r(0.0, 0.0, 0.0),
    Quaternionr(1.0, 0.0, 0.0, 0.0),
    vd, mesh,
    Vector3r(1.0, 1.0, 1.0)); // scale
  rb->setMass(0.0); // static

  std::vector<unsigned int> &faces = mesh.getFaces();
      
  // Discregrid::TriangleMesh sdfMesh(&vd.getPosition(0)[0], faces.data(), vd.size(), nFaces);
  // if type is float, copy vector to double vector
  std::vector<double> doubleVec;
  doubleVec.resize(3 * vd.size());
  for (unsigned int i = 0; i < vd.size(); i++)
    for (unsigned int j = 0; j < 3; j++)
      doubleVec[3 * i + j] = vd.getPosition(i)[j];
  Discregrid::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), vd.size(), nFaces);

  Discregrid::MeshDistance md(sdfMesh);
  Eigen::AlignedBox3d domain;
  for (auto const& x : sdfMesh.vertices()) {
    domain.extend(x);
  }
  Eigen::Vector3d extents = domain.max() - domain.min();
  auto resolution = std::array<unsigned int, 3>({
    static_cast<unsigned int>(float(extents[0]) * m_settings.sdf_samples_per_m)+1,
    static_cast<unsigned int>(float(extents[1]) * m_settings.sdf_samples_per_m)+1,
    static_cast<unsigned int>(float(extents[2]) * m_settings.sdf_samples_per_m)+1
  });
  domain.max() += 0.01 * Eigen::Vector3d::Ones();
  domain.min() -= 0.01 * Eigen::Vector3d::Ones();
  // Debug::log("res",resolution[0],resolution[1],resolution[2]);

  m_rbsdfs.push_back(std::make_shared<PBD::CubicSDFCollisionDetection::Grid>(domain, resolution));

  auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
  func = [&md](Eigen::Vector3d const& xi) {return md.signedDistanceCached(xi); };
  m_rbsdfs.back()->addFunction(func, true);
  // std::string cachePath = filepath + "_cache";
  // if (FileSystem::makeDir(cachePath) == 0)
  // {
  //   LOG_INFO << "Save SDF: " << sdfFileName;
  //   distanceFields[sdfFileName]->save(sdfFileName);
  //   FileSystem::writeMD5File(rbd.m_modelFile, md5FileName);
  // }

  const std::vector<Vector3r> *vertices = rb->getGeometry().getVertexDataLocal().getVertices();
  const unsigned int nVert = static_cast<unsigned int>(vertices->size());
    
  Vector3r scale = Vector3r::Ones();
  bool testMesh = true;
  bool invertSDF = false;
  m_cd2.addCubicSDFCollisionObject(m_rbsdfs.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, m_rbsdfs.back(), scale, testMesh, invertSDF);

  return true;
}

void PBDSimulation::update() {
  m_indicesDirty = false;

  auto ts = PBD::Simulation::getCurrent()->getTimeStep();
  
	PBD::Simulation::getCurrent()->getTimeStep()->setValue(PBD::TimeStepController::MAX_ITERATIONS, static_cast<unsigned int>(m_settings.iterations));
  for (unsigned int i = 0; i < m_settings.substeps; i++)
    ts->step(*m_model.get());
  

  // sync mesh update
  auto& particles = m_model->getParticles();
  auto& X = m_mesh.X.cpu();
  for (int i = 0; i < particles.size(); i++) {
    auto& p = particles.getPosition(i);
    X[i].map() << p[0],p[1],p[2];
  }
}



void PBDSimulation::applyForce(float fx , float fy, float fz) {
  Vector3s force;
  force << fx, fy, fz;
  auto *tm = PBD::TimeManager::getCurrent();
  const Real h = tm->getTimeStepSize();
  auto &pd = m_model->getParticles();
  // float factor = 1.0f/pd.size();
  // for (unsigned int j = 0; j < pd.size(); j++)
  // {
  //   const Real mass = pd.getMass(j);
  //   if (mass != 0.0f)
  //   {  pd.getVelocity(j) += factor*pd.getInvMass(j)*h*force;
  // // Debug::log("applying force: ",(pd.getInvMass(j)*force*h).transpose());
  //   }
  // }

  float factor = 1.0f/m_selected.size();
  for (auto i : m_selected) {
    pd.getVelocity(i) += factor * pd.getInvMass(i) * h * force; // note setting, not adding
    // pd.getPosition(i) += Vector3s(0.0,0.01,0.0); 
  }
}
#include "PBDSimulation.h"

// PBD LIBRARY NEEDS THESE DEFINITIONS SOMEWHERE,
// OTHERWISE IT HAS LOTS OF UNDEFINEDS SYMBOLS !
#include <Utils/Timing.h>
INIT_LOGGING
INIT_TIMING

#include "../utils/debug_logging.h"
#include "../utils/threadutils.h"
// #include <Common/Common.h>
#include <Simulation/Simulation.h>
#include <Simulation/TimeManager.h>
#include <Simulation/TimeStepController.h>

PBDSimulation::PBDSimulation(const Settings& settings) : m_settings(settings) {
	// TODO CONTINUE HERE: in any order
	// load cloth mesh
	// load obstacle mesh / add collision constraints
	// add pull method
	m_model = std::make_shared<PBD::SimulationModel>();
	m_model->init();
	PBD::Simulation::getCurrent()->setModel(m_model.get());
	PBD::TimeManager::getCurrent()->setTimeStepSize(m_settings.timestep);
	loadMesh();

	// initial flags
	m_indicesDirty = true;
}

void PBDSimulation::loadMesh() {
	const int nRows = 30;
	const int nCols = 30;
	const float width = 0.2;
	const float height = 0.2;
	
	PBD::TriangleModel::ParticleMesh::UVs uvs;
	uvs.resize(nRows*nCols);

	const float dy = width / float(nCols - 1);
	const float dx = height / float(nRows - 1);

	Vector3r points[nRows*nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			const float y = float(dy*j) - height/2;
			const float x = float(dx*i) - width/2;
			points[i*nCols + j] = Vector3r(x, 0.0, y);

			uvs[i*nCols + j][0] = x;
			uvs[i*nCols + j][1] = y;
		}
	}
	const int nIndices = 6 * (nRows - 1)*(nCols - 1);

	PBD::TriangleModel::ParticleMesh::UVIndices uvIndices;
	uvIndices.resize(nIndices);

	unsigned int indices[nIndices];
	int index = 0;
	for (int i = 0; i < nRows - 1; i++)
	{
		for (int j = 0; j < nCols - 1; j++)
		{
			int helper = 0;
			if (i % 2 == j % 2)
				helper = 1;

			indices[index] = i*nCols + j;
			indices[index + 1] = i*nCols + j + 1;
			indices[index + 2] = (i + 1)*nCols + j + helper;

			uvIndices[index] = i*nCols + j;
			uvIndices[index + 1] = i*nCols + j + 1;
			uvIndices[index + 2] = (i + 1)*nCols + j + helper;
			index += 3;

			indices[index] = (i + 1)*nCols + j + 1;
			indices[index + 1] = (i + 1)*nCols + j;
			indices[index + 2] = i*nCols + j + 1 - helper;

			uvIndices[index] = (i + 1)*nCols + j + 1;
			uvIndices[index + 1] = (i + 1)*nCols + j;
			uvIndices[index + 2] = i*nCols + j + 1 - helper;
			index += 3;
		}
	}

	m_model->addTriangleModel(nRows*nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);
	
	// TODO DEFINE MASS ? DENSITY ?
	PBD::ParticleData &pd = m_model->getParticles();
	float trisize = 0.5f * width / float(nCols - 1) * height / float(nRows - 1);
	for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
	{
		pd.setMass(i, trisize * m_settings.density);
	}
	Debug::log("Setting mass:", trisize * m_settings.density);

	// TODO SET FIXED VERTS
	// Set mass of points to zero => make it static
	pd.setMass(0, 0.0);
	pd.setMass((nRows-1)*nCols, 0.0);

	// init constraints
	for (unsigned int cm = 0; cm < m_model->getTriangleModels().size(); cm++)
	{
		if (m_settings.simulationMethod == 1)
		{
			const unsigned int offset = m_model->getTriangleModels()[cm]->getIndexOffset();
			const unsigned int nEdges = m_model->getTriangleModels()[cm]->getParticleMesh().numEdges();
			const Utilities::IndexedFaceMesh::Edge *edges = m_model->getTriangleModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;

				m_model->addDistanceConstraint(v1, v2);
			}
		}
		else if (m_settings.simulationMethod == 2)
		{
			const unsigned int offset = m_model->getTriangleModels()[cm]->getIndexOffset();
			PBD::TriangleModel::ParticleMesh &mesh = m_model->getTriangleModels()[cm]->getParticleMesh();
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
			const unsigned int offset = m_model->getTriangleModels()[cm]->getIndexOffset();
			PBD::TriangleModel::ParticleMesh &mesh = m_model->getTriangleModels()[cm]->getParticleMesh();
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

		if (m_settings.bendingMethod != 0)
		{
			const unsigned int offset = m_model->getTriangleModels()[cm]->getIndexOffset();
			PBD::TriangleModel::ParticleMesh &mesh = m_model->getTriangleModels()[cm]->getParticleMesh();
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
	Debug::log("Num mesh triangles:", nIndices/3);
	Debug::log("Num mesh verts:", nRows*nCols);

	// sync initial mesh
	auto& F = m_mesh.F.cpu();
	auto& Fms = m_mesh.Fms.cpu();
	F.reserve(nIndices/3);
	Fms.reserve(nIndices/3);
	for (int i = 0; i < nIndices/3; i++) {
		F.push_back({indices[i*3+0],indices[i*3+1],indices[i*3+2]});
		Fms.push_back({uvIndices[i*3+0],uvIndices[i*3+1],uvIndices[i*3+2]});
	}
	auto& X = m_mesh.X.cpu();
	auto& U = m_mesh.U.cpu();
	X.reserve(nRows*nCols);
	for (int i = 0; i < nRows*nCols; i++) {
		X.push_back({points[i][0],points[i][1],points[i][2]});
	}
	U.reserve(X.size());
	for (size_t i = 0; i < X.size(); i++)
		U.push_back({X[i].x,X[i].z});

}

void PBDSimulation::update() {
	m_indicesDirty = false;

	for (unsigned int i = 0; i < m_settings.substeps; i++)
    PBD::Simulation::getCurrent()->getTimeStep()->step(*m_model.get());

	// sync mesh update
	auto& particles = m_model->getParticles();
	auto& X = m_mesh.X.cpu();
	for (int i = 0; i < particles.size(); i++) {
		auto& p = particles.getPosition(i);
		X[i].map() << p[0],p[1],p[2];
	}
}
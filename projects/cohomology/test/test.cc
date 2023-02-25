/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#include <cohomology.h>
#include <verification.h>
#include "examples.h"

#include <gtest/gtest.h>

#include <unordered_set>
#include <vector>
#include <algorithm>

/**
 * @brief Contains the unit tests.
 */
namespace projects::cohomology::test
{

TEST(projects_cohomology, FaceVertices)
{
	FaceVertices a, b, c;
	a.vertexIndices = { 1, 2, 3, 4 };
	b.vertexIndices = { 2, 3, 4, 1 };
	c.vertexIndices = { 4, 3, 2, 1 };

	for(FaceVertices const& fv : { a, b, c })
		for(quiver::vertex_index_t v = 0; v < 10; ++v)
			EXPECT_EQ(fv.contains(v), 1 <= v && v <= 4);

	for(quiver::vertex_index_t from = 1; from < 4; ++from)
	{
		EXPECT_EQ(a.orderOf(from, from + 1),  1);
		EXPECT_EQ(b.orderOf(from, from + 1),  1);
		EXPECT_EQ(c.orderOf(from, from + 1), -1);
	}

	EXPECT_EQ(a.relativeOrientation(a), RelativeFaceOrientation::inAgreement);
	EXPECT_EQ(a.relativeOrientation(b), RelativeFaceOrientation::inAgreement);
	EXPECT_EQ(a.relativeOrientation(c), RelativeFaceOrientation::inDisagreement);

	EXPECT_EQ(b.relativeOrientation(a), RelativeFaceOrientation::inAgreement);
	EXPECT_EQ(b.relativeOrientation(b), RelativeFaceOrientation::inAgreement);
	EXPECT_EQ(b.relativeOrientation(c), RelativeFaceOrientation::inDisagreement);

	EXPECT_EQ(c.relativeOrientation(a), RelativeFaceOrientation::inDisagreement);
	EXPECT_EQ(c.relativeOrientation(b), RelativeFaceOrientation::inDisagreement);
	EXPECT_EQ(c.relativeOrientation(c), RelativeFaceOrientation::inAgreement);
}

TEST(projects_cohomology, error_non_orientable_mesh)
{
	const auto lfMeshCube = examples::generateCube();
	Mesh meshCube;
	EXPECT_NO_THROW(meshCube = Mesh(*lfMeshCube));

	const auto lfMeshMoebiusStrip = examples::generateMoebiusStrip();
	Mesh meshMoebiusStrip;
	EXPECT_THROW(meshMoebiusStrip = Mesh(*lfMeshMoebiusStrip), NonOrientableMeshError);
}

TEST(projects_cohomology, error_connected_ports)
{
	const auto lfMeshCube = examples::generateCube();
	Mesh meshCube;
	for(PortType portType1 : { electricPort, magneticPort })
	{
		for(PortType portType2 : { electricPort, magneticPort })
		{
			const unsigned portFace1 = 1;
			unsigned portFace2;
			const auto portMap = [&](unsigned f) -> PortInformation
			{
				if(f == portFace1)
					return { portType1, 0 };
				else if(f == portFace2)
					return { portType2, 1 };
				else
					return {};
			};
			for(portFace2 = (portFace1 + 1); portFace2 < (portFace1 + 1) + 4; ++portFace2)
				EXPECT_THROW(meshCube = Mesh(*lfMeshCube, portMap), ConnectedPortsError);
			portFace2 = 0;
			EXPECT_NO_THROW(meshCube = Mesh(*lfMeshCube, portMap));
		}
	}
}

TEST(projects_cohomology, error_invalid_edge)
{
	// Test for the exception InvalidEdgeError but this is really trivial and not really needed.
}

TEST(projects_cohomology, error_no_non_excluded_arc)
{
	const auto lfMeshCube = examples::generateCube();
	const auto portMap = [&](unsigned f) -> PortInformation
	{
		if(f == 0 || f == 1)
			return { electricPort, PortInformation::PortIdType(f) };
		else
			return {};
	};
	const Mesh meshCube = Mesh(*lfMeshCube, portMap);
	for(unsigned v = 2; v < 6; ++v)
	{
		Walk walk(&meshCube, meshCube.faceGraph.V[v].vertexIndices);
		EXPECT_EQ(walk.isLoop(), false);
		EXPECT_EQ(walk.vertices.size(), 4);

		walk.vertices.push_back(walk.vertices.front());
		EXPECT_EQ(walk.isLoop(), true);
		EXPECT_EQ(walk.vertices.size(), 5);

		EXPECT_THROW(pushLoopsOutOfPorts(meshCube, walk), NoNonExcludedArcError);
	}
}

TEST(projects_cohomology, excluded_electric)
{
	const PortInformation portInformation(electricPort, 0);
	const auto lfMesh = examples::generateCube(examples::Shape::quads);
	const auto portMap = [=](unsigned f) -> PortInformation
	{
		if(f == 0)
			return portInformation;
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);

	const std::unordered_set<quiver::vertex_index_t> excludedFaces = { 0 };
	const std::unordered_set<quiver::vertex_index_t> excludedEdges = { 1, 2, 6, 9 };
	const std::unordered_set<quiver::vertex_index_t> excludedVertices = { 0, 2, 4, 6 };

	for(quiver::vertex_index_t f = 0; f < mesh.faces(); ++f)
	{
		EXPECT_EQ(mesh.isExcludedFace(f), excludedFaces.contains(f));
	}

	for(quiver::vertex_index_t e = 0; e < mesh.edges(); ++e)
	{
		EXPECT_EQ(mesh.isExcludedEdge(e), excludedEdges.contains(e));
		EXPECT_EQ(mesh.isPortBoundaryEdge(e), excludedEdges.contains(e) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isElectricalBoundaryEdge(e), excludedEdges.contains(e));
	}

	for(quiver::vertex_index_t v = 0; v < mesh.vertices(); ++v)
	{
		EXPECT_EQ(mesh.isExcludedVertex(v), excludedVertices.contains(v));
		EXPECT_EQ(mesh.isPortVertex(v), excludedVertices.contains(v) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isPortBoundaryVertex(v), excludedVertices.contains(v) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isPortInteriorVertex(v), std::nullopt);
	}
}

TEST(projects_cohomology, excluded_magnetic)
{
	const PortInformation portInformation(magneticPort, 0);
	const auto lfMesh = examples::generateCube(examples::Shape::quads);
	const auto portMap = [=](unsigned f) -> PortInformation
	{
		if(f == 0 || f == 2 || f == 4)
			return portInformation;
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);

	const std::unordered_set<quiver::vertex_index_t> excludedFaces = { 0, 2, 4 };
	const std::unordered_set<quiver::vertex_index_t> excludedEdges = { 0, 1, 2 };
	const std::unordered_set<quiver::vertex_index_t> boundaryEdges = { 3, 4, 5, 6, 8, 9 };
	const std::unordered_set<quiver::vertex_index_t> excludedVertices = { 0, 1, 2, 4 };
	const std::unordered_set<quiver::vertex_index_t> portVertices = { 0, 1, 2, 3, 4, 5, 6 };

	for(quiver::vertex_index_t f = 0; f < mesh.faces(); ++f)
	{
		EXPECT_EQ(mesh.isExcludedFace(f), excludedFaces.contains(f));
	}

	for(quiver::vertex_index_t e = 0; e < mesh.edges(); ++e)
	{
		EXPECT_EQ(mesh.isExcludedEdge(e), excludedEdges.contains(e));
		EXPECT_EQ(mesh.isPortBoundaryEdge(e), boundaryEdges.contains(e) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isElectricalBoundaryEdge(e), false);
	}

	for(quiver::vertex_index_t v = 0; v < mesh.vertices(); ++v)
	{
		EXPECT_EQ(mesh.isExcludedVertex(v), excludedVertices.contains(v));
		EXPECT_EQ(mesh.isPortVertex(v), portVertices.contains(v) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isPortBoundaryVertex(v), portVertices.contains(v) && !excludedVertices.contains(v) ? std::optional(portInformation) : std::nullopt);
		EXPECT_EQ(mesh.isPortInteriorVertex(v), excludedVertices.contains(v) ? std::optional(portInformation) : std::nullopt);
	}
}

TEST(projects_cohomology, walk)
{
	const auto lfMesh = examples::generateCube();
	const Mesh mesh = Mesh(*lfMesh);

	using VertexPair = std::array<quiver::vertex_index_t, 2>;
	constexpr auto sameElements = [](VertexPair const& lhs, VertexPair const& rhs) constexpr noexcept
	{
		return lhs[0] == rhs[0] && lhs[1] == rhs[1]
			|| lhs[0] == rhs[1] && lhs[1] == rhs[0];
	};

	const Walk walk(&mesh, { 0, 2, 6, 4 });
	EXPECT_EQ(std::distance(walk.begin(), walk.end()), 3);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(walk.begin()[0]).edgeIndex].vertices, { 0, 2 }), true);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(walk.begin()[1]).edgeIndex].vertices, { 2, 6 }), true);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(walk.begin()[2]).edgeIndex].vertices, { 6, 4 }), true);
	EXPECT_EQ(walk.isEmpty(), false);
	EXPECT_EQ(walk.isClosed(), false);
	EXPECT_EQ(walk.realVertices(), 4);
	EXPECT_EQ(walk.edges(), 3);

	Walk walkRotated = walk;
	for(unsigned r = 1; r <= walk.realVertices(); ++r)
	{
		walkRotated.rotate(1);
		Walk walkRotatedDirectly = walk;
		walkRotatedDirectly.rotate(r);
		EXPECT_EQ(walkRotated, walkRotatedDirectly);
	}

	const Walk loop(&mesh, { 0, 2, 6, 4, 0 });
	EXPECT_EQ(std::distance(loop.begin(), loop.end()), 4);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(loop.begin()[0]).edgeIndex].vertices, { 0, 2 }), true);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(loop.begin()[1]).edgeIndex].vertices, { 2, 6 }), true);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(loop.begin()[2]).edgeIndex].vertices, { 6, 4 }), true);
	EXPECT_EQ(sameElements(mesh.edgeConnectivity[std::get<2>(loop.begin()[3]).edgeIndex].vertices, { 4, 0 }), true);
	EXPECT_EQ(loop.isEmpty(), false);
	EXPECT_EQ(loop.isClosed(), true);
	EXPECT_EQ(loop.realVertices(), 4);
	EXPECT_EQ(loop.edges(), 4);

	Walk loopRotated = loop;
	for(unsigned r = 1; r <= loop.realVertices(); ++r)
	{
		loopRotated.rotate(1);
		Walk loopRotatedDirectly = loop;
		loopRotatedDirectly.rotate(r);
		EXPECT_EQ(loopRotated, loopRotatedDirectly);
	}
}

TEST(projects_cohomology, handlebodies)
{
	for(unsigned genus = 0; genus <= 20; ++genus)
	{
		for(examples::Shape shape : { examples::Shape::tris, examples::Shape::quads })
		{
			const auto lfMesh = examples::generate3Handlebody(genus, shape);
			const Mesh mesh(*lfMesh);

			const RelativeCohomology relativeCohomology(mesh);

			EXPECT_EQ(relativeCohomology.topological().size(), genus * 2);
			EXPECT_EQ(relativeCohomology.electric().size(), 0);
			EXPECT_EQ(relativeCohomology.magnetic().size(), 0);
		}
	}
}

TEST(projects_cohomology, verification)
{
	using Scenario = std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>>;

	std::vector<Scenario> scenarios;
	scenarios.emplace_back(examples::getCubeExampleWithoutPorts());
	scenarios.emplace_back(examples::getCubeExampleWithTwoPorts());
	scenarios.emplace_back(examples::getSmallTorusExample1());
	scenarios.emplace_back(examples::getSmallTorusExample2());
	scenarios.emplace_back(examples::getBigTorusExample());
	scenarios.emplace_back(examples::getBigElectricTorusExample());
	scenarios.emplace_back(examples::getNearFailTorusExample());
	scenarios.emplace_back(examples::getCounterexampleTorusExample());
	for(unsigned genus = 0; genus <= 20; ++genus)
		scenarios.emplace_back(examples::getHandlebodyExample(genus));

	for(auto const& [mesh, lfMesh] : scenarios)
	{
		const RelativeCohomology relativeCohomology(mesh);

		const auto portCount = mesh.countPorts();
		EXPECT_EQ(portCount.electricPorts == 0 ? 0 : portCount.electricPorts - 1, relativeCohomology.electric().size());
		EXPECT_EQ(portCount.magneticPorts == 0 ? 0 : portCount.magneticPorts - 1, relativeCohomology.magnetic().size());

		EXPECT_EQ(verification::verifyGenerators(mesh, relativeCohomology), true);
	}
}

}  // namespace projects::cohomology::test

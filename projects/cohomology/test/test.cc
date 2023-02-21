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

#include <vector>
#include <algorithm>

/**
 * @brief Contains the unit tests.
 */
namespace projects::cohomology::test
{

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
	for(PortType portType1 : { 	electricPort, magneticPort })
	{
		for(PortType portType2 : { 	electricPort, magneticPort })
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
	for(unsigned i = 2; i < 6; ++i)
	{
		Walk walk(&meshCube, meshCube.faceGraph.V[i].vertexIndices);
		EXPECT_EQ(walk.isLoop(), false);
		EXPECT_EQ(walk.vertices.size(), 4);

		walk.vertices.push_back(walk.vertices.front());
		EXPECT_EQ(walk.isLoop(), true);
		EXPECT_EQ(walk.vertices.size(), 5);

		EXPECT_THROW(pushLoopsOutOfPorts(meshCube, walk), NoNonExcludedArcError);
	}
}

TEST(projects_cohomology, handlebodies)
{
	for(unsigned genus = 0; genus <= 10; ++genus)
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
	scenarios.emplace_back(examples::getCubeExample());
	scenarios.emplace_back(examples::getSmallTorusExample1());
	scenarios.emplace_back(examples::getSmallTorusExample2());
	scenarios.emplace_back(examples::getBigTorusExample());
	scenarios.emplace_back(examples::getBigElectricTorusExample());
	scenarios.emplace_back(examples::getNearFailTorusExample());
	scenarios.emplace_back(examples::getCounterexampleTorusExample());

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

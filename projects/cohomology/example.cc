/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#include <cohomology.h>
#include <verification.h>
#include <visualization.h>
#include <test/examples.h>

#include <iostream>
#include <fstream>
#include <vector>

int main()
{
	using namespace projects::cohomology;

	const auto [mesh, lfMesh] = examples::getBigTorusExample();

	auto topologicalCycles = getTopologicalCycles(mesh);
	auto electricConnectors = getElectricConnectors(mesh);
	auto magneticPortCycles = getMagneticPortCycles(mesh);

	{
		std::ofstream meshObj("mesh.obj", std::ofstream::trunc);
		visualization::exportMeshObj(meshObj, *lfMesh, mesh);

		std::ofstream topologicalCyclesObj("topologicalCycles.obj", std::ofstream::trunc);
		visualization::exportWalksObj(topologicalCyclesObj, *lfMesh, topologicalCycles);

		std::ofstream electricConnectorsObj("electricConnectors.obj", std::ofstream::trunc);
		visualization::exportWalksObj(electricConnectorsObj, *lfMesh, electricConnectors);

		std::ofstream magneticPortCyclesObj("magneticPortCycles.obj", std::ofstream::trunc);
		visualization::exportEdgeVectorsObj(magneticPortCyclesObj, *lfMesh, magneticPortCycles);

		std::vector<Walk> excludedEdges;
		for(EdgeIndex::Index i = 0; i < mesh.edgeConnectivity.size(); ++i)
		{
			if(mesh.isExcludedEdge(i))
			{
				auto const& vertices = mesh.edgeConnectivity[i].vertices;
				Walk walk(&mesh);
				walk.vertices.assign(vertices.begin(), vertices.end());
				excludedEdges.emplace_back(std::move(walk));
			}
		}
		std::ofstream excludedEdgesObj("excludedEdges.obj", std::ofstream::trunc);
		visualization::exportWalksObj(excludedEdgesObj, *lfMesh, excludedEdges);
	}

	std::cout << "vertices: " << mesh.vertices() << "\tedges: " << mesh.edges() << "\tfaces: " << mesh.faces() << '\n';

	const auto ports = mesh.countPorts();
	std::cout << "distinct electric ports: " << ports.electricPorts << "\tdistinct magnetic ports: " << ports.magneticPorts << '\n';

	std::cout << "topological cycles: " << topologicalCycles.size() << '\n';
	std::cout << "electric connectors: " << electricConnectors.size() << '\n';
	std::cout << "magnetic port cycles: " << magneticPortCycles.size() << '\n';

	RelativeCohomology relativeCohomology(mesh, topologicalCycles, electricConnectors, magneticPortCycles);
	std::cout << "verification: " << std::boolalpha << verification::verifyGenerators(mesh, relativeCohomology) << '\n';
}

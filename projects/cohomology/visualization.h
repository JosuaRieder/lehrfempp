/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#ifndef COHOMOLOGY_VISUALIZATION_H_INCLUDED
#define COHOMOLOGY_VISUALIZATION_H_INCLUDED

#include <cohomology.h>
#include <utility.h>

#include <Eigen/Dense>
#include <quiver.hpp>

#include <ostream>
#include <type_traits>
#include <stdexcept>
#include <string>

/**
 * @brief Contains the verification functionality of the cohomology subproject.
 */
namespace projects::cohomology::visualization
{

/**
 * @brief This is an internal function, do not use directly.
 */
void exportVerticesObj(std::ostream& stream, lf::mesh::Mesh const& lfMesh)
{
	for(std::size_t v = 0; v < lfMesh.NumEntities(2); ++v)
	{
		Eigen::Vector3d const& pos = lfMesh.EntityByIndex(2, v)->Geometry()->Global(Eigen::MatrixXd::Zero(0, 1));
		stream << 'v';
		for(unsigned i = 0; i < 3; ++i)
			stream << ' ' << pos[i];
		stream << '\n';
	}
}

/**
 * @brief This is an internal function, do not use directly.
 */
void printPortTypeMaterial(std::ostream& stream, PortType portType)
{
	stream << "usemtl ";
	switch(portType)
	{
	case electricPort:
		stream << "electricPort";
		break;
	case magneticPort:
		stream << "magneticPort";
		break;
	default:
		stream << "noPort";
	}
	stream << '\n';
}

/**
 * @brief This is an internal function, do not use directly.
 */
void exportMaterials(std::ostream& stream)
{
	// ParaView assigns IDs to materials based on the first occurence in the obj file.
	// For the IDs to be (semantically) consistent across different files, we should always list all materials first.
	printPortTypeMaterial(stream, noPort);
	printPortTypeMaterial(stream, electricPort);
	printPortTypeMaterial(stream, magneticPort);
}

/**
 * @brief Exports a mesh to a stream in the OBJ file format.
 * As this overload only takes an `lf::mesh::Mesh`, no ports will be exported.
 * @param stream The stream to write to.
 * @param lfMesh The mesh to export.
 */
void exportMeshObj(std::ostream& stream, lf::mesh::Mesh const& lfMesh)
{
	exportVerticesObj(stream, lfMesh);
	exportMaterials(stream);

	printPortTypeMaterial(stream, noPort);
	for(lf::mesh::Entity const* face : lfMesh.Entities(0))
	{
		stream << 'f';
		for(lf::mesh::Entity const* vertex : face->SubEntities(2))
			stream << ' ' << lfMesh.Index(*vertex) + 1;
		stream << '\n';
	}
}

/**
 * @brief Exports a mesh to a stream in the OBJ file format.
 * This overload takes both a `lf::mesh::Mesh` as well as a `Mesh` and thus also exports the ports found in `mesh`.
 * @param stream The stream to write to.
 * @param lfMesh The mesh whose vertices are exported.
 * @param mesh The mesh whose faces are exported.
 */
void exportMeshObj(std::ostream& stream, lf::mesh::Mesh const& lfMesh, Mesh const& mesh)
{
	exportVerticesObj(stream, lfMesh);
	exportMaterials(stream);

	PortType previousMat = static_cast<PortType>(-1);
	for(auto const& face : mesh.faceGraph.V)
	{
		if(previousMat != face.portType)
		{
			printPortTypeMaterial(stream, face.portType);
			previousMat = face.portType;
		}

		stream << 'f';
		for(const auto v : face.vertexIndices)
			stream << ' ' << v + 1;
		stream << '\n';
	}
}

/**
 * @brief Exports a range of `Walk`s to a stream in the OBJ file format.
 * @param stream The stream to write to.
 * @param lfMesh The mesh whose vertices are exported.
 * @param walks A range of `Walk`s.
 */
void exportWalksObj(std::ostream& stream, lf::mesh::Mesh const& lfMesh, InputRangeOf<std::is_convertible, Walk> auto const& walks)
{
	exportVerticesObj(stream, lfMesh);

	const std::size_t vertexCount = lfMesh.NumEntities(2);
	for(Walk const& walk : walks)
	{
		stream << 'l';
		assert(walk.vertices.size() >= 2);
		for(const auto vertexIndex : walk.vertices)
		{
			if(vertexIndex >= vertexCount)
				throw std::runtime_error("exportWalksObj: invalid vertex index: " + std::to_string(vertexIndex) + "; limit is " + std::to_string(vertexCount));
			stream << ' ' << vertexIndex + 1;
		}
		stream << '\n';
	}
}

/**
 * @brief Exports a range of `EdgeVector`s to a stream in the OBJ file format.
 * @throws std::runtime_error if the edges' vertices returned by `lfMesh` are invalid.
 * @param stream The stream to write to.
 * @param lfMesh The mesh whose vertices are exported.
 * @param edgeVectors A range of `EdgeVector`s.
 */
void exportEdgeVectorsObj(std::ostream& stream, lf::mesh::Mesh const& lfMesh, InputRangeOf<std::is_convertible, EdgeVector> auto const& edgeVectors)
{
	exportVerticesObj(stream, lfMesh);

	const std::size_t vertexCount = lfMesh.NumEntities(2);
	for(EdgeVector const& edgeVector : edgeVectors)
	{
		for(std::size_t e = 0; e < edgeVector.size(); ++e)
		{
			const int edgeCoeff = edgeVector[e];
			if(edgeCoeff != 0)
			{
				const auto endpoints = lfMesh.EntityByIndex(1, e)->SubEntities(1);
				if(endpoints.size() != 2)
					throw std::runtime_error("exportEdgeVectorsObj: invalid number of endpoints: expected 2, found " + std::to_string(endpoints.size()));
				unsigned from = lfMesh.Index(*endpoints[0]);
				unsigned to = lfMesh.Index(*endpoints[1]);
				for(const auto check : { from, to })
					if(check >= vertexCount)
						throw std::runtime_error("exportEdgeVectorsObj: invalid vertex index: " + std::to_string(check) + "; limit is " + std::to_string(vertexCount));

				if(edgeCoeff < 0)
					std::swap(from, to);

				stream << "l " << from + 1 << ' ' << to + 1 << '\n';
			}
		}
	}
}

} // namespace projects::cohomology

#endif // !defined(COHOMOLOGY_VISUALIZATION_H_INCLUDED)

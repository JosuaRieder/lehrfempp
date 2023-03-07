/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#ifndef COHOMOLOGY_COHOMOLOGY_H_INCLUDED
#define COHOMOLOGY_COHOMOLOGY_H_INCLUDED

#include <utility.h>

#include <lf/mesh/mesh.h>
#include <Eigen/Dense>
#include <quiver.hpp>

#include <type_traits>
#include <functional>
#include <iterator>
#include <vector>
#include <array>
#include <queue>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <string>
#include <sstream>
#include <tuple>
#include <optional>
#include <stack>
#include <ranges>

/**
 * @brief Contains the core functionality of the cohomology subproject.
 */
namespace projects::cohomology
{

/**
 * @brief An error for the case where different ports are connected.
 * @see Mesh::Mesh
 */
struct ConnectedPortsError : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};

/**
 * @brief An error for the case where the provided mesh is not orientable.
 * @see Mesh::Mesh
 */
struct NonOrientableMeshError : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};

/**
 * @brief An error for the case where a nonexistent edge is referenced.
 * @see EdgeVector::EdgeVector
 */
struct InvalidEdgeError : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};

/**
 * @brief An error for the case where no non-excluded arc can be found.
 * @see constructOutlineGraphMap
 */
struct NoNonExcludedArcError : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};

/**
 * @brief An enum type that describes which kind of port primitives are associated with.
 */
enum PortType
{
	noPort,
	electricPort,
	magneticPort,
};

/**
 * @brief Serializes the given `PortType` enumerator.
 * The function's name is deliberately spelled in snake case in order to be ADL-compatible with `std::to_string`.
 * @param portType A valid object of type `PortType`.
 * @return A pointer to a statically stored C-string containing the serialized text.
 */
char const* to_string(PortType portType)
{
	switch(portType)
	{
	case noPort:
		return "no port";
	case electricPort:
		return "electric port";
	case magneticPort:
		return "magnetic port";
	}
	assert(false);
	return nullptr;
}
/**
 * @brief Serializes a `PortType` object to a stream in plain text.
 */
std::ostream& operator<<(std::ostream& stream, PortType portType)
{
	return stream << to_string(portType);
}

/**
 * @brief An annotation type used to inject the member `.edgeIndex` into the graph's edge class.
 * @see Mesh::vertexGraph
 * @see Mesh::faceGraph
 */
struct EdgeIndex
{
	/**
	 * @brief An unsigned index type used to index graph edges.
	 */
	using Index = std::size_t;
	Index edgeIndex;

	/**
	 * @brief Default constructor that initializes `.edgeIndex` with all bits set (`~Index{}`).
	 */
	[[nodiscard]] constexpr EdgeIndex() noexcept
		: EdgeIndex(~Index{})
	{
	}
	/**
	 * @brief Construct from a given value.
	 */
	[[nodiscard]] constexpr EdgeIndex(Index edgeIndex) noexcept
		: edgeIndex(edgeIndex)
	{
	}

	friend constexpr bool operator==(EdgeIndex const&, EdgeIndex const&) = default;
};

/**
 * @brief An enumerator capturing whether the internal orientation of an edge points from the vertex with the larger index to the vertex with the smaller index or the other way around.
 */
enum EdgeOrientationEnum
{
	largeToSmall = -1, //!< -1: Internal edge orientation points from the vertex with the larger index to the vertex with the smaller index.
	smalltoLarge = 1, //!< 1: Internal edge orientation points from the vertex with the smaller index to the vertex with the larger index.
};

/**
 * @brief An annotation type used to inject the member `.orientation` into the graph's edge class.
 * @see Mesh::vertexGraph
 */
struct EdgeOrientation
{
	EdgeOrientationEnum orientation;

	/**
	 * @brief Construct from a given value. No default initialization possible.
	 */
	constexpr EdgeOrientation(EdgeOrientationEnum orientation) noexcept
		: orientation(orientation)
	{
	}

	/**
	 * @brief Find out whether an edge (from, to) is in agreement with this orientation.
	 * Given an edge with a specified internal orientation (`*this`) and an edge (here represented by the arguments `from` and `to`, the edge's vertices),
	 * `factor` computes whether the given vertex pair is in concurrence or in disagreement with this internal orientation.
	 * @param from The start point of the edge under examination.
	 * @param to The end point of the edge under examination.
	 * @return 1 if the internal orientation and the edge (from, to) agree, -1 otherwise.
	 */
	[[nodiscard]] constexpr int factor(quiver::vertex_index_t from, quiver::vertex_index_t to) const noexcept
	{
		return (orientation == smalltoLarge) == (from < to) ? 1 : -1;
	}
};

/**
 * @brief A type completely capturing a primitive's port affiliation (if any).
 * @see Mesh::vertexGraph
 * @see Mesh::faceGraph
 */
struct PortInformation
{
	/**
	 * @brief A signed type to identify ports.
	 */
	using PortIdType = int;

	/**
	 * @brief Which type of port the primitive is affiliated with (if any).
	 */
	PortType portType;
	/**
	 * @brief The ID of the port.
	 * This is necessary n a given problem, there may be multiple electric or multiple magnetic ports.
	 */
	PortIdType id;

	/**
	 * @brief Default constructor. Equivalent to `PortInformation(noPort, 0)`.
	 */
	[[nodiscard]] constexpr PortInformation() noexcept
		: PortInformation(noPort, 0)
	{
	}
	/**
	 * @brief Construct from a given type and ID.
	 */
	[[nodiscard]] constexpr PortInformation(PortType portType, PortIdType id) noexcept
		: portType(portType), id(id)
	{
		assert(portType != noPort || id == 0);
	}

	/**
	 * @brief Returns whether the information in this object indicate that the primitive is affiliated with a port.
	 * @return `true` iff `this->portType != noPort`.
	 */
	[[nodiscard]] constexpr bool isAffiliated() const noexcept
	{
		return portType != noPort;
	}
	/**
	 * @brief Returns whether the information in this object indicate that the primitive is unaffiliated with a port.
	 * @return `true` iff `this->portType == noPort`.
	 */
	[[nodiscard]] constexpr bool isUnaffiliated() const noexcept
	{
		return !isAffiliated();
	}

	friend constexpr bool operator==(PortInformation const&, PortInformation const&) = default;
};
/**
 * @brief Serializes a `PortInformation` object to a stream in plain text.
 * The output text is formatted as follows: portType#id where the #id part is only displayed in affiliated situations.
 */
std::ostream& operator<<(std::ostream& stream, PortInformation portInformation)
{
	stream << portInformation.portType;
	if(portInformation.isAffiliated())
		stream << '#' << portInformation.id;
	return stream;
}
/**
 * @brief Serializes a `PortInformation` object to a string in plain text.
 * The output text is formatted as follows: portType#id where the #id part is only displayed in affiliated situations.
 */
std::string to_string(PortInformation portInformation)
{
	std::ostringstream oss;
	oss << portInformation;
	return std::move(oss).str();
}

} // namespace projects::cohomology

namespace std
{
	/**
	 * @brief A template specialization of `std::hash` to work with objects of type `PortInformation`.
	 */
	template<>
	struct hash<projects::cohomology::PortInformation>
	{
		size_t operator()(projects::cohomology::PortInformation const& portInformation) const
		{
			return hashTogether(portInformation.portType, portInformation.id);
		}
	};
}

namespace projects::cohomology
{
/**
 * @brief A concept for types that are invocable with `quiver::vertex_index_t` and return `PortInformation`.
 * Such mappings are used to communicate the affiliations of faces with ports.
 * @see Mesh::Mesh
 */
template<typename Invocable>
concept PortInformationInvocable = std::is_invocable_r_v<PortInformation, Invocable, quiver::vertex_index_t>;

/**
 * @brief A constant mapping from any index to a `PortInformation` object that indicates no affiliation.
 * Satisfies `PortInformationInvocable`.
 * @see PortInformationInvocable
 * @see Mesh::Mesh
 */
[[nodiscard]] constexpr PortInformation noPorts(quiver::vertex_index_t face) noexcept
{
	(void)face; // -Wunused-parameter
	return {};
}
static_assert(PortInformationInvocable<decltype(noPorts)>);

/**
 * @brief An enumerator type that captures the relation of two faces with respect to their relative internal orientation.
 * The orientation of interest is the one that is induced on the common edge by the faces' internal orientation.
 */
enum class RelativeFaceOrientation
{
	vertexDisjoint, //!< The two faces don't share any vertices.
	edgeDisjoint, //!< The two faces share a vertex but no edges.
	inAgreement, //!< The two faces share an edge and are in agreement with regard to the orientation induced on the shared edge.
	inDisagreement, //!< The two faces share an edge and are in disagreement with regard to the orientation induced on the shared edge.
};

/**
 * @brief An annotation type used to inject the member `.vertexIndices` into the face graph's vertex class.
 * @see Mesh::faceGraph
 */
struct FaceVertices
{
	/**
	 * @brief The vertices belonging to a face in the order induced by the internal orientation of the face.
	 */
	std::vector<quiver::vertex_index_t> vertexIndices;

	/**
	 * @brief Default constructor, leaves `.vertexIndices` empty.
	 */
	[[nodiscard]] FaceVertices() noexcept
	{
	}

	/**
	 * @brief Returns whether a given vertex belongs to the face.
	 * @param vertex The index of the vertex to be checked.
	 * @return `true` iff the vertex with index `vertex` belongs to this face, i.e. whether it occurs in `.vertexIndices`.
	 */
	[[nodiscard]] bool contains(quiver::vertex_index_t vertex) const noexcept
	{
		return std::find(vertexIndices.begin(), vertexIndices.end(), vertex) != vertexIndices.end();
	}

	/**
	 * @brief Returns the index within `.vertexIndices` of a given vertex index.
	 * @param vertex The index of the vertex to be found.
	 * @return The index of `vertex` within `.vertexIndices` if it exists, `.vertexIndices.size()` otherwise.
	 */
	[[nodiscard]] std::size_t indexOf(quiver::vertex_index_t vertex) const noexcept
	{
		return std::distance(vertexIndices.begin(), std::find(vertexIndices.begin(), vertexIndices.end(), vertex));
	}

	/**
	 * @brief For a given edge (represented by its vertices `from` and `to`) that is known to belong to this face, returns whether the orientations agree.
	 * In other words, it checks which one of the sequences `(from, to)` or `(to, from)` occurs in the face's vertex list (modulo wrap-around).
	 * @param from The start point of the edge under examination.
	 * @param to The end point of the edge under examination.
	 * @return 1 if the internal orientation of the face induces an orientation on the edge that points from `from` to `to`, -1 otherwise.
	 */
	[[nodiscard]] int orderOf(quiver::vertex_index_t from, quiver::vertex_index_t to) const noexcept
	{
		assert(contains(from));
		assert(contains(to));
		assert(from != to);
		const std::size_t fromIndex = indexOf(from);
		const std::size_t nextIndex = (fromIndex + 1) % vertexIndices.size();
		assert(vertexIndices[nextIndex] == to || vertexIndices[(fromIndex + vertexIndices.size() - 1) % vertexIndices.size()] == to);
		return vertexIndices[nextIndex] == to ? 1 : -1;
	}

	/**
	 * @brief For two faces' vertex sequences, finds the first shared vertex and returns its indices (within the sequences).
	 * @param other The other face's vertex sequence.
	 * @return `std::nullopt` if `*this` and `other` share no vertices, otherwise a pair of indices, the first being the index of the shared vertex in `this->vertexIndices` and the second being its index in `other.vertexIndices`.
	 */
	[[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>> getSharedVertex(FaceVertices const& other) const noexcept
	{
		// This algorithm is O(n*m) but probably the fastest for n,m <= 4.
		// TODO: For larger vertex sets, a O(max(n,m)) algorithm may be implemented.
		for(std::size_t i = 0; i < vertexIndices.size(); ++i)
			for(std::size_t j = 0; j < other.vertexIndices.size(); ++j)
				if(vertexIndices[i] == other.vertexIndices[j])
					return std::make_pair(i, j);
		return std::nullopt;
	}

	/**
	 * @brief Categorizes two faces according to `RelativeFaceOrientation`.
	 * @param other The other face's vertex sequence.
	 * @return The appropriate value of `RelativeFaceOrientation` describing the relationship between `*this` and `other`.
	 */
	[[nodiscard]] RelativeFaceOrientation relativeOrientation(FaceVertices const& other) const noexcept
	{
		const auto indices = getSharedVertex(other);
		if(!indices)
			return RelativeFaceOrientation::vertexDisjoint;
		const auto [thisFirst, otherFirst] = *indices;
		const std::size_t thisPrev = (thisFirst + vertexIndices.size() - 1) % vertexIndices.size();
		const std::size_t thisNext = (thisFirst + 1) % vertexIndices.size();
		const std::size_t otherPrev = (otherFirst + other.vertexIndices.size() - 1) % other.vertexIndices.size();
		const std::size_t otherNext = (otherFirst + 1) % other.vertexIndices.size();
		if(vertexIndices[thisNext] == other.vertexIndices[otherNext] || vertexIndices[thisPrev] == other.vertexIndices[otherPrev])
			return RelativeFaceOrientation::inAgreement;
		if(vertexIndices[thisNext] == other.vertexIndices[otherPrev] || vertexIndices[thisPrev] == other.vertexIndices[otherNext])
			return RelativeFaceOrientation::inDisagreement;
		return RelativeFaceOrientation::edgeDisjoint;
	}

	/**
	 * @brief For two faces that share an edge, this function reverses `other`'s vertices iff the two faces are in agreement with respect to the shared edge (see `RelativeFaceOrientation::inAgreement`).
	 * After this function has been called, `relativeOrientation(other)` will always return `RelativeFaceOrientation::inDisagreement`.
	 * @see RelativeFaceOrientation::inAgreement
	 * @param other The other face's vertex sequence.
	 * @return Whether the vertices of `other` have been reversed.
	 */
	bool orientDisagreeably(FaceVertices& other) const noexcept
	{
		const RelativeFaceOrientation relative = relativeOrientation(other);
		assert(relative == RelativeFaceOrientation::inAgreement || relative == RelativeFaceOrientation::inDisagreement);
		if(relative == RelativeFaceOrientation::inAgreement)
			std::reverse(other.vertexIndices.begin(), other.vertexIndices.end());
		return relative == RelativeFaceOrientation::inAgreement;
	}

	friend constexpr bool operator==(FaceVertices const&, FaceVertices const&) = default;
};

/**
 * @brief The information stored per edge.
 * @see Mesh::edgeConnectivity
 */
struct EdgeInformation
{
	/**
	 * @brief A placeholder value for the entries of `.vertices` and `.faces` if they have not been set yet.
	 */
	static constexpr quiver::vertex_index_t notSet = ~0;

	/**
	 * @brief The indices of the two vertices contained in this edge.
	 */
	std::array<quiver::vertex_index_t, 2> vertices = { notSet, notSet };
	/**
	 * @brief The indices of the two faces that contain this edge.
	 */
	std::array<quiver::vertex_index_t, 2> faces = { notSet, notSet };

	/**
	 * @brief Default constructor that sets all elements of `.vertices` and `.faces` to `.notSet`.
	 */
	[[nodiscard]] constexpr EdgeInformation() noexcept
	{
	}

	/**
	 * @brief Returns whether two faces have been set for this edge.
	 * @return `true` iff two faces have been set for this edge.
	 */
	[[nodiscard]] constexpr bool hasTwoFaces() const noexcept
	{
		return faces[1] != notSet;
	}
	/**
	 * @brief If less than two faces have been set, this member function adds another face to `.faces`.
	 * @param face The index of the face to be added.
	 */
	constexpr void addFace(quiver::vertex_index_t face) noexcept
	{
		assert(!hasTwoFaces());
		faces[faces[0] == notSet ? 0 : 1] = face;
	}

	/**
	 * @brief Returns whether this edge contains a specific vertex.
	 * @param vertex The index of the vertex in question.
	 * @return `true` iff this edge contains `vertex`.
	 */
	[[nodiscard]] constexpr bool containsVertex(quiver::vertex_index_t vertex) const noexcept
	{
		return vertices[0] == vertex || vertices[1] == vertex;
	}
	/**
	 * @brief Returns whether this edge is contained by a specific face.
	 * @param face The index of the face in question.
	 * @return `true` iff this edge contains `face`.
	 */
	[[nodiscard]] constexpr bool containsFace(quiver::vertex_index_t face) const noexcept
	{
		return faces[0] == face || faces[1] == face;
	}

	/**
	 * @brief Returns whether this edge shares a vertex with another edge.
	 * @param edge The other edge.
	 * @return `true` iff this edge `*this` and the other edge `edge` have at least one vertex in common.
	 */
	[[nodiscard]] constexpr bool sharesVertex(EdgeInformation const& edge) const noexcept
	{
		return containsVertex(edge.vertices[0]) || containsVertex(edge.vertices[1]);
	}
	/**
	 * @brief Returns whether this edge shares a face with another edge.
	 * @param edge The other edge.
	 * @return `true` iff this edge `*this` and the other edge `edge` have at least one face in common.
	 */
	[[nodiscard]] constexpr bool sharesFace(EdgeInformation const& edge) const noexcept
	{
		return containsFace(edge.faces[0]) || containsFace(edge.faces[1]);
	}

	friend constexpr bool operator==(EdgeInformation const&, EdgeInformation const&) = default;
};

/**
 * @brief A product type to capture the number of ports of each type.
 * @see PortType
 */
struct PortCount
{
	/**
	 * @brief The number of electric ports.
	 */
	std::size_t electricPorts = 0;
	/**
	 * @brief The number of magnetic ports.
	 */
	std::size_t magneticPorts = 0;

	friend constexpr bool operator==(PortCount const&, PortCount const&) = default;
};

/**
 * @brief The main data type that wholly captures the problem's input.
 * It additionally contains precomputed data (refer to the thesis).
 * After an instance of this class has been created, the data structures used to construct it don't have to be preserved.
 */
struct Mesh
{
	/**
	 * @brief The vertex graph, implemented as an undirected graph. Its vertex set represents the mesh vertices, its edge set represents the mesh edges.
	 * The vertices (i.e. mesh vertices) are annotated with `PortInformation`.
	 * The edges (i.e. mesh edges) are annotated with `EdgeIndex` and `EdgeOrientation`.
	 */
	quiver::adjacency_list<quiver::undirected, quiver::cmb<EdgeIndex, EdgeOrientation>, PortInformation> vertexGraph; // V = mesh vertices, E = mesh edges
	/**
	 * @brief The face graph, implemented as an undirected graph. Its vertex set represents the mesh faces, its edge set represents the mesh face adjacencies.
	 * The vertices (i.e. mesh faces) are annotated with `PortInformation` and `FaceVertices`.
	 * The edges (i.e. mesh face adjacencies) are annotated with `EdgeIndex`. This annotation (which is attached to face adjacencies) captures the index of the edge that lies between the two adjacent faces.
	 */
	quiver::adjacency_list<quiver::undirected, EdgeIndex, quiver::cmb<PortInformation, FaceVertices>> faceGraph; // V = mesh faces, E = face adjacency
	/**
	 * @brief Auxiliary information stored for each (mesh) edge.
	 * @see EdgeInformation
	 */
	std::vector<EdgeInformation> edgeConnectivity;

	/**
	 * @brief Tries to establish a topological orientation on the mesh by selectively reversing the order of vertices of the faces.
	 * @return `true` iff either the mesh has already been oriented or the mesh has successfully been oriented.
	 * @see FaceVertices::orientDisagreeably
	 */
	bool tryOrient()
	{
		// DFS might be preferable over BFS here in the hopes that it may detect unorientability more quickly.
		// "Worst case", they both take O(|F|) time.

		if(faceGraph.V.empty())
			return true;

		std::vector<bool> oriented(faceGraph.V.size(), false);
		std::size_t left = faceGraph.V.size();
		const auto nowOriented = [&](quiver::vertex_index_t i)
		{
			assert(oriented[i] == false);
			oriented[i] = true;
			--left;
		};

		std::vector<bool> enqueued(faceGraph.V.size(), false);
		std::stack<std::pair<quiver::vertex_index_t, quiver::vertex_index_t>> neighbors; // .first = oriented face, .second = unoriented, adjacent face
		while(left != 0) // This is only needed if the mesh consists of multiple connected components.
		{
			const quiver::vertex_index_t start = std::distance(oriented.begin(), std::find(oriented.begin(), oriented.end(), false));
			nowOriented(start);

			// Start of manual DFS
			enqueued[start] = true;
			neighbors.emplace(start, start);
			do
			{
				const auto [predIndex, index] = neighbors.top();
				neighbors.pop();
				assert(enqueued[index]);
				auto& vertex = faceGraph.V[index];
				auto const& predVertex = faceGraph.V[predIndex];

				if(index != predIndex) // see above .emplace(start, start)
				{
					predVertex.orientDisagreeably(vertex);
					nowOriented(index);
				}

				// reverse here so that first inserted is first visited
				const auto begin = vertex.out_edges.rbegin();
				const auto end = vertex.out_edges.rend();
				for(auto iter = begin; iter != end; ++iter)
					if(!enqueued[iter->to])
					{
						neighbors.emplace(index, iter->to);
						enqueued[iter->to] = true;
					}
					else if(oriented[iter->to])
					{
						if(vertex.relativeOrientation(faceGraph.V[iter->to]) != RelativeFaceOrientation::inDisagreement)
							return false;
					}
			}
			while(!neighbors.empty());
			// End of manual DFS
		}
		return true;
	}

	/**
	 * @brief Default constructor that leaves all data structures empty.
	 */
	[[nodiscard]] Mesh()
	{
	}
	/**
	 * @brief Initializes from a mesh and assumes no ports. Equivalent to `Mesh(mesh, noPorts)`. Performs all precomputations.
	 * @throws ConnectedPortsError iff a vertex is shared among faces with different port associations; iff faces of different ports touch.
	 * @throws NonOrientableMeshError iff the mesh is not orientable.
	 * @see noPorts
	 */
	[[nodiscard]] Mesh(lf::mesh::Mesh const& mesh)
		: Mesh(mesh, noPorts)
	{
	}
	/**
	 * @brief Initializes from a mesh and a mapping that indicates port associations. Performs all precomputations.
	 * @throws ConnectedPortsError iff a vertex is shared among faces with different port associations; iff faces of different ports touch.
	 * @throws NonOrientableMeshError iff the mesh is not orientable.
	 * @param mesh The mesh of the problem.
	 * @param portInformationInvocable A mapping from face indices to `PortInformation` to indicate port associations.
	 * @see PortInformationInvocable
	 */
	[[nodiscard]] Mesh(lf::mesh::Mesh const& mesh, PortInformationInvocable auto&& portInformationInvocable)
		: vertexGraph(mesh.NumEntities(2)), faceGraph(mesh.NumEntities(0)), edgeConnectivity(mesh.NumEntities(1)) // number = codimension
	{
		using std::begin, std::end, std::to_string;

		assert(mesh.DimMesh() == 2);

		// Construct all the edges for the vertex graph.
		// Fill in the vertex ownership data to the edge list.
		for(lf::mesh::Entity const* const edge : mesh.Entities(1))
		{
			auto const& endpoints = edge->SubEntities(1);

			assert(endpoints.size() == 2);
			assert(std::all_of(begin(endpoints), end(endpoints), [](auto* endpoint){ return endpoint->Codim() == 2; }));

			std::array<quiver::vertex_index_t, 2> vertices;
			for(int i = 0; i < 2; ++i)
				vertices[i] = mesh.Index(*endpoints[i]);
			const EdgeIndex::Index edgeIndex = mesh.Index(*edge);

			assert(std::all_of(begin(vertices), end(vertices), [this](quiver::vertex_index_t vertexIndex){ return vertexIndex < this->vertexGraph.V.size(); }));
			assert(edgeIndex < edgeConnectivity.size());

			vertexGraph.E.emplace(vertices[0], vertices[1], edgeIndex, vertices[0] < vertices[1] ? smalltoLarge : largeToSmall);
			edgeConnectivity[edgeIndex].vertices = vertices;
		}

		// Complete the edge list (add face ownership information)
		for(lf::mesh::Entity const* const face : mesh.Entities(0))
			for(lf::mesh::Entity const* const edge : face->SubEntities(1))
				edgeConnectivity[mesh.Index(*edge)].addFace(mesh.Index(*face));

		// Assign all the properties to the faces, i.e. port information and vertex indices
		// Assign port information to the vertices.
		for(quiver::vertex_index_t v = 0; v < faceGraph.V.size(); ++v)
		{
			const PortInformation portInformation = portInformationInvocable(v);
			faceGraph.V[v] = portInformation;

			lf::mesh::Entity const* const face = mesh.EntityByIndex(0, v);
			auto const& vertices = face->SubEntities(2);
			assert(vertices.size() >= 3);
			faceGraph.V[v].vertexIndices.reserve(vertices.size());
			for(lf::mesh::Entity const* const vertex : vertices)
			{
				const quiver::vertex_index_t vertexIndex = mesh.Index(*vertex);
				faceGraph.V[v].vertexIndices.emplace_back(vertexIndex);

				// Check that every vertex belongs to only one port.
				PortInformation& vertexPortInformation = vertexGraph.V[vertexIndex];
				if(portInformation.portType != noPort && portInformation != vertexPortInformation)
				{
					if(vertexGraph.V[vertexIndex].portType != noPort)
						throw ConnectedPortsError("Mesh::Mesh: vertex " + to_string(vertexIndex)
							+ " is shared among ports " + to_string(vertexPortInformation)
							+ " and " + to_string(portInformation));
					vertexPortInformation = portInformation;
				}
			}
		}

		// Construct the face graph
		for(EdgeIndex::Index e = 0; e < edgeConnectivity.size(); ++e)
			if(edgeConnectivity[e].hasTwoFaces())
				faceGraph.E.emplace(edgeConnectivity[e].faces[0], edgeConnectivity[e].faces[1], e);

		// Calculate an orientation of the mesh
		if(!tryOrient())
			throw NonOrientableMeshError("Mesh::Mesh: mesh is not orientable");
	}

	/**
	 * @brief Returns the number of vertices of the mesh.
	 * @return The number of vertices of the mesh.
	 */
	[[nodiscard]] constexpr std::size_t vertices() const noexcept
	{
		return vertexGraph.V.size();
	}
	/**
	 * @brief Returns the number of edges of the mesh.
	 * @return The number of edges of the mesh.
	 */
	[[nodiscard]] constexpr std::size_t edges() const noexcept
	{
		return vertexGraph.E.size();
	}
	/**
	 * @brief Returns the number of faces of the mesh.
	 * @return The number of faces of the mesh.
	 */
	[[nodiscard]] constexpr std::size_t faces() const noexcept
	{
		return faceGraph.V.size();
	}

	/**
	 * @brief Computes the number of distinct ports of each type.
	 * @return A filled-out object of type `PortCount`.
	 * @see PortCount
	 */
	[[nodiscard]] PortCount countPorts() const noexcept
	{
		std::unordered_set<int> ports[2];
		for(auto const& face : faceGraph.V)
			if(face.portType != noPort)
				ports[face.portType == electricPort ? 0 : 1].emplace(face.id);
		return { ports[0].size(), ports[1].size() };
	}

	/**
	 * @brief Returns which type of port a vertex is affiliated with (if any).
	 * @param vertexIndex The index of the vertex in question.
	 * @return `std::nullopt` if the vertex does not belong to a port, the type of the port it belongs to otherwise.
	 * @see PortType
	 */
	[[nodiscard]] std::optional<PortInformation> isPortVertex(quiver::vertex_index_t vertexIndex) const noexcept
	{
		auto const& vertex = vertexGraph.V[vertexIndex];
		if(vertex.isAffiliated())
			return static_cast<PortInformation const&>(vertex);
		else
			return {};
	}
	/**
	 * @brief Returns which type of port a vertex is affiliated with if the vertex is on a port boundary.
	 * @param vertexIndex The index of the vertex in question.
	 * @return `std::nullopt` if the vertex does not belong to a port or is not located on the port's boundary, the type of the port it belongs to otherwise.
	 * @see PortType
	 */
	[[nodiscard]] std::optional<PortInformation> isPortBoundaryVertex(quiver::vertex_index_t vertexIndex) const noexcept
	{
		auto const& vertex = vertexGraph.V[vertexIndex];
		auto const& pi = static_cast<PortInformation const&>(vertex);
		if(vertex.isAffiliated())
		{
			using std::begin, std::end;
			for(auto const& outEdge : vertex.out_edges)
				if(pi != vertexGraph.V[outEdge.to])
					return pi;
		}
		return {};
	}
	/**
	 * @brief Returns which type of port a vertex is affiliated with if the vertex is in the interior of a port.
	 * @param vertexIndex The index of the vertex in question.
	 * @return `std::nullopt` if the vertex is not in the interior of a port, the type of the port it belongs to otherwise.
	 * @see PortType
	 */
	[[nodiscard]] std::optional<PortInformation> isPortInteriorVertex(quiver::vertex_index_t vertexIndex) const noexcept
	{
		if(std::optional<PortInformation> result = isPortVertex(vertexIndex))
			if(!isPortBoundaryVertex(vertexIndex))
				return result;
		return {};
	}
	/**
	 * @brief Returns whether a vertex is excluded (see the thesis for an explanation of excluded primitives).
	 * @param vertexIndex The index of the vertex in question.
	 * @return `true` iff the vertex with index `vertexIndex` is an excluded vertex.
	 * @see isExcludedEdge
	 * @see isExcludedFace
	 */
	[[nodiscard]] bool isExcludedVertex(quiver::vertex_index_t vertexIndex) const noexcept
	{
		const auto isPort = isPortVertex(vertexIndex);
		return isPort && (isPort->portType == electricPort || isPortInteriorVertex(vertexIndex));
	}

	/**
	 * @brief Returns how many of the faces the edge is contained in are affiliated with a port and, if any, what type of port it is.
	 * @param edgeIndex The index of the edge in question.
	 * @return A pair consisting of the number of affiliated faces and their `PortInformation`.
	 * @see PortInformation::isAffiliated
	 */
	[[nodiscard]] std::pair<unsigned, PortInformation> affiliatedPortOfEdge(EdgeIndex::Index edgeIndex) const noexcept
	{
		using std::begin, std::end;
		PortInformation portInformation;
		auto const& faces = edgeConnectivity[edgeIndex].faces;
		const unsigned affiliatedPorts = std::count_if(begin(faces), end(faces), [this, &portInformation](quiver::vertex_index_t faceIndex)
		{
			auto const& face = this->faceGraph.V[faceIndex];
			if(portInformation.portType == noPort)
				portInformation = face.properties();
			return face.isAffiliated();
		});
		assert(in(affiliatedPorts, 0, 2));
		return { affiliatedPorts, portInformation };
	}
	/**
	 * @brief Returns which type of port an edge is affiliated with if the edge is on a port boundary.
	 * @param edgeIndex The index of the edge in question.
	 * @return `std::nullopt` if the edge does not belong to a port or is not located on the port's boundary, the port information it belongs to otherwise.
	 * @see PortInformation
	 */
	[[nodiscard]] std::optional<PortInformation> isPortBoundaryEdge(EdgeIndex::Index edgeIndex) const noexcept
	{
		const auto [count, portInformation] = affiliatedPortOfEdge(edgeIndex);
		if(count == 1)
			return portInformation;
		else
			return std::nullopt;
	}
	/**
	 * @brief Returns whether an edge is on the boundary of an electrical port.
	 * @param edgeIndex The index of the edge in question.
	 * @return `true` iff the edge with index `edgeIndex` lies on the boundary of an electrical port, i.e. between an affiliated, electrical face and an unaffiliated face.
	 * @see PortInformation::isAffiliated
	 */
	[[nodiscard]] bool isElectricalBoundaryEdge(EdgeIndex::Index edgeIndex) const noexcept
	{
		const auto [count, portInformation] = affiliatedPortOfEdge(edgeIndex);
		return portInformation.portType == electricPort && count == 1;
	}
	/**
	 * @brief Returns whether an edge is excluded (see the thesis for an explanation of excluded primitives).
	 * @param edgeIndex The index of the edge in question.
	 * @return `true` iff the edge with index `edgeIndex` is an excluded edge.
	 * @see isExcludedVertex
	 * @see isExcludedFace
	 */
	[[nodiscard]] bool isExcludedEdge(EdgeIndex::Index edgeIndex) const noexcept
	{
		const auto [count, portInformation] = affiliatedPortOfEdge(edgeIndex);
		switch(portInformation.portType)
		{
		case electricPort:
			assert(in(count, 1, 2));
			return true;
		case magneticPort:
			assert(in(count, 1, 2));
			return count == 2;
		default:
			return false;
		}
	}

	/**
	 * @brief Returns whether a face is excluded (see the thesis for an explanation of excluded primitives).
	 * @param faceIndex The index of the face in question.
	 * @return `true` iff the face with index `faceIndex` is an excluded face.
	 * @see isExcludedVertex
	 * @see isExcludedEdge
	 */
	[[nodiscard]] bool isExcludedFace(quiver::vertex_index_t faceIndex) const noexcept
	{
		return faceGraph.V[faceIndex].isAffiliated();
	}

	/**
	 * @brief Constructs the edge-cell incidence matrix where the edges correspond to columns and the faces to rows.
	 * @return The edge-cell incidence matrix, with dimensions `.rows() == mesh.faces()` and `.cols() == mesh.edges()`.
	 */
	[[nodiscard]] Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> getEdgeCellIncidence() const
	{
		using Matrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
		Matrix C = Matrix::Zero(faces(), edges());
		std::size_t matrixFaceEntry = 0;
		for(quiver::vertex_index_t faceIndex = 0; faceIndex < faceGraph.V.size(); ++faceIndex)
		{
			assert(faceGraph.V[faceIndex].vertexIndices.size() >= 3);
			auto first = faceGraph.V[faceIndex].vertexIndices.begin();
			const auto last = faceGraph.V[faceIndex].vertexIndices.end();
			const unsigned front = *first;
			unsigned current = front;
			do
			{
				unsigned next;
				if(++first == last)
					next = front;
				else
					next = *first;

				const auto outEdge = vertexGraph.E.get(current, next);
				assert(outEdge != nullptr);
				const int direction = outEdge->factor(current, next);
				C(matrixFaceEntry, outEdge->edgeIndex) = direction;

				current = next;
			}
			while(first != last);

			++matrixFaceEntry;
		}
		return C;
	}
	/**
	 * @brief Constructs the node-edge incidence matrix where the vertices correspond to columns and the edges to rows.
	 * @return The node-edge incidence matrix, with dimensions `.rows() == mesh.edges()` and `.cols() == mesh.vertices()`.
	 */
	[[nodiscard]] Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> getNodeEdgeIncidence() const
	{
		using Matrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
		Matrix G = Matrix::Zero(edges(), vertices());
		for(quiver::vertex_index_t vertex = 0; vertex < vertexGraph.V.size(); ++vertex)
			for(auto& outEdge : vertexGraph.V[vertex].out_edges)
				if(outEdge.factor(vertex, outEdge.to) == 1)
				{
					G(outEdge.edgeIndex, vertex) = -1;
					G(outEdge.edgeIndex, outEdge.to) = 1;
				}

		return G;
	}
};

/**
 * @brief Removes a row from a matrix.
 * @param matrix The matrix to remove the row from.
 * @param rowToRemove The index of the row to remove.
 */
void removeRow(auto& matrix, std::size_t rowToRemove)
{
	const std::size_t rows = matrix.rows() - 1;
	const std::size_t cols = matrix.cols();

	assert(rowToRemove <= rows);

	matrix.block(rowToRemove, 0, rows - rowToRemove, cols) = matrix.block(rowToRemove + 1, 0, rows - rowToRemove, cols).eval();
	matrix.conservativeResize(rows, cols);
}

/**
 * @brief Removes a column from a matrix.
 * @param matrix The matrix to remove the column from.
 * @param colToRemove The index of the column to remove.
 */
void removeColumn(auto& matrix, std::size_t colToRemove)
{
	const std::size_t rows = matrix.rows();
	const std::size_t cols = matrix.cols() - 1;

	assert(colToRemove <= cols);

	matrix.block(0, colToRemove, rows, cols - colToRemove) = matrix.block(0, colToRemove + 1, rows, cols - colToRemove).eval();
	matrix.conservativeResize(rows, cols);
}

// TODO: This function can be massively improved (performance-wise) using Eigen's ArithmeticSequence but the version of Eigen available in LehrFEM does not support this feature.
/**
 * @brief Removes slices from a matrix according to whether these slices correspond to excluded primitives.
 * The slices along the axis `matrixAxis` are interpreted as corresponding to primitives of dimension `primitiveDimensions`.
 * The matrix endows these slices with a natural index which is interpreted as the respective primitive's index.
 * This primitive index is then used to query whether the primitive in question is excluded. If it is, the slice is removed from the matrix.
 * The typical use case is to remove rows and columns of incidence matrices.
 * @tparam primitiveDimensions The number of dimensions of the primitives in question (0 = vertices, 1 = edges, 2 = faces).
 * @tparam matrixAxis Along which axis the slices are to be taken.
 * @param mesh The mesh whose `Mesh::isExcludedVertex`, `Mesh::isExcludedEdge` and `Mesh::isExcludedFace` methods are invoked.
 * @param matrix The matrix whose slices are to be removed, typically an incidence matrix generated by `Mesh::getEdgeCellIncidence` or `Mesh::getNodeEdgeIncidence`.
 * @return The number of removed slices.
 * @see Mesh::isExcludedVertex
 * @see Mesh::isExcludedEdge
 * @see Mesh::isExcludedFace
 */
template<std::size_t primitiveDimensions, std::size_t matrixAxis>
std::size_t removeExcludedSlices(Mesh const& mesh, auto& matrix)
{
	static_assert(primitiveDimensions < 3);
	static_assert(matrixAxis < 2);

	using Matrix = std::decay_t<decltype(matrix)>;
	constexpr auto primitiveCount = std::array<std::size_t (Mesh::*)() const noexcept, 3>{ &Mesh::vertices, &Mesh::edges, &Mesh::faces }[primitiveDimensions];
	constexpr auto isExcluded = std::array<bool (Mesh::*)(quiver::vertex_index_t) const noexcept, 3>{ &Mesh::isExcludedVertex, &Mesh::isExcludedEdge, &Mesh::isExcludedFace }[primitiveDimensions];
	constexpr auto sliceSize = std::array<Eigen::Index (Matrix::*)() const, 2>{ &Matrix::rows, &Matrix::cols }[matrixAxis];
	constexpr auto sliceDelete = std::array<void (*)(Matrix&, std::size_t), 2>{ &removeRow<Matrix>, &removeColumn<Matrix> }[matrixAxis];

	const auto primitives = (mesh.*primitiveCount)();
	const auto size = (matrix.*sliceSize)();
	assert(primitives == size);

	std::size_t matrixIndex = 0;
	for(std::size_t primitiveIndex = 0; primitiveIndex < size; ++primitiveIndex)
	{
		if((mesh.*isExcluded)(primitiveIndex))
			sliceDelete(matrix, matrixIndex);
		else
			++matrixIndex;
	}
	return size - matrixIndex;
}

/**
 * @brief A walk on a graph represented by a list of vertices.
 * The edge `{.vertices[0], .vertices[.size()-1]}` is conceptually not thought of as being included.
 * This means that, in order to represent closed walks, loops, the last vertex of the walk should be identical to the first vertex.
 * If one of the edges purported to exist by their inclusion in the walk does not exist, the behavior is undefined.
 */
class Walk
{
	/**
	 * @brief An iterator type that returns, among others, a reference to the corresponding `out_edge_t`.
	 * Using this iterator makes it so that, conceptually, one operates on edges even though the walk is represented by vertices.
	 * The value type is `const std::tuple<quiver::vertex_index_t, quiver::vertex_index_t, decltype(Mesh::vertexGraph)::out_edge_t const&>`.
	 */
	class ConstEdgeIterator
	{
		Walk const* walk = nullptr;
		std::ptrdiff_t i = 0;

		[[nodiscard]] bool isValid() const noexcept
		{
			return walk == nullptr || (0 <= i && i <= walk->edges());
		}
		[[nodiscard]] bool isValidNonsingular() const noexcept
		{
			return walk != nullptr && (0 <= i && i <= walk->edges());
		}
		[[nodiscard]] bool isDereferenceable() const noexcept
		{
			return walk != nullptr && (0 <= i && i < walk->edges());
		}

	public:
		using difference_type = std::ptrdiff_t;
		using value_type = const std::tuple<quiver::vertex_index_t, quiver::vertex_index_t, decltype(Mesh::vertexGraph)::out_edge_t const&>;
		using pointer = value_type*;
		using reference = value_type;
		using iterator_category = std::random_access_iterator_tag;

		[[nodiscard]] ConstEdgeIterator() noexcept
		{
		}
		[[nodiscard]] ConstEdgeIterator(Walk const* walk, std::ptrdiff_t i) noexcept
			: walk(walk), i(i)
		{
			assert(isValid());
		}

		value_type operator*() const noexcept
		{
			assert(isDereferenceable());
			const quiver::vertex_index_t from = walk->vertices[i];
			const quiver::vertex_index_t to = walk->vertices[i + 1];
			const auto outEdge = walk->mesh->vertexGraph.E(from, to);
			assert(outEdge != nullptr);
			return { from, to, *outEdge };
		}

		ConstEdgeIterator& operator++() noexcept
		{
			assert(isValidNonsingular());
			++i;
			assert(isValidNonsingular());
			return *this;
		}
		ConstEdgeIterator operator++(int) noexcept
		{
			ConstEdgeIterator result = *this;
			operator++();
			return result;
		}
		ConstEdgeIterator& operator--() noexcept
		{
			assert(isValidNonsingular());
			--i;
			assert(isValidNonsingular());
			return *this;
		}
		ConstEdgeIterator operator--(int) noexcept
		{
			ConstEdgeIterator result = *this;
			operator--();
			return result;
		}

		ConstEdgeIterator& operator+=(difference_type n) noexcept
		{
			assert(isValidNonsingular());
			i += n;
			assert(isValidNonsingular());
			return *this;
		}
		ConstEdgeIterator& operator-=(difference_type n) noexcept
		{
			assert(isValidNonsingular());
			i -= n;
			assert(isValidNonsingular());
			return *this;
		}

		friend /*[[nodiscard]]*/ ConstEdgeIterator operator+(ConstEdgeIterator iter, difference_type n) noexcept
		{
			return iter += n;
		}
		friend /*[[nodiscard]]*/ ConstEdgeIterator operator+(difference_type n, ConstEdgeIterator iter) noexcept
		{
			return iter += n;
		}
		friend /*[[nodiscard]]*/ difference_type operator-(ConstEdgeIterator const& lhs, ConstEdgeIterator const& rhs) noexcept
		{
			assert(lhs.walk == rhs.walk);
			return lhs.i - rhs.i;
		}
		friend /*[[nodiscard]]*/ ConstEdgeIterator operator-(ConstEdgeIterator lhs, difference_type n) noexcept
		{
			return lhs -= n;
		}

		value_type operator[](difference_type n) const noexcept
		{
			return *(*this + n);
		}

		[[nodiscard]] bool operator==(ConstEdgeIterator const& rhs) const noexcept
		{
			return walk == rhs.walk && i == rhs.i;
		}
		[[nodiscard]] bool operator!=(ConstEdgeIterator const& rhs) const noexcept
		{
			return !operator==(rhs);
		}
		[[nodiscard]] bool operator<(ConstEdgeIterator const& rhs) const noexcept
		{
			return walk == rhs.walk && i < rhs.i;
		}
		[[nodiscard]] bool operator<=(ConstEdgeIterator const& rhs) const noexcept
		{
			return walk == rhs.walk && i <= rhs.i;
		}
		[[nodiscard]] bool operator>(ConstEdgeIterator const& rhs) const noexcept
		{
			return walk == rhs.walk && i > rhs.i;
		}
		[[nodiscard]] bool operator>=(ConstEdgeIterator const& rhs) const noexcept
		{
			return walk == rhs.walk && i >= rhs.i;
		}
	};
	static_assert(std::random_access_iterator<ConstEdgeIterator>);

private:
	Mesh const* mesh;

public:
	/**
	 * @brief The vertices of the walk.
	 */
	std::vector<quiver::vertex_index_t> vertices;

	/**
	 * @brief Initializes the walk with just the mesh. The walk will be empty.
	 * The pointer to the `Mesh` is necessary for `ConstEdgeIterator::operator*`.
	 * @param mesh A non-zero pointer to the mesh on which the walk exists.
	 */
	[[nodiscard]] Walk(Mesh const* mesh) noexcept
		: mesh(mesh)
	{
	}
	/**
	 * @brief Initializes the walk with the mesh and a vector of vertices.
	 * The pointer to the `Mesh` is necessary for `ConstEdgeIterator::operator*`.
	 * @param mesh A non-zero pointer to the mesh on which the walk exists.
	 * @param vertices A vector of vertices constituting the walk.
	 */
	[[nodiscard]] Walk(Mesh const* mesh, std::vector<quiver::vertex_index_t> vertices) noexcept
		: mesh(mesh), vertices(std::move(vertices))
	{
	}

	/**
	 * @brief Returns the `Mesh` with which the walk was constructed.
	 * @return A pointer to the `Mesh` with which the walk was constructed.
	 */
	[[nodiscard]] constexpr Mesh const* getMesh() const noexcept
	{
		return mesh;
	}

	/**
	 * @brief Returns an iterator pointing to the first edge of the walk.
	 * @see ConstEdgeIterator
	 */
	[[nodiscard]] ConstEdgeIterator begin() const noexcept
	{
		return { this, 0 };
	}
	[[nodiscard]] std::reverse_iterator<ConstEdgeIterator> rbegin() const noexcept
	{
		return std::reverse_iterator<ConstEdgeIterator>(end());
	}
	/**
	 * @brief Returns an iterator pointing behind the last edge of the walk.
	 * @see ConstEdgeIterator
	 */
	[[nodiscard]] ConstEdgeIterator end() const noexcept
	{
		return { this, std::ptrdiff_t(edges()) };
	}
	[[nodiscard]] std::reverse_iterator<ConstEdgeIterator> rend() const noexcept
	{
		return std::reverse_iterator<ConstEdgeIterator>(begin());
	}

	/**
	 * @brief Returns whether the walk is empty (contains no vertices).
	 * @return `true` iff the walk is empty (contains no vertices).
	 */
	[[nodiscard]] bool isEmpty() const noexcept
	{
		return vertices.empty();
	}
	/**
	 * @brief Returns whether the walk is closed (i.e. the last vertex is equal to the first vertex).
	 * Same as `isLoop`.
	 * @return `true` iff the walk is closed (i.e. the last vertex is equal to the first vertex).
	 * @see isLoop
	 */
	[[nodiscard]] bool isClosed() const noexcept
	{
		return vertices.empty() || vertices.front() == vertices.back();
	}
	/**
	 * @brief Returns whether the walk is open (i.e. the last vertex is not equal to the first vertex).
	 * @return `true` iff the walk is open (i.e. the last vertex is not equal to the first vertex).
	 */
	[[nodiscard]] bool isOpen() const noexcept
	{
		return !isClosed();
	}
	/**
	 * @brief Same as `isClosed`.
	 * @return `true` iff the walk is closed (i.e. the last vertex is equal to the first vertex).
	 * @see isClosed
	 */
	[[nodiscard]] bool isLoop() const noexcept
	{
		return isClosed();
	}
	/**
	 * @brief Returns how many real vertices are contained in the walk (i.e. if the walk is a loop, the double inclusion of the first (and last) vertex is counted as only one).
	 */
	[[nodiscard]] std::size_t realVertices() const noexcept
	{
		if(vertices.empty())
			return 0;
		else
			return vertices.size() - isLoop();
	}
	/**
	 * @brief Returns how many edges are contained in the walk.
	 * @return The number of edges that are contained in the walk.
	 */
	[[nodiscard]] std::size_t edges() const noexcept
	{
		if(vertices.empty())
			return 0;
		else
			return vertices.size() - 1;
	}

	/**
	 * @brief Rotates the vertices in the vertex list of the walk to the left.
	 * For closed walks (loops), the end vertices are handled correctly.
	 * The first (== the last) vertex will not be doubly contained in the rotated walk but the closedness (loopness) of the walk is preserved.
	 * To rotate to the right by `k` vertices, rotate to the left by `realVertices() - k` vertices.
	 * @param shift By how many vertices the vertex list is to be rotated to the left.
	 */
	// left-rotate
	void rotate(std::size_t shift) noexcept
	{
		if(isEmpty())
			return;

		const bool loop = isLoop();
		const std::size_t v = realVertices();
		const std::size_t s = shift % v;

		if(s != 0)
		{
			const auto begin = vertices.begin();
			std::rotate(begin, begin + s, begin + v);
			if(loop)
				vertices.back() = vertices.front();
		}
	}

	friend bool operator==(Walk const&, Walk const&) noexcept = default;
};

/**
 * @brief Serializes a `Walk` object to a stream in plain text.
 */
std::ostream& operator<<(std::ostream& stream, Walk const& walk)
{
	if(walk.isEmpty())
	{
		stream << "empty walk";
	}
	else
	{
		const auto vertexPrinter = [&](quiver::vertex_index_t v)
		{
			stream << "(v:" << v;
			auto const& vertex = walk.getMesh()->vertexGraph.V[v];
			if(vertex.isAffiliated())
			{
				stream << " " << static_cast<PortInformation const&>(vertex);
			}
			stream << ")";
		};
		const auto edgePrinter = [&](quiver::vertex_index_t from, quiver::vertex_index_t to)
		{
			auto const& outEdgePtr = walk.getMesh()->vertexGraph.E(from, to);
			assert(outEdgePtr != nullptr);
			auto const& outEdge = *outEdgePtr;
			stream << "--{e:" << outEdge.edgeIndex;
			if(walk.getMesh()->isExcludedEdge(outEdge.edgeIndex))
			{
				const auto [count, portInformation] = walk.getMesh()->affiliatedPortOfEdge(outEdge.edgeIndex);
				stream << " excl. by " << portInformation;
			}
			stream << "}->";
		};
		quiver::vertex_index_t last = walk.vertices[0];
		vertexPrinter(last);
		for(std::size_t v = 1; v < walk.vertices.size(); ++v)
		{
			const quiver::vertex_index_t current = walk.vertices[v];
			edgePrinter(last, current);
			vertexPrinter(current);
			last = current;
		}
	}
	return stream;
}

/**
 * @brief An edge vector is a function that assigns a value to every edge of a graph.
 * Due to the circumstances in which the type is used, `EdgeVector` is only `int`-valued and will typically only assume the values `-1`, `0` or `1`.
 * @see IsConvertibleToEdgeVector
 * @see ToEdgeVector
 */
class EdgeVector
{
public:
	/**
	 * @brief The data type that we are interested in storing for each edge.
	 */
	using Data = int;
	/**
	 * @brief The underlying data type for storing all the values..
	 */
	using Underlying = Eigen::Matrix<Data, Eigen::Dynamic, 1>;

private:
	Underlying data;

public:
	/**
	 * @brief Constructs an edge vector for `edges` edges.
	 * @param edges The number of edges for which a value should be stored.
	 */
	[[nodiscard]] EdgeVector(std::size_t edges)
		: data(decltype(data)::Zero(edges, 1))
	{
	}
	/**
	 * @brief Constructs an edge vector for as many edges as there are edges in `mesh`.
	 * Equivalent to `EdgeVector(mesh.edges())`.
	 * @param mesh The mesh for which an edge vector should be constructed.
	 */
	[[nodiscard]] EdgeVector(Mesh const& mesh)
		: EdgeVector(mesh.edges())
	{
	}
	/**
	 * @brief Constructs an edge vector for as many edges as there are edges in `mesh` and populates that edge vector according to the given vertex index range `vertices`.
	 * The orientation of the edges described by this vertex index range is taken into account (by using `EdgeOrientation::factor` whose data is stored in `Mesh::vertexGraph`).
	 * @throws InvalidEdgeError iff the vertex index range `vertices` describes an edge not present in `mesh`.
	 * @param mesh The mesh for which an edge vector should be constructed.
	 * @param vertices A range of vertex indices.
	 */
	[[nodiscard]] EdgeVector(Mesh const& mesh, InputRangeOf<std::is_convertible, quiver::vertex_index_t> auto const& vertices)
		: EdgeVector(mesh)
	{
		auto first = std::ranges::begin(vertices);
		auto last = std::ranges::end(vertices);

		if(first == last)
			return;

		quiver::vertex_index_t previous = *first++;
		quiver::vertex_index_t current;
		while(first != last)
		{
			current = *first++;
			const auto outEdge = mesh.vertexGraph.E.get(previous, current);
			if(outEdge == nullptr)
				throw InvalidEdgeError("EdgeVector::EdgeVector: walk contains a nonexistent edge (" + std::to_string(previous) + "--" + std::to_string(current) + ")");

			const EdgeIndex::Index edgeIndex = outEdge->edgeIndex;
			const int orientation = outEdge->factor(previous, current);
			data[edgeIndex] += orientation;

			previous = current;
		}
	}
	/**
	 * @brief Constructs an edge vector from a mesh and a walk.
	 * Equivalent to `EdgeVector(mesh, walk.vertices)`.
	 * @throws InvalidEdgeError iff the vertex sequence in `walk` describes an edge not present in `mesh`.
	 * @param mesh The mesh for which an edge vector should be constructed.
	 * @param walk A walk on this mesh.
	 */
	[[nodiscard]] EdgeVector(Mesh const& mesh, Walk const& walk)
		: EdgeVector(mesh, walk.vertices)
	{
	}

	/**
	 * @brief Returns the edge vector's size.
	 * @return The edge vector's size.
	 */
	[[nodiscard]] std::size_t size() const noexcept
	{
		return data.size();
	}

	/**
	 * @brief Returns a reference to an element of the edge vector.
	 * @param i The index of the element.
	 * @return A reference to the element at the `i`th position.
	 */
	[[nodiscard]] Data& operator[](std::size_t i) noexcept
	{
		return data[i];
	}
	/**
	 * @brief Returns a reference to an element of the edge vector.
	 * @param i The index of the element.
	 * @return A reference to the element at the `i`th position.
	 */
	[[nodiscard]] Data const& operator[](std::size_t i) const noexcept
	{
		return data[i];
	}

	/**
	 * @brief Returns the underlying Eigen matrix.
	 * @return The underlying Eigen matrix.
	 */
	[[nodiscard]] Underlying const& vector() const noexcept
	{
		return data;
	}

	friend bool operator==(EdgeVector const&, EdgeVector const&) noexcept = default;

	EdgeVector& operator+=(EdgeVector const& rhs) noexcept
	{
		data += rhs.data;
		return *this;
	}
	EdgeVector& operator-=(EdgeVector const& rhs) noexcept
	{
		data -= rhs.data;
		return *this;
	}

	friend EdgeVector operator+(EdgeVector lhs, EdgeVector const& rhs) noexcept
	{
		return lhs += rhs;
	}
	friend EdgeVector operator-(EdgeVector lhs, EdgeVector const& rhs) noexcept
	{
		return lhs -= rhs;
	}
};

/**
 * @brief A trait that indicates that a type is convertible to either `EdgeVector` or `Walk`.
 * Useful in conjunction with `ToEdgeVector`.
 * @see ToEdgeVector
 */
template<typename T>
struct IsConvertibleToEdgeVector : public std::bool_constant<std::is_convertible_v<T, EdgeVector> || std::is_convertible_v<T, Walk>>
{
};

/**
 * @brief An overloaded function that guarantees to produce an `EdgeVector`, irrespective of whether it is given a `Walk` or an `EdgeVector`.
 * This is the `EdgeVector` to `EdgeVector` overload which does nothing except for forwarding the given reference.
 * @see IsConvertibleToEdgeVector
 */
[[nodiscard]] constexpr EdgeVector const& ToEdgeVector(Mesh const& mesh, EdgeVector const& edgeVector) noexcept
{
	(void)mesh; // -Wunused-parameter
	return edgeVector;
}
/**
 * @brief An overloaded function that guarantees to produce an `EdgeVector`, irrespective of whether it is given a `Walk` or an `EdgeVector`.
 * This is the `Walk` to `EdgeVector` overload which invokes `EdgeVector`'s conversion constructor.
 * @see IsConvertibleToEdgeVector
 */
[[nodiscard]] EdgeVector ToEdgeVector(Mesh const& mesh, Walk const& walk)
{
	return EdgeVector(mesh, walk);
}

/**
 * @brief A product type that stores a correspondence between indices in two different index systems (a local and a global system).
 * Refer to the thesis for more information.
 * @see OutlineGraphMap
 */
struct LocalGlobalPair
{
	quiver::vertex_index_t local, global;

	friend constexpr bool operator==(LocalGlobalPair, LocalGlobalPair) noexcept = default;
};

/**
 * @brief An annotation type used to inject the member `.globalIndex` into a graph's vertex class.
 * @see LocalGlobalPair
 * @see OutlineGraphMap::OutlineGraph
 * @see OutlineGraphMap::outlineGraph
 */
struct GlobalIndex
{
	quiver::vertex_index_t globalIndex;

	[[nodiscard]] constexpr GlobalIndex(quiver::vertex_index_t globalIndex) noexcept
		: globalIndex(globalIndex)
	{
	}
};

/**
 * @brief A type that stores a vertex mapping (between a local and a global indexing system) and a graph.
 * Used to construct outline graphs as described in the thesis.
 */
struct OutlineGraphMap
{
	using VertexMapping = std::unordered_multimap<quiver::vertex_index_t, quiver::vertex_index_t>; // maps global vertex indices (from mesh) to local vertex indices
	using OutlineGraph = quiver::adjacency_list<quiver::undirected, void, GlobalIndex>;

	/**
	 * @brief An unordered multimap that maps from global to local indices.
	 */
	VertexMapping vertexMapping;
	/**
	 * @brief An undirected graph (the local theater) whose vertices are annotated with their respective global index.
	 * @see GlobalIndex
	 */
	OutlineGraph outlineGraph;

	/**
	 * @brief Returns the number of vertices contained in the outline graph.
	 * @return The number of vertices contained in the outline graph.
	 */
	[[nodiscard]] std::size_t vertices() const noexcept
	{
		return outlineGraph.V.size();
	}

	/**
	 * @brief Inserts a new vertex into the local outline graph.
	 * This new vertex has to belong at least one global vertex as specified by `globalIndex`.
	 * While insertion always takes place and while the new vertex is always locally unique (meaning that it is distinct within `outlineGraph`),
	 * it may correspond to a global vertex for which there already are local correspondences in `outlineGraph`.
	 * Refer to the thesis for a more detailed explanation.
	 * @param globalIndex The index of the global vertex to which the newly inserted local vertex ought to belong.
	 * @return The index pair of the newly inserted vertex.
	 */
	LocalGlobalPair insertVertex(quiver::vertex_index_t globalIndex)
	{
		const quiver::vertex_index_t localIndex = outlineGraph.V.emplace(globalIndex);
		vertexMapping.emplace(globalIndex, localIndex);
		return { localIndex, globalIndex };
	}

	/**
	 * @brief For a given index of a global vertex, this function returns the indices of all local vertices corresponding to that global vertex.
	 * @param globalIndex The index of the global vertex in question.
	 * @return A view (`std::ranges`) of the indices of all local vertices corresponding to the global vertex with index `globalIndex`.
	 */
	[[nodiscard]] auto getLocalIndices(quiver::vertex_index_t globalIndex) const noexcept
	{
		const auto [from, to] = vertexMapping.equal_range(globalIndex);
		return std::ranges::subrange(from, to) | std::ranges::views::transform([](VertexMapping::value_type const& pair) constexpr noexcept
		{
			return pair.second;
		});
	}
	/**
	 * @brief For a given index of a local vertex, this function returns the index of the unique global vertex corresponding to that local vertex.
	 * @param localIndex The index of the local vertex in question.
	 * @return The index of the unique global vertex corresponding to the local vertex with index `localIndex`.
	 */
	[[nodiscard]] quiver::vertex_index_t getGlobalIndex(quiver::vertex_index_t localIndex) const noexcept
	{
		assert(localIndex < outlineGraph.V.size());
		return outlineGraph.V[localIndex].globalIndex;
	}
};

/**
 * @brief A type that stores (indices to) a vertex and its predecessor (as used in graph search algorithms).
 */
struct VertexPredecessorPair
{
	quiver::vertex_index_t vertex;
	quiver::vertex_index_t predecessor;
};

/**
 * @brief Tries to find a walk (an arc) of non-excluded edges between `startVertex` and `endVertex`.
 * If the return value is not `std::nullopt`, a walk can be extracted from such a predecessor vector by simply backtracking,
 * starting from `predecessors[endVertex]`.
 * @param mesh The mesh on which the search takes place.
 * @param startVertex The start point of the search.
 * @param endVertex The end point of the search.
 * @return A vector of predecessors if a walk was found, `std::nullopt` otherwise.
 * @see Mesh::isExcludedEdge
 */
std::optional<std::vector<quiver::vertex_index_t>> getNonexcludedArc(Mesh const& mesh, quiver::vertex_index_t startVertex, quiver::vertex_index_t endVertex)
{
	// Adapted from quiver::bfs
	// TODO: A bidirectional search would be faster.

	constexpr quiver::vertex_index_t noPredecessor = ~quiver::vertex_index_t{};

	std::vector<bool> enqueued(mesh.vertexGraph.V.size(), false);
	std::vector<quiver::vertex_index_t> predecessors(mesh.vertexGraph.V.size(), noPredecessor);
	std::queue<VertexPredecessorPair> neighbors;

	enqueued[startVertex] = true;
	neighbors.emplace(startVertex, noPredecessor);

	do
	{
		const VertexPredecessorPair front = neighbors.front();
		assert(enqueued[front.vertex]);
		auto const& vertex = mesh.vertexGraph.V[front.vertex];

		predecessors[front.vertex] = front.predecessor;
		if(front.vertex == endVertex)
			return std::move(predecessors);

		neighbors.pop();
		for(auto const& edge : vertex.out_edges)
			if(!mesh.isExcludedEdge(edge.edgeIndex) && !enqueued[edge.to])
			{
				neighbors.emplace(edge.to, front.vertex);
				enqueued[edge.to] = true;
			}
	}
	while(!neighbors.empty());
	return std::nullopt;
}

/**
 * @brief For a given vertex on a mesh affiliated to a port, this function computes an `OutlineGraphMap` for the entire port.
 * @throws NoNonExcludedArcError iff somewhere along the port's boundary, no arc (detour of non-excluded edges) could be found.
 * @param mesh The mesh on which the vertex and the port reside.
 * @param startVertex An affiliated vertex.
 * @return A pair containing the `PortInformation` of the port as well as a `OutlineGraphMap` along the boundary of the port.
 * @see PortInformation::isAffiliated
 */
std::pair<PortInformation, OutlineGraphMap> constructOutlineGraphMap(Mesh const& mesh, quiver::vertex_index_t startVertex)
{
	assert(startVertex < mesh.vertexGraph.V.size());

	const PortInformation pi = mesh.vertexGraph.V[startVertex];
	assert(pi.isAffiliated());
	assert(mesh.isPortBoundaryVertex(startVertex));

	OutlineGraphMap outlineGraphMap;

	// In this section, global vertices and edges are those belonging to mesh.
	// On the other hand, local vertices and edges are those belonging to outlineGraphMap.
	// In general they coincide but in some select circumstances, a single global vertex can give rise to
	// multiple local vertices. A local vertex can, however, only belong to a single global vertex.
	// outlineGraphMap maintains the necessary data structures to translate between local and global indices.

	constexpr LocalGlobalPair noneVertex{ ~quiver::vertex_index_t{}, ~quiver::vertex_index_t{} };
	LocalGlobalPair prevVertex = noneVertex;
	LocalGlobalPair currVertex = outlineGraphMap.insertVertex(startVertex);
	const LocalGlobalPair endVertex = currVertex;
	for(;;)
	{
		// In this loop, we traverse along the port's boundary and create the corresponding vertices in outlineGraphMap.
		// The exact procedure by which we connect these boundary vertices by edges depends on whether we're dealing with
		// a magnetic (simple) or an electric (more complicated) port.

		LocalGlobalPair nextVertex = noneVertex;
		for(auto const& outEdge : mesh.vertexGraph.V[currVertex.global].out_edges)
		{
			// We search for an edge that is both not our previous edge as well as a port boundary edge.
			if(outEdge.to != prevVertex.global && mesh.isPortBoundaryEdge(outEdge.edgeIndex))
			{
				nextVertex.global = outEdge.to;
				break;
			}
		}

		// Check that we have found something.
		assert(nextVertex != noneVertex);

		// Insert a new vertex into outlineGraphMap except if we have reached the starting point again.
		const bool finished = (nextVertex.global == endVertex.global);
		if(finished)
			nextVertex = endVertex;
		else
			nextVertex = outlineGraphMap.insertVertex(nextVertex.global);

		if(pi.portType == magneticPort)
		{
			// While, for magnetic ports, we can simply insert the edge on the boundary of the port...
			outlineGraphMap.outlineGraph.E.emplace(currVertex.local, nextVertex.local);
		}
		else // if(pi.portType == electricPort)
		{
			// ...this is not possible for electric ports because the edges on electric ports' boundaries are excluded.
			// Instead, we find the the shortest "replacement" for this impermissible boundary edge.
			// These replacements are referred to as "arcs" because they look like little arcs coming out of and going back into ports.
			// We perform a BFS search from the beginning vertex to the end vertex on the global mesh but with all
			// excluded edges removed. We don't actually create such a mesh, we simply ignore excluded edges during our search.

			// Global predecessors.
			const auto arcPredecessors = getNonexcludedArc(mesh, currVertex.global, nextVertex.global);

			// If this check fails, it means that the arc over the boundary could not be found successfully.
			// Generally, this means that one of the preconditions has not been met, e.g. the precondition that the mesh shall
			// have no boundaries. However, it may also fail for well-formed meshes if two electic ports are too close to each other.
			// Concretely, two electric ports being "too close" means that there is a portion where no walk of non-excluded edges
			// can pass between them.
			// The solution is to refine the mesh until the ports in question are not too close to each other anymore.
			// This cannot be done unilaterally within in the cohomology subproject because, if we were to alter the mesh,
			// there would be a discrepancy between the mesh and the outside data structures that were used to construct the
			// mesh (i.e. the lf::mesh::Mesh instance with which projects::cohomology::Mesh::Mesh was invoked).
			// The user should refine the outside data structures (usually the lf::mesh::Mesh instance) and then rerun the
			// entire cohomology pipeline.
			if(!arcPredecessors)
				throw NoNonExcludedArcError("constructOutlineGraphMap: failed to find non-excluded arc between vertices #"
					+ std::to_string(currVertex.global) + " and #" + std::to_string(nextVertex.global)
					+ " while constructing the boundary of " + to_string(pi));

			// Backtrack along the arc and insert missing vertices and edges into outlineGraphMap.
			// This is the place where a single global vertex can give rise to multiple local vertices.
			// This is an indisposable mechanism because if we were to naïvely copy the local structure of
			// the global graph into the boundary graph, the boundary could become homologically nontrivial,
			// for instance in the case that a mesh vertex is in the outline of the same port twice.
			for(LocalGlobalPair arcVertex = nextVertex;;)
			{
				const quiver::vertex_index_t globalPredecessor = (*arcPredecessors)[arcVertex.global];
				const bool arrived = (globalPredecessor == currVertex.global);
				const LocalGlobalPair predecessor = arrived ? currVertex : outlineGraphMap.insertVertex(globalPredecessor);

				outlineGraphMap.outlineGraph.E.emplace(predecessor.local, arcVertex.local);

				// Iteration condition.
				if(arrived)
					break;

				// Iteration step: advance to the next vertex in the arc.
				arcVertex = predecessor;
			}
		}

		// Iteration condition.
		if(finished)
			break;

		// Iteration step: advance to the next vertex on the boundary.
		prevVertex = currVertex;
		currVertex = nextVertex;
	}

	// Optional: Verify that the generated boundary graph is indeed homologically equivalent to a single loop.
	//   We first verify that the graph has only one connected components.
	//   Then we verify that the graph's cyclomatic number is 1. This is done by checking whether
	//   r = m - n + c holds where r is the cyclomatic number (a.k.a. the circuit rank) and where c is
	//   the number of connected components which we already verified to be 1. The equality thus simplifies to m = n.
	assert(quiver::ccs(outlineGraphMap.outlineGraph) == 1);
	assert(outlineGraphMap.outlineGraph.E.size() == outlineGraphMap.outlineGraph.V.size());

	return { pi, std::move(outlineGraphMap) };
};

/**
 * @brief For a walk on a mesh, reroute the walk so as to not contain excluded edges anymore.
 * Refer to the thesis for more details.
 * @throws std::runtime_error iff rerouting is not possible.
 * @param mesh The mesh on which the walk resides.
 * @param walk A walk on that mesh.
 * @see Mesh::isExcludedEdge
 */
void pushLoopsOutOfPorts(Mesh const& mesh, Walk& walk)
{
	// Could theoretically also be extended to non-loops, but not necessary within the scope of this work.
	assert(walk.isLoop());

	// Step 1: Check if there is at least one excluded edge in the generated walk.
	constexpr std::size_t notFound = ~0;
	std::size_t firstVertexOfFirstExcludedEdge = notFound;
	for(std::size_t i = 0; auto&& [from, to, outEdge] : walk)
	{
		const EdgeIndex::Index edgeIndex = outEdge.edgeIndex;
		if(mesh.isExcludedEdge(edgeIndex))
		{
			firstVertexOfFirstExcludedEdge = i;
			break;
		}
		++i;
	}

	if(firstVertexOfFirstExcludedEdge != notFound)
	{
		// Step 2: Rotate the walk until the wraparound isn't an excluded edge anymore.
		//         This is to ensure that we can detour all the port crossings "left to right", including the first one.
		const auto getWraparound = [](Walk const& walk){ return std::get<2>(*walk.rbegin()); };
		if(mesh.isExcludedEdge(getWraparound(walk).edgeIndex))
		{
			std::size_t rotation = walk.realVertices();
			for(auto iter = walk.rbegin(); iter != walk.rend(); ++iter, --rotation)
				if(!mesh.isExcludedEdge(std::get<2>(*iter).edgeIndex))
					break;
			assert(rotation <= walk.realVertices()); // no underflow
			walk.rotate(rotation);

			firstVertexOfFirstExcludedEdge = 0;
			assert(mesh.isExcludedEdge(std::get<2>(*walk.begin()).edgeIndex));
		}
		assert(!mesh.isExcludedEdge(getWraparound(walk).edgeIndex));

		// Step 3: Create a table of boundary graphs for all (necessary) ports.
		std::unordered_map<PortInformation, OutlineGraphMap> outlineGraphMaps;
		const auto populate = [&](quiver::vertex_index_t startVertex)
		{
			auto [pi, outlineGraphMap] = constructOutlineGraphMap(mesh, startVertex);
			[[maybe_unused]] const auto [iter, inserted] = outlineGraphMaps.emplace(pi, std::move(outlineGraphMap));
			assert(inserted);
			return iter;
		};

		// Step 4: Find detours around all intersected ports and weave them into the topological cycle.
		const auto getDetour = [&](quiver::vertex_index_t from, quiver::vertex_index_t to) -> Walk
		{
			const PortInformation portInformation = mesh.vertexGraph.V[from];
			auto finding = outlineGraphMaps.find(portInformation);
			if(finding == outlineGraphMaps.end())
				finding = populate(from);
			OutlineGraphMap const& outlineGraphMap = finding->second;

			std::size_t length;
			quiver::vertex_index_t localEnd;
			constexpr quiver::vertex_index_t noPredecessor = ~quiver::vertex_index_t{};
			std::vector<quiver::vertex_index_t> predecessors(outlineGraphMap.vertices(), noPredecessor);
			const auto bfsVisitor = [&](quiver::vertex_index_t vertexIndex, quiver::bfs_path_length_t distance, quiver::vertex_index_t predecessor)
			{
				predecessors[vertexIndex] = predecessor;
				const bool reachedEnd = (outlineGraphMap.getGlobalIndex(vertexIndex) == to);
				if(reachedEnd)
				{
					length = distance;
					localEnd = vertexIndex;
				}
				return reachedEnd;
			};
			const bool found = quiver::bfs_shortest_path(outlineGraphMap.outlineGraph, outlineGraphMap.getLocalIndices(from), bfsVisitor);
			if(!found)
				throw std::runtime_error("pushLoopsOutOfPorts: cannot resolve intersection between walk and " + to_string(portInformation));

			Walk detour(&mesh);
			detour.vertices.reserve(length + 1);
			for(quiver::vertex_index_t v = localEnd; outlineGraphMap.getGlobalIndex(v) != from; v = predecessors[v])
				detour.vertices.emplace_back(outlineGraphMap.getGlobalIndex(v));
			detour.vertices.emplace_back(from);
			assert(detour.vertices.size() == length + 1);

			return detour;
		};

		{
			quiver::vertex_index_t portStart = firstVertexOfFirstExcludedEdge;
			bool inPort = true;
			assert(static_cast<PortInformation>(mesh.vertexGraph.V[walk.vertices[firstVertexOfFirstExcludedEdge]]) == static_cast<PortInformation>(mesh.vertexGraph.V[walk.vertices[firstVertexOfFirstExcludedEdge + 1]]));
			PortInformation portInformation = mesh.vertexGraph.V[walk.vertices[firstVertexOfFirstExcludedEdge]];
			for(std::size_t i = firstVertexOfFirstExcludedEdge + 1; i < walk.vertices.size(); ++i)
			{
				const quiver::vertex_index_t previous = walk.vertices[i - 1];
				const quiver::vertex_index_t current = walk.vertices[i];

				EdgeIndex::Index edgeIndex;
				const auto edge = mesh.vertexGraph.E.get(previous, current);
				assert(edge != nullptr);
				edgeIndex = edge->edgeIndex;

				if(mesh.isExcludedEdge(edgeIndex))
				{
					if(!inPort)
					{
						portStart = i - 1;
						inPort = true;
						assert(static_cast<PortInformation>(mesh.vertexGraph.V[previous]) == static_cast<PortInformation>(mesh.vertexGraph.V[current]));
						portInformation = mesh.vertexGraph.V[current];
					}
					else
					{
						assert(static_cast<PortInformation>(mesh.vertexGraph.V[current]) == portInformation);
					}
				}
				else
				{
					if(inPort)
					{
						const quiver::vertex_index_t portEnd = i - 1;

						// assert(portStart < portEnd);
						Walk detour = getDetour(walk.vertices[portStart], walk.vertices[portEnd]);
						const std::size_t growth = detour.vertices.size() - (portEnd - portStart + 1);
						// assert(detour.vertices.size() >= portEnd - portStart + 1);
						// assert(*detour.vertices.rbegin() == walk.vertices[portStart]);
						// assert(*std::prev(detour.vertices.rend()) == walk.vertices[portEnd]);

						const auto mid = std::prev(detour.vertices.rend(), growth);
						std::copy(std::next(detour.vertices.rbegin()), mid, std::next(walk.vertices.begin(), portStart + 1));
						walk.vertices.insert(walk.vertices.begin() + portEnd + 1, mid, detour.vertices.rend());

						i += growth;
						inPort = false;
						portInformation = {};
					}
				}
			}
			assert(!inPort);
		}
	}
}

/**
 * @brief Computes all topological cycles of a mesh.
 * Refer to the thesis for more details.
 * @param mesh The mesh.
 * @see Walk
 */
[[nodiscard]] std::vector<Walk> getTopologicalCycles(Mesh const& mesh)
{
	// Step 1: Compute a MST of the face graph.
	auto faceGraphMST = quiver::kruskal(mesh.faceGraph);

	// Step 2: Remove all edges in the vertex graph that the MST of the face graph crosses.
	std::unordered_set<EdgeIndex::Index> removedEdges;
	for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : faceGraphMST.V)
	{
		for(auto const& outEdge : vertex.out_edges)
			if(vertexIndex <= outEdge.to)
				removedEdges.insert(outEdge.edgeIndex);
		++vertexIndex;
	}

	auto culledVertexGraph = mesh.vertexGraph.strip_edges();
	for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : mesh.vertexGraph.V)
	{
		for(auto const& outEdge : vertex.out_edges)
			if(vertexIndex < outEdge.to)
				if(removedEdges.count(outEdge.edgeIndex) == 0)
					culledVertexGraph.E.emplace(vertexIndex, outEdge.to, outEdge.properties());
		++vertexIndex;
	}

	// Step 3: Compute a MST of the reduced vertex graph.
	// TODO: Due to the small amount of edges to be removed here, maybe a top-down algorithm like Reverse-delete is superior to bottom-up algorithms like Kruskal's.
	auto culledVertexGraphMST = quiver::kruskal(culledVertexGraph);

	// Step 4: Find the "buckles": The only edges removed in the previous step.
	std::vector<Walk> topologicalCycles;
	topologicalCycles.reserve(culledVertexGraph.E.size() - culledVertexGraphMST.E.size());
	for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : culledVertexGraph.V)
	{
		for(auto const& outEdge : vertex.out_edges)
			if(vertexIndex <= outEdge.to)
				if(!culledVertexGraphMST.E(vertexIndex, outEdge.to))
					topologicalCycles.push_back(Walk(&mesh, {vertexIndex, outEdge.to}));
		++vertexIndex;
	}

	// Step 5: For each buckle, find the associated topological cycle by (temporarily) removing the buckle and performing a search from the buckle's vertices.
	for(Walk& walk : topologicalCycles)
	{
		assert(walk.vertices.size() == 2);
		const quiver::vertex_index_t buckleStart = walk.vertices[0];
		const quiver::vertex_index_t buckleEnd = walk.vertices[1];
		constexpr quiver::vertex_index_t noPredecessor = ~quiver::vertex_index_t{};
		std::vector<quiver::vertex_index_t> predecessors(culledVertexGraph.V.size(), noPredecessor);
		std::size_t length;
		const auto bfsVisitor = [&](quiver::vertex_index_t v, quiver::bfs_path_length_t distance, quiver::vertex_index_t predecessor)
		{
			predecessors[v] = predecessor;
			if(v == buckleEnd)
			{
				length = distance;
				return true;
			}
			else
				return false;
		};

		auto properties = culledVertexGraph.E.get(buckleStart, buckleEnd)->properties();
		culledVertexGraph.E.erase(buckleStart, buckleEnd);
		[[maybe_unused]] const bool found = quiver::bfs_shortest_path(culledVertexGraph, buckleStart, bfsVisitor);
		assert(found);
		culledVertexGraph.E.emplace(buckleStart, buckleEnd, std::move(properties));

		walk.vertices.clear();
		walk.vertices.reserve(length + 2);
		for(quiver::vertex_index_t v = buckleEnd; v != buckleStart; v = predecessors[v])
			walk.vertices.emplace_back(v);
		walk.vertices.emplace_back(buckleStart);
		walk.vertices.emplace_back(buckleEnd);
		assert(walk.vertices.size() == length + 2);

		// Step 6: Make sure that the cycles don't intersect any ports.
		pushLoopsOutOfPorts(mesh, walk);
	}

	return topologicalCycles;
}

/**
 * @brief Computes all electric connectors of a mesh.
 * Refer to the thesis for more details.
 * @param mesh The mesh.
 * @see Walk
 */
[[nodiscard]] std::vector<Walk> getElectricConnectors(Mesh const& mesh)
{
	// Step 1: Construct an undirected vertex graph containing no interior edges.
	quiver::adjacency_list<quiver::undirected, EdgeIndex, void> exteriorVertexGraph(mesh.vertexGraph.V.size());
	for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : mesh.vertexGraph.V)
	{
		for(auto const& outEdge : vertex.out_edges)
			if(vertexIndex < outEdge.to)
				if(!mesh.isExcludedEdge(outEdge.edgeIndex))
					exteriorVertexGraph.E.emplace(vertexIndex, outEdge.to, outEdge.edgeIndex);
		++vertexIndex;
	}

	// Step 2: Compute the set of boundary vertices for each electric port.
	// port id -> {boundary vertices}
	std::unordered_map<PortInformation::PortIdType, std::unordered_set<quiver::vertex_index_t>> boundaryVertices;
	for(auto const& face : mesh.faceGraph.V)
		if(face.portType == electricPort)
			for(auto const& outEdge : face.out_edges)
				if(mesh.isElectricalBoundaryEdge(outEdge.edgeIndex))
					for(quiver::vertex_index_t vertex : mesh.edgeConnectivity[outEdge.edgeIndex].vertices)
						boundaryVertices[face.id].emplace(vertex);

	// N_E
	const std::size_t electricPortCount = boundaryVertices.size();

	// Step 3: Compute the shortest distances between the ports.
	// Some types that help us with the BFS.
	constexpr quiver::vertex_index_t noPredecessor = ~quiver::vertex_index_t{};
	struct Predecessors
	{
		std::vector<quiver::vertex_index_t> predecessors;
	};
	struct Endpoint
	{
		quiver::vertex_index_t endpoint;
	};

	// Unfortunately, we have to employ this awkward index translation table here because the port ids are not guaranteed to be a permutation of {0, ..., N_E - 1}.
	// port id -> {0, ..., N_E - 1}
	std::unordered_map<PortInformation::PortIdType, std::size_t> portId2portIndex;
	{
		std::size_t portIndex = 0;
		for(auto const& [portId, vertices] : boundaryVertices)
		{
			portId2portIndex[portId] = portIndex;
			++portIndex;
		}
	}

	// V = ports, E = walks between ports
	quiver::adjacency_list<quiver::undirected, quiver::cmb<quiver::wt<quiver::bfs_path_length_t>, Endpoint>, quiver::cmb<PortInformation, Predecessors>> portReachability(electricPortCount);
	{
		quiver::vertex_index_t portIndex = 0;
		for(auto const& [portId, vertices] : boundaryVertices)
		{
			portReachability.V[portIndex].portType = electricPort;
			portReachability.V[portIndex].id = portId;
			portReachability.V[portIndex].predecessors.resize(exteriorVertexGraph.V.size(), noPredecessor);
			++portIndex;
		}
	}
	{
		// We use this bool array to keep track of which ports we need to reach before we can move on to the next loop iteration.
		// At the beginning of every iteration, we set all relevant bits to true and we set them to false as soon as that
		// specific port has been reached.
		// Example:
		//   For port 0, we abort as soon as we have reached ports 1, 2, ..., N_E - 1
		//   For port 1, we abort as soon as we have reached ports 2, 3, ..., N_E - 1
		//     There is no need to reach port 0 again because the shortest path from port 0 to port 1 is equivalent
		//     to the shortest path from port 1 to port 0, which has already been computed in the first iteration.
		std::vector<bool> unreachedPorts(electricPortCount); // deliberately uninitialized

		for(quiver::vertex_index_t vertexIndex = 0; auto& port : portReachability.V)
		{
			// Skip the last port entirely because all shortest distances to it have already been calculated.
			if(vertexIndex + 1 >= electricPortCount)
				break;

			// The lower bits are never accessed again so we don't need to set them to a defined value.
			// The higher bits are all set to true and whenever we reach that port, we set it to false.
			const auto unreachedPortsBegin = unreachedPorts.begin() + (vertexIndex + 1);
			std::fill(unreachedPortsBegin, unreachedPorts.end(), true);

			// We also keep track of the number of ports remaining so that we don't have to check the entire
			// bool array every time.
			int remainingPorts = unreachedPorts.size() - (vertexIndex + 1);

			const auto bfsVisitor = [&](quiver::vertex_index_t v, quiver::bfs_path_length_t distance, quiver::vertex_index_t predecessor)
			{
				port.predecessors[v] = predecessor;
				auto const& vertex = mesh.vertexGraph.V[v];

				// If we have reached a vertex that belongs to an electric port...
				if(vertex.portType == electricPort)
				{
					// ...we check whether we wanted to find it at all and whether it has been found already.
					const unsigned portIndex = portId2portIndex[vertex.id];
					if(portIndex > vertexIndex && unreachedPorts[portIndex])
					{
						// If we wanted to find it but it hasn't been found before, we discovered a new
						// shortest path between two ports and we record that accordingly.
						portReachability.E.emplace(vertexIndex, portIndex, distance, v);
						unreachedPorts[portIndex] = false;
						--remainingPorts;
						// As soon as we have found all ports that we wanted to find, we return true which lets the
						// BFS algorithm know that we found what we were looking for which makes it terminate early.
						if(remainingPorts == 0)
							return true;
					}
				}

				return false;
			};
			[[maybe_unused]] const bool found = quiver::bfs_shortest_path(exteriorVertexGraph, boundaryVertices.find(port.id)->second, bfsVisitor);
			// Optional check and exception in case a port was not reachable from this particular port.
			// This check should be disabled if disjoint unions of polyhedra should also be supported (i.e. the input mesh consisting of multiple connected components).
			/*
			if(!found)
			{
				assert(remainingPorts != 0);
				const auto unreachablePort = std::find(unreachedPortsBegin, unreachedPorts.end(), true);
				assert(unreachablePort != unreachedPorts.end());
				const unsigned unreachablePortIndex = std::distance(unreachedPorts.begin(), unreachablePort);
				const PortInformation unreachablePortInformation(electricPort, std::next(boundaryVertices.begin(), unreachablePortIndex)->first);
				throw std::runtime_error("getElectricConnectors: " + to_string(static_cast<PortInformation const&>(port)) + " cannot reach " + to_string(unreachablePortInformation));
			}
			*/

			++vertexIndex;
		}
	}

	// Step 4: Compute the MST of the electric port connectivity graph.
	auto portReachabilityMST = quiver::kruskal(std::move(portReachability));
	// Optional check. If there are less MST edges than MST vertices - 1, it means the graph had multiple connected components.
	// This check should be disabled if disjoint unions of polyhedra should also be supported (i.e. the input mesh consisting of multiple connected components).
	/*
	if(portReachabilityMST.E.size() < portReachabilityMST.V.size() - 1)
		throw std::runtime_error("getElectricConnectors: The port reachability graph has multiple connected components.");
	*/

	// Step 5: Reconstruct the shortest paths between the electric ports.
	std::vector<Walk> electricConnectors(portReachabilityMST.E.size(), Walk(&mesh));
	{
		// The index in the electricConnectors vector.
		std::size_t i = 0;

		for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : portReachabilityMST.V)
		{
			for(auto const& outEdge : vertex.out_edges)
			{
				// As both portReachability as well as portReachabilityMST are undirected graphs, this makes
				// a vertex 'x' contain an out edge to 'y' while 'y' simultaneously contains an out edge to 'x'.
				// To iterate over all edges only once, we use this inequality condition below.
				if(vertexIndex < outEdge.to)
				{
					auto& walk = electricConnectors[i++].vertices;
					const auto isVertexOfStartingPort = [&](quiver::vertex_index_t v)
					{
						return boundaryVertices[vertex.id].contains(v);
					};

					// We step backwards through the predecessor array until we reach the original boundary vertex.
					quiver::vertex_index_t v;
					for(v = outEdge.endpoint; !isVertexOfStartingPort(v); v = vertex.predecessors[v])
						walk.emplace_back(v);
					walk.emplace_back(v);
				}
			}
			++vertexIndex;
		}
	}

	return electricConnectors;
}

/**
 * @brief Computes all magnetic port cycles of a mesh.
 * Refer to the thesis for more details.
 * @param mesh The mesh.
 * @see Walk
 */
[[nodiscard]] std::vector<EdgeVector> getMagneticPortCycles(Mesh const& mesh)
{
	// Just like in getElectricConnectors, we unfortunately have to make use of a hash table here
	// because the port ids are not guaranteed to be a permutation of {0, ..., N_M - 1}.
	std::unordered_set<PortInformation::PortIdType> finishedPortIds;

	// Step 1: Find all magnetic port boundaries.
	std::vector<EdgeVector> magneticPortCycles;
	std::vector<std::size_t> magneticPortCycleEdges;

	// Iterate over all vertex graph edges.
	for(quiver::vertex_index_t vertexIndex = 0; auto const& vertex : mesh.vertexGraph.V)
	{
		for(auto const& outEdge : vertex.out_edges)
		{
			// Ensure uniqueness.
			if(!(vertexIndex < outEdge.to))
				continue;

			// Ensure that the edge is a port boundary edge.
			const auto optionalPortInformation = mesh.isPortBoundaryEdge(outEdge.edgeIndex);
			if(!optionalPortInformation)
				continue;

			// Ensure that the port is a magnetic port.
			PortInformation const& pi = *optionalPortInformation;
			if(pi.portType != PortType::magneticPort)
				continue;

			// Ensure that the respective port has not already been treated.
			[[maybe_unused]] const auto [iter, inserted] = finishedPortIds.emplace(pi.id);
			if(!inserted)
				continue;

			// Prepare for a new edge vector.
			magneticPortCycles.emplace_back(mesh.edges());
			EdgeVector& edgeVector = magneticPortCycles.back();
			magneticPortCycleEdges.emplace_back(0);
			std::size_t& edgeVectorSize = magneticPortCycleEdges.back();

			// Traverse along the port's boundary and update the edge vector accordingly.
			constexpr quiver::vertex_index_t noneVertex{ ~quiver::vertex_index_t{} };
			quiver::vertex_index_t prevVertex = vertexIndex;
			quiver::vertex_index_t currVertex = outEdge.to;
			decltype(Mesh::vertexGraph)::out_edge_t const* edge = &outEdge;
			const quiver::vertex_index_t endVertex = currVertex;
			for(;;)
			{
				// Register the edge {prevVertex, currVertex} in edgeVector and edgeVectorSize while paying attention to the direction.
				const int orientation = edge->factor(prevVertex, currVertex);
				int& edgeValue = edgeVector[edge->edgeIndex];
				edgeVectorSize -= (edgeValue != 0);
				edgeValue += orientation;
				edgeVectorSize += (edgeValue != 0);

				// Search for the next edge and vertex along the port's boundary.
				quiver::vertex_index_t nextVertex = noneVertex;
				for(auto const& outEdge : mesh.vertexGraph.V[currVertex].out_edges)
				{
					// It cannot be our previous vertex.
					if(outEdge.to == prevVertex)
						continue;

					// It has to be a port boundary edge.
					const auto optionalPortInformation = mesh.isPortBoundaryEdge(outEdge.edgeIndex);
					if(!optionalPortInformation)
						continue;

					// Check that we haven't somehow left our port for another port.
					assert(*optionalPortInformation == pi);

					// We found the successor edge.
					nextVertex = outEdge.to;
					edge = &outEdge;
					break;
				}

				// Check that we have found something.
				assert(nextVertex != noneVertex);

				// Iteration condition.
				if(nextVertex == endVertex)
					break;

				// Iteration step: advance to the next vertex on the boundary.
				prevVertex = currVertex;
				currVertex = nextVertex;
			}

			// Check whether the edge count was correctly updated.
			assert(edgeVectorSize == (edgeVector.vector().array() != 0).count());
		}
		++vertexIndex;
	}

	// Return if there's less than two magnetic ports.
	if(magneticPortCycles.size() <= 1)
		return {};

	// This function returns for a boundary path whether the paths agrees (1) or disagrees (-1)
	// with the internal order of (the faces of) the port it bounds.
	const auto getBoundaryOrientation = [&mesh](EdgeVector const& edgeVector) -> int
	{
		// Find the first edge that appears in this path.
		for(quiver::vertex_index_t e = 0; e < edgeVector.size(); ++e)
		{
			const int edgeDirection = edgeVector[e];
			if(edgeDirection != 0)
			{
				// Find the port face that this edge touches.
				const EdgeInformation& edgeInformation = mesh.edgeConnectivity[e];
				const auto& faceIndices = edgeInformation.faces;
				assert(mesh.faceGraph.V[faceIndices[0]].isAffiliated()
					 ^ mesh.faceGraph.V[faceIndices[1]].isAffiliated());
				const auto& portFace = mesh.faceGraph.V[faceIndices[mesh.faceGraph.V[faceIndices[0]].isAffiliated() ? 0 : 1]];

				// Now multiply the two edge direction prescriptions:
				//   The path prescribes one direction while the internal order of the face prescribes another.
				// If they agree about the edge direction, the result is 1, otherwise it is -1.
				// Because the mesh has been oriented in the precomputation step, it doesn't matter
				// for which edge-face pair we compute this. The result is always the same.
				return edgeDirection * portFace.orderOf(edgeInformation.vertices[0], edgeInformation.vertices[1]);
			}
		}
		assert(false);
		return 0;
	};

	// Step 2: Construct the actual cycles from the boundaries by adding/subtracting one magnetic port boundary from all the others.

	// Single out the shortest boundary.
	const auto shortestBoundaryIndex = std::distance(magneticPortCycleEdges.begin(), std::min_element(magneticPortCycleEdges.begin(), magneticPortCycleEdges.end()));
	const EdgeVector shortestBoundary = std::move(magneticPortCycles[shortestBoundaryIndex]);
	magneticPortCycles.erase(magneticPortCycles.begin() + shortestBoundaryIndex);
	const bool shortestBoundaryOrientation = std::signbit(getBoundaryOrientation(shortestBoundary));

	// Add/subtract it to/from all other boundaries.
	for(auto iter = magneticPortCycles.begin(); iter != magneticPortCycles.end(); ++iter)
	{
		// The current and the short boundary cannot both agree or both disagree with their respective ports.
		// One has to disagree while the other has to agree. Hence we check the orientations for inequality.
		// If the orientations are equal, we subtract the short boundary which it flip its orientation.
		if(std::signbit(getBoundaryOrientation(*iter)) != shortestBoundaryOrientation)
			*iter += shortestBoundary;
		else
			*iter -= shortestBoundary;
	}

	return magneticPortCycles;
}

/**
 * @brief A type that holds all cycles of the relative cohomology.
 */
class RelativeCohomology
{
	std::size_t topologicalN, electricN;

public:
	/**
	 * @brief The cycles, represented as `EdgeVector`s.
	 */
	std::vector<EdgeVector> cycles;

	/**
	 * @brief Copies/converts the three given sets of cycles into the data structure.
	 * @param mesh The mesh whose relative cohomology is being described.
	 * @param topologicalCycles The range of topological cycles to copy.
	 * @param electricConnectors The range of electric connectors to copy.
	 * @param magneticPortCycles The range of magnetic port cycles to copy.
	 */
	[[nodiscard]] explicit RelativeCohomology(Mesh const& mesh, InputRangeOf<IsConvertibleToEdgeVector> auto const& topologicalCycles, InputRangeOf<IsConvertibleToEdgeVector> auto const& electricConnectors, InputRangeOf<IsConvertibleToEdgeVector> auto const& magneticPortCycles)
	{
		for(auto&& topologicalCycle : topologicalCycles)
			cycles.emplace_back(ToEdgeVector(mesh, topologicalCycle));
		topologicalN = cycles.size();

		for(auto&& electricConnector : electricConnectors)
			cycles.emplace_back(ToEdgeVector(mesh, electricConnector));
		electricN = cycles.size() - topologicalN;

		for(auto&& magneticPortCycles : magneticPortCycles)
			cycles.emplace_back(ToEdgeVector(mesh, magneticPortCycles));
	}
	/**
	 * @brief Computes a generating set of cycles for the given `Mesh`.
	 * @param mesh The mesh for which to compute the generating set of cycles.
	 */
	[[nodiscard]] explicit RelativeCohomology(Mesh const& mesh)
		: RelativeCohomology(mesh, getTopologicalCycles(mesh), getElectricConnectors(mesh), getMagneticPortCycles(mesh))
	{
	}

	/**
	 * @brief Returns a view to the portion of the cycles generated by `getTopologicalCycles`.
	 * @return A view to the portion of the cycles generated by `getTopologicalCycles`.
	 */
	[[nodiscard]] auto topological() const
	{
		return std::ranges::subrange(cycles.begin(), cycles.begin() + topologicalN);
	}
	/**
	 * @brief Returns a view to the portion of the cycles generated by `getElectricConnectors`.
	 * @return A view to the portion of the cycles generated by `getElectricConnectors`.
	 */
	[[nodiscard]] auto electric() const
	{
		return std::ranges::subrange(cycles.begin() + topologicalN, cycles.begin() + topologicalN + electricN);
	}
	/**
	 * @brief Returns a view to the portion of the cycles generated by `getMagneticPortCycles`.
	 * @return A view to the portion of the cycles generated by `getMagneticPortCycles`.
	 */
	[[nodiscard]] auto magnetic() const
	{
		return std::ranges::subrange(cycles.begin() + topologicalN + electricN, cycles.end());
	}
};

} // namespace projects::cohomology

#endif // !defined(COHOMOLOGY_COHOMOLOGY_H_INCLUDED)

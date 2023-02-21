/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#ifndef COHOMOLOGY_TEST_EXAMPLES_H_INCLUDED
#define COHOMOLOGY_TEST_EXAMPLES_H_INCLUDED

#include <cohomology.h>
#include <utility.h>

#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/base/span.h>
#include <Eigen/Dense>

#include <array>
#include <memory>
#include <functional>
#include <utility>
#include <cmath>

/**
 * @brief Contains a number of mesh generators and example problems to be used for tests and/or illustrations.
 */
namespace projects::cohomology::examples
{

enum class Shape { tris, quads };
char const* to_string(Shape shape)
{
	switch(shape)
	{
	case Shape::tris:
		return "tris";
	case Shape::quads:
		return "quads";
	}
	assert(false);
	return nullptr;
}

using MeshFactory = lf::mesh::hybrid2d::MeshFactory;
using PointIndex = MeshFactory::size_type;
using QuadrilateralIndices = std::array<PointIndex, 4>;

template<typename T, typename Map>
struct CastMap
{
	// L-Ref or value.
	Map map;

	template<typename... Args>
	constexpr CastMap(Args&&... args)
	: map(std::forward<Args>(args)...)
	{
	}

	template<typename... Args>
	constexpr void operator()(Args&&... args)
	{
		return map(T{ std::forward<Args>(args)... });
	}
};

template<typename T, typename Map>
[[nodiscard]] constexpr auto castMap(Map&& map) noexcept
{
    return CastMap<T, Map>(std::forward<Map>(map));
}

[[nodiscard]] MeshFactory getGridlikeMeshFactory(unsigned xMax, unsigned yMax, unsigned zMax)
{
	MeshFactory meshFactory(3);
	Eigen::VectorXd p(3);
	for(unsigned z = 0; z < zMax; ++z)
	{
		for(unsigned y = 0; y < yMax; ++y)
		{
			for(unsigned x = 0; x < xMax; ++x)
			{
				p << x, y, z;
				meshFactory.AddPoint(p);
			}
		}
	}
	return meshFactory;
}

[[nodiscard]] constexpr auto getGridlikePointMapping(unsigned xMax, unsigned yMax, unsigned zMax) noexcept
{
	return [=](unsigned x, unsigned y, unsigned z) constexpr noexcept -> PointIndex
	{
		assert(x < xMax);
		assert(y < yMax);
		assert(z < zMax);
		return (z * yMax + y) * xMax + x;
	};
}

[[nodiscard]] constexpr auto getToroidalPointMapping(unsigned thetaSamples, unsigned phiSamples) noexcept
{
	return [=](unsigned theta, unsigned phi) constexpr noexcept -> PointIndex
	{
		return (theta % thetaSamples) * phiSamples + (phi % phiSamples);
	};
}

[[nodiscard]] constexpr auto getInverseToroidalPointMapping(unsigned phiSamples) noexcept
{
	return [=](unsigned vertex) constexpr noexcept -> std::pair<unsigned, unsigned>
	{
		const unsigned theta = vertex / phiSamples;
		const unsigned phi = vertex % phiSamples;
		return { theta, phi };
	};
}

void addQuadrilateral(MeshFactory& meshFactory, Shape shape, QuadrilateralIndices indices)
{
	if(shape == Shape::tris)
	{
		using Span = nonstd::span<const PointIndex>;
		for(Span span : { Span(indices.data(), 3), Span(indices.data() + 1, 3) })
		{
			Eigen::Matrix<double, Eigen::Dynamic, 3> coords(3, 3);
			for(int i = 0; i < 3; ++i)
				coords.col(i) = getPointPosition(meshFactory.GetPoint(span[i]));
			meshFactory.AddEntity(lf::base::RefEl::kTria(), span, std::make_unique<lf::geometry::TriaO1>(std::move(coords)));
		}
	}
	else
	{
		// Swap because indices are provided row-wise.
		QuadrilateralIndices swappedIndices = { indices[0], indices[1], indices[3], indices[2] };
		Eigen::Matrix<double, Eigen::Dynamic, 4> coords(3, 4);
		for(int i = 0; i < 4; ++i)
			coords.col(i) = getPointPosition(meshFactory.GetPoint(swappedIndices[i]));
		meshFactory.AddEntity(lf::base::RefEl::kQuad(), swappedIndices, std::make_unique<lf::geometry::QuadO1>(std::move(coords)));
	}
}

[[nodiscard]] std::shared_ptr<lf::mesh::Mesh> generateCube(Shape shape = Shape::quads)
{
	constexpr unsigned xMax = 2, yMax = 2, zMax = 2;

	auto meshFactory = getGridlikeMeshFactory(xMax, yMax, zMax);
	auto pointMapping = getGridlikePointMapping(xMax, yMax, zMax);
	auto add = castMap<QuadrilateralIndices>(std::bind_front(addQuadrilateral, std::ref(meshFactory), shape));

	add(pointMapping(0, 0, 0), pointMapping(0, 0, 1), pointMapping(0, 1, 0), pointMapping(0, 1, 1));
	add(pointMapping(1, 0, 0), pointMapping(1, 0, 1), pointMapping(1, 1, 0), pointMapping(1, 1, 1));
	add(pointMapping(0, 0, 0), pointMapping(0, 0, 1), pointMapping(1, 0, 0), pointMapping(1, 0, 1));
	add(pointMapping(0, 1, 0), pointMapping(0, 1, 1), pointMapping(1, 1, 0), pointMapping(1, 1, 1));
	add(pointMapping(0, 0, 0), pointMapping(0, 1, 0), pointMapping(1, 0, 0), pointMapping(1, 1, 0));
	add(pointMapping(0, 0, 1), pointMapping(0, 1, 1), pointMapping(1, 0, 1), pointMapping(1, 1, 1));

	return meshFactory.Build();
}

[[nodiscard]] std::shared_ptr<lf::mesh::Mesh> generateMoebiusStrip(Shape shape = Shape::quads)
{
	constexpr unsigned xMax = 2, yMax = 2, zMax = 2;

	auto meshFactory = getGridlikeMeshFactory(xMax, yMax, zMax);
	auto pointMapping = getGridlikePointMapping(xMax, yMax, zMax);
	auto add = castMap<QuadrilateralIndices>(std::bind_front(addQuadrilateral, std::ref(meshFactory), shape));

	add(pointMapping(0, 0, 0), pointMapping(0, 0, 1), pointMapping(0, 1, 0), pointMapping(0, 1, 1));
	add(pointMapping(1, 0, 0), pointMapping(1, 0, 1), pointMapping(1, 1, 0), pointMapping(1, 1, 1));
	add(pointMapping(0, 0, 0), pointMapping(0, 0, 1), pointMapping(1, 0, 0), pointMapping(1, 0, 1));
	// add(pointMapping(0, 1, 0), pointMapping(0, 1, 1), pointMapping(1, 1, 0), pointMapping(1, 1, 1)); // cube
	add(pointMapping(0, 1, 0), pointMapping(0, 1, 1), pointMapping(1, 1, 1), pointMapping(1, 1, 0)); // twisted (last two vertices swapped)

	return meshFactory.Build();
}

[[nodiscard]] std::shared_ptr<lf::mesh::Mesh> generate3Handlebody(unsigned genus, Shape shape = Shape::quads)
{
	const unsigned xMax = 2 + 2 * genus, yMax = 4, zMax = 2;

	auto meshFactory = getGridlikeMeshFactory(xMax, yMax, zMax);
	auto pointMapping = getGridlikePointMapping(xMax, yMax, zMax);
	auto add = castMap<QuadrilateralIndices>(std::bind_front(addQuadrilateral, std::ref(meshFactory), shape));

	// construct the faces
	{
		QuadrilateralIndices indices;

		// x-extremal sides
		for(int x = 0; x < xMax; x += xMax - 1)
		{
			for(int y = 0; y < yMax - 1; ++y)
			{
				for(int zd = 0; zd < 2; ++zd)
				{
					for(int yd = 0; yd < 2; ++yd)
					{
						indices[zd * 2 + yd] = pointMapping(x, y + yd, zd);
					}
				}
				add(indices);
			}
		}

		// y-extremal sides
		for(int x = 0; x < xMax - 1; ++x)
		{
			for(int y = 0; y < yMax; y += yMax - 1)
			{
				for(int xd = 0; xd < 2; ++xd)
				{
					for(int zd = 0; zd < 2; ++zd)
					{
						indices[xd * 2 + zd] = pointMapping(x + xd, y, zd);
					}
				}
				add(indices);
			}
		}

		// z-extremal sides
		for(int x = 0; x < xMax - 1; ++x)
		{
			for(int y = 0; y < yMax - 1; ++y)
			{
				const bool makeHole = x % 2 == 1 && y == 1;
				if(makeHole)
				{
					for(int xd = 0; xd < 2; ++xd)
					{
						for(int zd = 0; zd < 2; ++zd)
						{
							for(int yd = 0; yd < 2; ++yd)
							{
								indices[zd * 2 + yd] = pointMapping(x + xd, y + yd, zd);
							}
						}
						add(indices);
					}
					for(int yd = 0; yd < 2; ++yd)
					{
						for(int xd = 0; xd < 2; ++xd)
						{
							for(int zd = 0; zd < 2; ++zd)
							{
								indices[xd * 2 + zd] = pointMapping(x + xd, y + yd, zd);
							}
						}
						add(indices);
					}
				}
				else
				{
					for(int z = 0; z < zMax; ++z)
					{
						for(int xd = 0; xd < 2; ++xd)
						{
							for(int yd = 0; yd < 2; ++yd)
							{
								indices[xd * 2 + yd] = pointMapping(x + xd, y + yd, z);
							}
						}
						add(indices);
					}
				}
			}
		}
	}

	return meshFactory.Build();
}

[[nodiscard]] std::shared_ptr<lf::mesh::Mesh> generateTriangle()
{
	lf::mesh::hybrid2d::MeshFactory factory(3);

	Eigen::VectorXd p(3);
	p << 1, 0, 0;
	factory.AddPoint(p);
	p << 0, 1, 0;
	factory.AddPoint(p);
	p << 0, 0, 1;
	factory.AddPoint(p);

	PointIndex indices[3] = {0, 1, 2};
	nonstd::span<const PointIndex> span(indices, 3);

	Eigen::Matrix<double, Eigen::Dynamic, 3> coords(3, 3);
	for(int i = 0; i < 3; ++i)
		coords.col(i) = getPointPosition(factory.GetPoint(span[i]));

	factory.AddEntity(lf::base::RefEl::kTria(), span, std::make_unique<lf::geometry::TriaO1>(std::move(coords)));

	return factory.Build();
}

[[nodiscard]] std::shared_ptr<lf::mesh::Mesh> generateTorus(double R = 3.0, double r = 1.0, unsigned thetaSamples = 30, unsigned phiSamples = 100, Shape shape = Shape::quads)
{
	lf::mesh::hybrid2d::MeshFactory meshFactory(3);
	auto pointMapping = getToroidalPointMapping(thetaSamples, phiSamples);
	auto add = castMap<QuadrilateralIndices>(std::bind_front(addQuadrilateral, std::ref(meshFactory), shape));

	// add all points
	{
		Eigen::VectorXd c(3);
		for(unsigned theta = 0; theta < thetaSamples; ++theta)
		{
			for(unsigned phi = 0; phi < phiSamples; ++phi)
			{
				using std::sin, std::cos;

				const double t = 2 * M_PI * theta / thetaSamples;
				const double p = 2 * M_PI * phi / phiSamples;
				c[0] = (R + r * cos(t)) * cos(p);
				c[1] = (R + r * cos(t)) * sin(p);
				c[2] = r * sin(t);
				meshFactory.AddPoint(c);
			}
		}
	}

	// construct the faces
	for(unsigned theta = 0; theta < thetaSamples; ++theta)
	{
		for(unsigned phi = 0; phi < phiSamples; ++phi)
		{
			QuadrilateralIndices indices;
			for(unsigned thetaD = 0; thetaD < 2; ++thetaD)
				for(unsigned phiD = 0; phiD < 2; ++phiD)
					indices[thetaD * 2 + phiD] = pointMapping(theta + thetaD, phi + phiD);
			add(indices);
		}
	}

	return meshFactory.Build();
}


// Small cube with two magnetic ports.
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getCubeExample()
{
	const auto lfMesh = generateCube(Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		if(f == 0)
			return { magneticPort, 0 };
		else if(f == 1)
			return { magneticPort, 1 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Small torus with two electric and two magnetic ports.
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getSmallTorusExample1()
{
	constexpr unsigned thetaSamples = 10;
	constexpr unsigned phiSamples = 10;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		const auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		if(in(theta, 2, 3) && in(phi, 0, 1))
			return { electricPort, 0 };
		else if(in(theta, 8, 9) && in(phi, 0, 1))
			return { electricPort, 1 };
		else if(theta != 1 && phi == 4)
			return { magneticPort, 0 };
		else if(theta != 0 && phi == 8)
			return { magneticPort, 1 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Small torus with a tricky port-path overlap.
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getSmallTorusExample2()
{
	constexpr unsigned thetaSamples = 10;
	constexpr unsigned phiSamples = 10;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		const auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		if(in(theta, 2, 9) && in(phi, 8, 9))
			return { electricPort, 0 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Big torus with a number of ports and overlaps.
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getBigTorusExample()
{
	constexpr unsigned thetaSamples = 30;
	constexpr unsigned phiSamples = 100;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		const auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		if(in(theta, 1, 2) && in(phi, 2, 4))
			return { electricPort, 0 };
		else if(in(theta, 5, 7) && in(phi, 5, 7))
			return { electricPort, 1 };
		else if(in(theta, 10, 20) && in(phi, 40, 60))
			return { electricPort, 2 };
		else if(in(theta, 5, 7) && (phi > 97 || phi < 3))
			return { electricPort, 3 };
		else if(in(theta, 6, 8) && in(phi, 70, 80))
			return { magneticPort, 0 };
		else if(theta != 3 && phi == 20)
			return { magneticPort, 1 };
		else if(!in(theta, 3, 4) && phi == 85)
			return { electricPort, 4 };
		else if(theta != 2 && phi == 90)
			return { magneticPort, 2 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Big torus with many small electric ports (illustrated in the thesis).
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getBigElectricTorusExample()
{
	constexpr unsigned thetaSamples = 30;
	constexpr unsigned phiSamples = 100;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		theta /= 3;
		phi /= 3;
		if(in(theta, 2, 4) && theta % 2 == 0 && in(phi, 0, 18) && phi % 2 == 0)
			return { electricPort, int(theta * thetaSamples + phi) };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Big torus with one electric and one magnetic port with the topological cycle intersecting both.
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getNearFailTorusExample()
{
	constexpr unsigned thetaSamples = 15;
	constexpr unsigned phiSamples = 100;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [](unsigned f) -> PortInformation
	{
		auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		if(!in(theta, 11, 12) && in(phi, 80, 82))
			return { electricPort, 0 };
		else if(in(theta, 13, 15) && in(phi, 89, 91))
			return { magneticPort, 0 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

// Big torus with one electric port where the naive intersection resolution algorithm would fail (as explained in the thesis).
[[nodiscard]] std::pair<Mesh, std::shared_ptr<lf::mesh::Mesh>> getCounterexampleTorusExample()
{
	constexpr unsigned thetaSamples = 15;
	constexpr unsigned phiSamples = 100;
	const auto lfMesh = generateTorus(3.0, 1.0, thetaSamples, phiSamples, Shape::quads);
	const auto portMap = [=](unsigned f) -> PortInformation
	{
		auto [theta, phi] = getInverseToroidalPointMapping(phiSamples)(f);
		if(!in(theta, thetaSamples - 2, thetaSamples) && in(phi, phiSamples - 2, phiSamples))
			return { electricPort, 0 };
		else
			return {};
	};
	Mesh mesh(*lfMesh, portMap);
	return { std::move(mesh), std::move(lfMesh) };
}

} // namespace projects::cohomology

#endif // !COHOMOLOGY_TEST_EXAMPLES_H_INCLUDED

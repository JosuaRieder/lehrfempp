/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#ifndef COHOMOLOGY_VERIFICATION_H_INCLUDED
#define COHOMOLOGY_VERIFICATION_H_INCLUDED

#include <cohomology.h>
#include <utility.h>

#include <Eigen/Dense>
#include <quiver.hpp>

#include <type_traits>
#include <iterator>

/**
 * @brief Contains the verification functionality of the cohomology subproject.
 */
namespace projects::cohomology::verification
{

/**
 * @brief Verifies a tentative set of generators.
 * Refer to the thesis for more information.
 * @param mesh The mesh to which the generators pertain.
 * @param generators A forward range of tentative generators.
 * @return `true` iff the verification is successful.
 */
[[nodiscard]] bool verifyGenerators(Mesh const& mesh, ForwardRangeOf<std::is_convertible, EdgeVector> auto const& generators)
{
	using T = double;
	using MatT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	MatT D1 = mesh.getNodeEdgeIncidence().cast<T>().transpose();
	MatT D2 = mesh.getEdgeCellIncidence().cast<T>().transpose();

	const MatT D1D2 = D1 * D2;
	if(!D1D2.isZero())
		return false;

	MatT G(mesh.edges(), std::ranges::distance(generators));
	for(std::size_t i = 0; auto const& generator : generators)
		G.col(i++) = generator.vector().template cast<T>();

	removeExcludedSlices<0, 0>(mesh, D1);
	removeExcludedSlices<1, 1>(mesh, D1);

	removeExcludedSlices<1, 0>(mesh, D2);
	removeExcludedSlices<2, 1>(mesh, D2);

	removeExcludedSlices<1, 0>(mesh, G);

	if(!(D1 * D2).isZero())
		return false;

	auto const& B1 = D1;
	auto const& B2 = D2;

	assert(B2.rows() == G.rows());

	MatT X(B2.rows(), B2.cols() + G.cols());
	X << B2, G;

	MatT Y = B1.transpose() * B1 + X * X.transpose();

	if(!isPSD(Y))
		return false;

	return true;
}
/**
 * @brief Verifies a tentative set of generators.
 * Refer to the thesis for more information.
 * @param mesh The mesh to which the generators pertain.
 * @param relativeCohomology The computed relative cohomology.
 * @return `true` iff the verification is successful.
 */
[[nodiscard]] bool verifyGenerators(Mesh const& mesh, RelativeCohomology const& relativeCohomology)
{
	return verifyGenerators(mesh, relativeCohomology.cycles);
}

} // namespace projects::cohomology::verification

#endif // !defined(COHOMOLOGY_VERIFICATION_H_INCLUDED)

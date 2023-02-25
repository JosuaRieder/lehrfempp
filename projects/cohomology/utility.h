/**
 * @file
 * @brief Bachelor Thesis "Constructing Generators of Cohomology Classes on Surfaces"
 * @author Josua Rieder
 * @date 2022
 * @copyright ETH Zurich
 */

#ifndef COHOMOLOGY_UTILITY_H_INCLUDED
#define COHOMOLOGY_UTILITY_H_INCLUDED

#include <lf/mesh/mesh.h>
#include <Eigen/Dense>

#include <type_traits>
#include <functional>
#include <ranges>
#include <iterator>

namespace projects::cohomology
{

/**
 * @brief Checks whether a value is in the closed interval `[min, max]`.
 * @param v The value.
 * @param min The lower bound.
 * @param max The upper bound.
 * @return `true` iff `min <= v && v <= max`.
 */
template<typename T, typename U, typename V>
[[nodiscard]] constexpr bool in(T const& v, U const& min, V const& max) noexcept
{
	assert(min <= max);
	return min <= v && v <= max;
}

/**
 * @brief Extracts the position of an object of `lf::geometry::Point`.
 * @param v The point.
 * @return The point's position.
 */
[[nodiscard]] Eigen::Vector3d getPointPosition(lf::geometry::Point const& point)
{
	return point.Global(Eigen::MatrixXd::Zero(0, 1));
}

/**
 * @brief Combines the hashes of individual objects in an order-sensitive way, generating one single hash for the entire collection of objects.
 * Adapted from `boost::hash_combine`.
 * @return The combined hash for all arguments.
 */
template<typename... Args>
[[nodiscard]] std::size_t hashTogether(Args&&... args)
{
	std::size_t result = 0;
	((result = result ^ (std::hash<std::decay_t<Args>>{}(args) + 0x9e3779b9 + (result << 6) + (result >> 2))), ...);
	return result;
}

/**
 * @brief A type trait that is `true` for Eigen types (derived from `Eigen::EigenBase`.
 * @see EigenType
 */
template<typename T>
class IsEigenType
{
	static std::true_type test(Eigen::EigenBase<T> const*);
	static std::false_type test(...);

public:
	static constexpr bool value = decltype(test((T*)(nullptr)))::value;
};

/**
 * @brief A concept for Eigen types.
 * @see IsEigenType
 */
template<typename T>
concept EigenType = IsEigenType<T>::value;

/**
 * @brief A concept for input ranges that are restricted in their value types.
 * @tparam Predicate The meta predicate that must be fulfilled by the range's value type.
 * @tparam Trailing The (optional) trailing arguments to the predicate (the first one being the range type).
 */
template<typename Range, template<typename...> class Predicate, typename... Trailing>
concept InputRangeOf = std::ranges::input_range<Range> && Predicate<typename std::ranges::range_value_t<Range>, Trailing...>::value;
/**
 * @brief A concept for forward ranges that are restricted in their value types.
 * @tparam Predicate The meta predicate that must be fulfilled by the range's value type.
 * @tparam Trailing The (optional) trailing arguments to the predicate (the first one being the range type).
 */
template<typename Range, template<typename...> class Predicate, typename... Trailing>
concept ForwardRangeOf = std::ranges::forward_range<Range> && Predicate<typename std::ranges::range_value_t<Range>, Trailing...>::value;

/**
 * @brief Returns whether a given matrix is positive-semidefinite using an LLT decomposition.
 * @param matrix The matrix whose positive-semidefiniteness is to be verified.
 * @return `true` iff `matrix` is PSD.
 */
[[nodiscard]] bool isPSD(EigenType auto const& matrix)
{
	// Matrix has to be approximately equal to its transpose (symmetry).
	if(!matrix.isApprox(matrix.transpose()))
		return false;

	// Compute an LLT decomposition which fails if the matrix is not positive-definite.
	auto llt = matrix.llt();
	if(llt.info() == Eigen::NumericalIssue)
		return false;

	return true;
}

} // namespace projects::cohomology

#endif // !defined(COHOMOLOGY_UTILITY_H_INCLUDED)

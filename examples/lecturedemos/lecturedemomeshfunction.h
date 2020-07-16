#ifndef LF_LD_MF_H
#define LF_LD_MF_H
/**
 * @file
 * @brief Functions for simple LehrFEM++ demos + sample codes
 * @author Ralf Hiptmair
 * @date  July 2020
 * @copyright MIT License
 */

#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/refinement/refinement.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>

namespace lecturedemo {

/** @brief Computation of integral expression by means of mesh functions
 *
 * @tparam COEFFFUNCTION a functor type std::function<R(Eigen::Vector2d)>, where
 * R is a type that can left-multiply a vector, that is, either a scalar or a
 * matrix
 * @tparam VECTORFIELD a functor tye for a vector field in 2D
 * @param fe_space a Lagrangian finite-element space object
 * @param coeff_vec basis expansion coefficient vector for a finite element
 * function
 * @param coefficient coefficient function (either scalar or matrix valued
 * @param vectorfield vector field
 *
 *
 *
 */
template <typename COEFFFUNCTION, typename VECTORFIELD>
double integrateCoeffgradUhVf(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const Eigen::VectorXd& coeff_vec, COEFFFUNCTION&& coefficient,
    VECTORFIELD&& vectorfield) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{fe_space->Mesh()};
  // Coefficient function and weight function
  const lf::mesh::utils::MeshFunctionGlobal mf_coefficient(coefficient);
  const lf::mesh::utils::MeshFunctionGlobal mf_vectorfield(vectorfield);
  // Build a MeshFunction representing the gradient of the finite element
  // solution
  const lf::uscalfe::MeshFunctionGradFE mf_grad(fe_space, coeff_vec);
  // Mesh function representing the integrand
  const auto mf_itg{lf::mesh::utils::transpose(mf_coefficient * mf_grad) *
                    mf_vectorfield};
  const double s = lf::uscalfe::IntegrateMeshFunction(
      *mesh_p, mf_itg, [](const lf::mesh::Entity& e) {
        return lf::quad::make_QuadRule(e.RefEl(), 2);
      })(0, 0);
  return s;
}  // end stabFlux

void lecturedemomeshfunction();

}  // namespace lecturedemo

#endif

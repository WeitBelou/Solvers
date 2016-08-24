#ifndef UTILS_HPP
#define UTILS_HPP

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/symmetric_tensor.h>

#include "global.hpp"

namespace Utils {
//begin namespace Utils
using namespace dealii;

inline SymmetricTensor<2, DIM> get_strain(const FEValues<DIM> &fe_values, const size_t shape_func,
                                          const size_t q_point)
{
    SymmetricTensor<2, DIM> strain;

    for (size_t i = 0; i < DIM; ++i) {
        strain[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
    }

    for (size_t i = 0; i < DIM; ++i) {
        for (size_t j = i + 1; j < DIM; ++j) {
            strain[i][j] = (fe_values.shape_grad_component(shape_func, q_point, i)[j]
                            + fe_values.shape_grad_component(shape_func, q_point, j)[i]) / 2;
        }
    }

    return strain;
}

SymmetricTensor<4, DIM> get_stress_strain_tensor(const double lambda, const double mu);

inline SymmetricTensor<2, DIM> get_strain(const std::vector<Tensor<1, DIM>> &grad)
{
    SymmetricTensor<2, DIM> strain;

    Assert(grad.size() == DIM, ExcDimensionMismatch(grad.size(), DIM));

    for (size_t i = 0; i < DIM; ++i) {
        strain[i][i] = grad[i][i];
    }

    for (size_t i = 0; i < DIM; ++i) {
        for (size_t j = i + 1; j < DIM; ++j) {
            strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
        }
    }

    return strain;
}

Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u);

//end namespace Utils
}

#endif // UTILS_HPP

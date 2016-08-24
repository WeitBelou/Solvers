#include "Utils.hpp"

using namespace Utils;

Tensor<2, 3> Utils::get_rotation_matrix(const std::vector<Tensor<1, 3> > &grad_u)
{
    const Point<3> curl(grad_u[2][1] - grad_u[1][2],
                        grad_u[0][2] - grad_u[2][0],
                        grad_u[1][0] - grad_u[0][1]);

    const double angle = std::atan(curl.norm());

    if (angle < 1e-9) {
        static const double rotation[3][3] = {{1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
        static const Tensor<2, 3> rot(rotation);
        return rot;
    }

    const double s = std::sin(angle);
    const double c = std::cos(angle);
    const double t = 1 - c;

    const Point<3> axis = curl / curl.norm();

    const double rotation[3][3] = {
        {
            t *axis[0] *axis[0] + c,
            t *axis[0] *axis[1] + s *axis[2],
            t *axis[0] *axis[2] - s *axis[1]
        },
        {
            t *axis[0] *axis[1] - s *axis[2],
            t *axis[0] *axis[1] + c,
            t *axis[1] *axis[2] + s *axis[0]
        },
        {
            t *axis[0] *axis[2] + s *axis[1],
            t *axis[1] *axis[1] - s *axis[0],
            t *axis[2] *axis[2] + c
        }
    };

    return Tensor<2, 3, double>(rotation);
}

SymmetricTensor<4, DIM> Utils::get_stress_strain_tensor(const double lambda,
                                                                   const double mu)
{
    auto d = [](size_t i, size_t j) {
        return (i == j) ? (1.0) : (0.0);
    };

    SymmetricTensor<4, DIM> stress_strain;
    for (size_t i = 0; i < DIM; ++i) {
        for (size_t j = 0; j < DIM; ++j) {
            for (size_t k = 0; k < DIM; ++k) {
                for (size_t l = 0; l < DIM; ++l) {
                    stress_strain[i][j][k][l] = mu * (d(i, k) * d(j, l) + d(i, l) * d(j, k))
                                                + lambda * d(i, j) * d(k, l);
                }
            }
        }
    }
    return stress_strain;
}

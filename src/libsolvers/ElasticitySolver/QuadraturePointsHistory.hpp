#ifndef QUADRATUREPOINTSHISTORY_HPP
#define QUADRATUREPOINTSHISTORY_HPP

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>

#include "global.hpp"

namespace PointsHistory {

using namespace dealii;

struct PointHistory {
    SymmetricTensor<2, DIM> old_stress;
};

class QuadraturePointsHistory
{
public:
    QuadraturePointsHistory(const FiniteElement<DIM> &fe,
                            const Quadrature<DIM> &quadrature);

    void setup(Triangulation<DIM> &triangulation);
    void update(const DoFHandler<DIM> &dof_handler,
                Vector<double> incremental_displacement,
                const SymmetricTensor<4, DIM> &stress_strain_tensor);
    void get_cell_quadrature_points_data();
private:
    const SmartPointer<const FiniteElement<DIM>> fe;
    const SmartPointer<const Quadrature<DIM>> quadrature;

    std::vector<PointHistory> quadrature_points_history;
};

}
#endif // QUADRATUREPOINTSHISTORY_HPP

#ifndef QUADRATUREPOINTSHISTORY_HPP
#define QUADRATUREPOINTSHISTORY_HPP

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>

#include "src/libsolvers/global.hpp"

namespace ElasticityEquation
{

struct PointHistory
{
    dealii::SymmetricTensor<2, DIM> old_stress;
};

class QuadraturePointsHistory
{
public:
    QuadraturePointsHistory(const dealii::FiniteElement<DIM> &fe,
                            const dealii::Quadrature<DIM> &quadrature);

    void setup(dealii::Triangulation<DIM> &triangulation);
    void update(const dealii::DoFHandler<DIM> &dof_handler,
                dealii::Vector<double> incremental_displacement,
                const dealii::SymmetricTensor<4, DIM> &stress_strain_tensor);
    const std::vector<PointHistory> &
    get_cell_quadrature_points_data(const typename dealii::DoFHandler<DIM>::active_cell_iterator &cell) const;
private:
    const dealii::SmartPointer<const dealii::FiniteElement <DIM>> fe;
    const dealii::SmartPointer<const dealii::Quadrature <DIM>> quadrature;

    std::vector<std::vector<PointHistory>> quadrature_points_history;
};

}
#endif // QUADRATUREPOINTSHISTORY_HPP

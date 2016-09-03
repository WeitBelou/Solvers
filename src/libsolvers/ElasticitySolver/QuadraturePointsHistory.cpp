#include "QuadraturePointsHistory.hpp"
#include "Utils.hpp"

namespace ut = ElasticityEquation;
using namespace dealii;

ElasticityEquation::QuadraturePointsHistory::QuadraturePointsHistory(const FiniteElement<DIM> &fe,
                                                                     const Quadrature<DIM> &quadrature)
    :
    fe(&fe),
    quadrature(&quadrature)
{

}

void ElasticityEquation::QuadraturePointsHistory::setup(Triangulation<DIM> &triangulation)
{
    size_t n_cells = triangulation.n_active_cells();

    triangulation.clear_user_data();

    quadrature_points_history.clear();
    quadrature_points_history.resize(n_cells, std::vector<PointHistory>(quadrature->size()));

    size_t history_index = 0;
    for (auto &&cell : triangulation.active_cell_iterators())
    {
        cell->set_user_index(history_index);
        history_index++;
    }
    Assert(history_index == n_cells, ExcInternalError());
}

void ElasticityEquation::QuadraturePointsHistory::update(const DoFHandler<DIM> &dof_handler,
                                                         Vector<double> incremental_displacement,
                                                         const SymmetricTensor<4, DIM> &stress_strain_tensor)
{
    FEValues<DIM> fe_values(*fe, *quadrature,
                            update_values | update_gradients);
    std::vector<std::vector<Tensor<1, DIM>>>
        displacement_increment_grads(quadrature->size(), std::vector<Tensor<1, DIM>>(DIM));

    for (auto &&cell : dof_handler.active_cell_iterators())
    {
        std::vector<PointHistory> &local_quadrature_points_history =
            quadrature_points_history[cell->user_index()];

        fe_values.reinit(cell);
        fe_values.get_function_gradients(incremental_displacement,
                                         displacement_increment_grads);
        for (size_t q = 0; q < quadrature->size(); ++q)
        {
            const SymmetricTensor<2, DIM> new_stress = local_quadrature_points_history[q].old_stress
                                                       + stress_strain_tensor
                                                         * ut::get_strain(displacement_increment_grads[q]);

            const Tensor<2, DIM> rotation = ut::get_rotation_matrix(displacement_increment_grads[q]);
            const SymmetricTensor<2, DIM> rotated_new_stress = symmetrize(transpose(rotation)
                                                                          * static_cast<Tensor<2, DIM>>(new_stress)
                                                                          * rotation);
            local_quadrature_points_history[q].old_stress = rotated_new_stress;
        }
    }
}

const std::vector<ElasticityEquation::PointHistory>
&ElasticityEquation::QuadraturePointsHistory::get_cell_quadrature_points_data(
    const typename DoFHandler<DIM>::active_cell_iterator &cell) const
{
    return quadrature_points_history[cell->user_index()];
}

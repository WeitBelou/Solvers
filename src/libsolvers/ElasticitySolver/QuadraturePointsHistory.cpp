#include "QuadraturePointsHistory.hpp"
#include "Utils.hpp"

namespace ut = Utils;
PointsHistory::QuadraturePointsHistory::QuadraturePointsHistory(const FiniteElement<DIM> &fe,
                                                                const Quadrature<DIM> &quadrature)
    :
    fe(&fe),
    quadrature(&quadrature)
{

}

void PointsHistory::QuadraturePointsHistory::setup(Triangulation<DIM> &triangulation)
{
    size_t n_cells = triangulation.n_active_cells();

    triangulation.clear_user_data();

    quadrature_points_history.resize(n_cells * quadrature->size());

    size_t history_index = 0;
    for (auto && cell : triangulation.active_cell_iterators()) {
        cell->set_user_pointer(&quadrature_points_history[history_index]);
        history_index += quadrature->size();
    }
    Assert(history_index == quadrature_points_history.size(), ExcInternalError());
}

void PointsHistory::QuadraturePointsHistory::update(const DoFHandler<DIM> &dof_handler,
                                                    Vector<double> incremental_displacement,
                                                    const SymmetricTensor<4, DIM> &stress_strain_tensor)
{
    FEValues<DIM> fe_values(*fe, *quadrature,
                            update_values | update_gradients);
    std::vector<std::vector<Tensor<1, DIM>>>
    displacement_increment_grads(quadrature->size(), std::vector<Tensor<1, DIM>>(DIM));

    for (auto && cell : dof_handler.active_cell_iterators()) {
        PointHistory *local_quadrature_points_history = reinterpret_cast<PointHistory *>
                                                        (cell->user_pointer());
        Assert(local_quadrature_points_history >= &quadrature_points_history.front(), ExcInternalError());
        Assert(local_quadrature_points_history < &quadrature_points_history.back(), ExcInternalError());

        fe_values.reinit(cell);
        fe_values.get_function_gradients(incremental_displacement,
                                         displacement_increment_grads);
        for (size_t q = 0; q < quadrature->size(); ++q) {
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

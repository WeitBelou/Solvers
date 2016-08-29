//
// Created by ivan on 10.08.16.
//

#include "ElasticitySolver.hpp"

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>

using namespace ElasticitySolver;

const SymmetricTensor<4, DIM> TopLevel::stress_strain_tensor = ut::get_stress_strain_tensor(
                                                                   9.695e10,
                                                                   7.617e10);

TopLevel::TopLevel(Triangulation<DIM> &triangulation,
                   const FESystem<DIM> &fe,
                   const QGauss<DIM> &quadrature,
                   const Function<DIM> &body_force,
                   BoundaryConditions::FunctionTimeBoundaryConditions &boundary_conditions)
    :
    triangulation(&triangulation),
    fe(&fe),
    quadrature(&quadrature),
    dof_handler(triangulation),
    body_force(&body_force),
    boundary_conditions(&boundary_conditions),
    quadrature_points_history(fe, quadrature)
{

}

TopLevel::~TopLevel()
{
    dof_handler.clear();
}

void TopLevel::run(double timestep, double end_time)
{
    present_time = 0.0;
    present_timestep = timestep;
    this->end_time = end_time;

    timestep_no = 0;

    boundary_conditions->reinit(present_time);

    do_initial_timestep();

    while (present_time < end_time) {
        do_timestep();
    }
}

void TopLevel::do_initial_timestep()
{
    dof_handler.distribute_dofs(*fe);
    quadrature_points_history.setup(*triangulation);

    do_timestep();
}

void TopLevel::do_timestep()
{
    present_time += present_timestep;
    ++timestep_no;

    std::cout << "Timestep " << timestep_no << " at time " << present_time << std::endl;

    if (present_time > end_time) {
        present_timestep -= (present_time - end_time);
        present_time = end_time;
    }

    solve_timestep();
    move_mesh();
    output_results();

    std::cout << std::endl;
}

void TopLevel::solve_timestep()
{
    std::cout << "    Assembling system..." << std::flush;

    LinearSystem linear_system(dof_handler);
    assemble_linear_system(linear_system);
    std::cout << std::endl;

    std::cout << "    Solving linear system..." << std::flush;
    linear_system.solve(incremental_displacement);
    std::cout << std::endl;

    std::cout << "    Updating quadrature points history..." << std::flush;
    quadrature_points_history.update(dof_handler, incremental_displacement, stress_strain_tensor);
    std::cout << std::endl;
}

void TopLevel::assemble_linear_system(TopLevel::LinearSystem &linear_system)
{
    auto assemble_system_func = [this](typename DoFHandler<DIM>::active_cell_iterator cell,
                                       AssemblyScratchData & scratch_data,
    AssemblyCopyData & copy_data) {
        this->local_assemble_system(cell, scratch_data, copy_data);
    };

    auto copy_local_to_global_func = [this, &linear_system](const AssemblyCopyData & copy_data) {
        this->copy_local_to_global(copy_data, linear_system);
    };

    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    assemble_system_func,
                    copy_local_to_global_func,
                    AssemblyScratchData(*fe, *quadrature),
                    AssemblyCopyData());

    incremental_displacement.reinit(dof_handler.n_dofs());

    boundary_conditions->update(present_timestep);
    std::map<types::global_dof_index, double> boundary_values
        = boundary_conditions->interpolate(dof_handler);

    MatrixTools::apply_boundary_values(boundary_values,
                                       linear_system.matrix,
                                       incremental_displacement,
                                       linear_system.rhs,
                                       false);
}

void TopLevel::move_mesh()
{
    std::cout << "    Moving mesh..." << std::flush;

    std::vector<bool> vertex_touched(triangulation->n_vertices(), false);
    for (auto && cell : dof_handler.active_cell_iterators()) {
        for (size_t v = 0; v < GeometryInfo<DIM>::vertices_per_cell; ++v) {
            if (vertex_touched[cell->vertex_index(v)] == false) {
                vertex_touched[cell->vertex_index(v)] = true;

                Point<DIM> vertex_displacement;
                for (size_t d = 0; d < DIM; ++d) {
                    vertex_displacement[d] = incremental_displacement(cell->vertex_dof_index(v, d));
                }
                cell->vertex(v) += vertex_displacement;
            }
        }
    }

    std::cout << std::endl;
}

void TopLevel::local_assemble_system(const typename DoFHandler<DIM>::active_cell_iterator &cell,
                                     TopLevel::AssemblyScratchData &scratch_data,
                                     TopLevel::AssemblyCopyData &copy_data) const
{
    size_t dofs_per_cell = fe->dofs_per_cell;
    size_t n_q_points = quadrature->size();

    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    copy_data.cell_rhs.reinit(dofs_per_cell);

    copy_data.local_dofs_indices.resize(dofs_per_cell);
    scratch_data.fe_values.reinit(cell);

    for (size_t i = 0; i < dofs_per_cell; ++i) {
        for (size_t j = 0; j < dofs_per_cell; ++j) {
            for (size_t q = 0; q < n_q_points; ++q) {
                const SymmetricTensor<2, DIM>
                eps_phi_i = ut::get_strain(scratch_data.fe_values, i, q),
                eps_phi_j = ut::get_strain(scratch_data.fe_values, j, q);

                copy_data.cell_matrix(i, j) += (eps_phi_i * stress_strain_tensor * eps_phi_j
                                                * scratch_data.fe_values.JxW(q));
            }
        }
    }

    std::vector<Vector<double>> body_force_values(n_q_points, Vector<double>(DIM));
    body_force->vector_value_list(scratch_data.fe_values.get_quadrature_points(), body_force_values);

    const std::vector<ph::PointHistory> &local_quadrature_points_data =
        quadrature_points_history.get_cell_quadrature_points_data(cell);

    for (size_t i = 0; i < dofs_per_cell; ++i) {
        size_t component_i = fe->system_to_component_index(i).first;

        for (size_t q = 0; q < n_q_points; ++q) {
            const SymmetricTensor<2, DIM> &old_stress = local_quadrature_points_data[q].old_stress;
            copy_data.cell_rhs(i) += (body_force_values[q](component_i)
                                      * scratch_data.fe_values.shape_value(i, q)
                                      -
                                      old_stress *
                                      ut::get_strain(scratch_data.fe_values, i, q))
                                     * scratch_data.fe_values.JxW(q);
        }
    }

    cell->get_dof_indices(copy_data.local_dofs_indices);
}

void
TopLevel::copy_local_to_global(const TopLevel::AssemblyCopyData &copy_data,
                               TopLevel::LinearSystem &linear_system) const
{
    linear_system.hanging_node_constraints.distribute_local_to_global(copy_data.cell_matrix,
                                                                      copy_data.cell_rhs,
                                                                      copy_data.local_dofs_indices,
                                                                      linear_system.matrix,
                                                                      linear_system.rhs);
}

void TopLevel::output_results() const
{
    DataOut<DIM> data_out;
    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names;

    solution_names.push_back("delta_x");
    solution_names.push_back("delta_y");
    solution_names.push_back("delta_z");

    data_out.add_data_vector(incremental_displacement, solution_names);

    data_out.build_patches();
    std::string filename = "solution" + Utilities::int_to_string(timestep_no, 4) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
}

TopLevel::LinearSystem::LinearSystem(const DoFHandler<DIM> &dof_handler)
{
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, hanging_node_constraints);
    hanging_node_constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    hanging_node_constraints.condense(dsp);
    sparsity_pattern.copy_from(dsp);

    matrix.reinit(sparsity_pattern);
    rhs.reinit(dof_handler.n_dofs());
}

void TopLevel::LinearSystem::solve(Vector<double> &solution) const
{
    SolverControl solver_control(1000, 1e-12);
    SolverCG<> cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(matrix, 1.2);

    cg.solve(matrix, solution, rhs, preconditioner);
    hanging_node_constraints.distribute(solution);
}

TopLevel::AssemblyScratchData::AssemblyScratchData(const FiniteElement<DIM> &fe,
                                                   const Quadrature<DIM> &quadrature)
    :
    fe_values(fe, quadrature,
             update_values | update_gradients |
             update_quadrature_points | update_JxW_values)
{

}

TopLevel::AssemblyScratchData::AssemblyScratchData(const AssemblyScratchData &scratch)
    :
    fe_values(scratch.fe_values.get_fe(),
             scratch.fe_values.get_quadrature(),
             scratch.fe_values.get_update_flags())
{}

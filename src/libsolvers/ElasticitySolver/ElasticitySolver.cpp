//
// Created by ivan on 10.08.16.
//

#include "ElasticitySolver.hpp"

using namespace ElasticitySolver;

const SymmetricTensor<4, DIM> TopLevel::stress_strain_tensor = get_stress_strain_tensor(9.695e10, 7.617e10);

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
    boundary_conditions(&boundary_conditions)
{

}

TopLevel::~TopLevel()
{
    dof_handler.clear();
}

void TopLevel::run()
{
    present_time = 0.0;
    present_timestep = 1.0;
    end_time = 10.0;

    timestep_no = 0;

    boundary_conditions->reinit(present_time);

    do_initial_timestep();

    while (present_time < end_time)
    {
        do_timestep();
    }
}

void TopLevel::do_initial_timestep()
{
    present_time += present_timestep;
    ++timestep_no;

    std::cout << "Timestep " << timestep_no << " at time " << present_time << std::endl;

    for (size_t cycle = 0; cycle < 1; ++cycle)
    {
        std::cout << "    Cycle " << cycle << ':' << std::endl;

        dof_handler.distribute_dofs(*fe);
        std::cout << "    Number of dofs: " << dof_handler.n_dofs() << std::endl;

        if (cycle == 0)
        {
            setup_quadrature_point_history();
        }
        else
        {
            refine_grid();
        }

        std::cout << "    Number of active cells: " << triangulation->n_active_cells() << std::endl;
        solve_timestep();

    }

    move_mesh();
    output_results();

    std::cout << std::endl;
}

void TopLevel::do_timestep()
{
    present_time += present_timestep;
    ++timestep_no;

    std::cout << "Timestep " << timestep_no << " at time " << present_time << std::endl;

    if (present_time > end_time)
    {
        present_timestep -= (present_time - end_time);
        present_time = end_time;
    }

    solve_timestep();
    move_mesh();
    output_results();

    std::cout << std::endl;
}

void TopLevel::refine_grid()
{
    Vector<float> error_per_cell(triangulation->n_active_cells());
    KellyErrorEstimator<DIM>::estimate(dof_handler,
                                       QGauss<DIM - 1>(2),
                                       typename FunctionMap<DIM>::type(),
                                       incremental_displacement,
                                       error_per_cell,
                                       ComponentMask(),
                                       0);

    GridRefinement::refine_and_coarsen_fixed_number(*triangulation, error_per_cell,
                                                    0.35, 0.03);
    triangulation->execute_coarsening_and_refinement();
    setup_quadrature_point_history();
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
    update_quadrature_point_history();
    std::cout << std::endl;
}

void TopLevel::assemble_linear_system(TopLevel::LinearSystem &linear_system)
{
    auto assemble_system_func = [this](typename DoFHandler<DIM>::active_cell_iterator cell,
                                       AssemblyScratchData &scratch_data,
                                       AssemblyCopyData &copy_data)
    {
        this->local_assemble_system(cell, scratch_data, copy_data);
    };

    auto copy_local_to_global_func = [this, &linear_system](const AssemblyCopyData &copy_data)
    {
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

void TopLevel::setup_quadrature_point_history()
{
    size_t n_cells = triangulation->n_active_cells();

    triangulation->clear_user_data();
    {
        std::vector<PointHistory> tmp;
        tmp.swap(quadrature_point_history);
    }
    quadrature_point_history.resize(n_cells * quadrature->size());

    size_t history_index = 0;
    for (auto &&cell : triangulation->active_cell_iterators())
    {
        cell->set_user_pointer(&quadrature_point_history[history_index]);
        history_index += quadrature->size();
    }
    Assert(history_index == quadrature_point_history.size(), ExcInternalError());
}

void TopLevel::update_quadrature_point_history()
{
    FEValues<DIM> fe_values(*fe, *quadrature,
                            update_values | update_gradients);
    std::vector<std::vector<Tensor<1, DIM>>>
        displacement_increment_grads(quadrature->size(), std::vector<Tensor<1, DIM>>(DIM));
    for (auto &&cell : dof_handler.active_cell_iterators())
    {
        PointHistory *local_quadrature_points_history = reinterpret_cast<PointHistory *>(cell->user_pointer());
        Assert(local_quadrature_points_history >= &quadrature_point_history.front(), ExcInternalError());
        Assert(local_quadrature_points_history < &quadrature_point_history.back(), ExcInternalError());

        fe_values.reinit(cell);
        fe_values.get_function_gradients(incremental_displacement,
                                         displacement_increment_grads);
        for (size_t q = 0; q < quadrature->size(); ++q)
        {
            const SymmetricTensor<2, DIM> new_stress = local_quadrature_points_history[q].old_stress
                                                       + stress_strain_tensor
                                                         * get_strain(displacement_increment_grads[q]);

            const Tensor<2, DIM> rotation = get_rotation_matrix(displacement_increment_grads[q]);
            const SymmetricTensor<2, DIM> rotated_new_stress = symmetrize(transpose(rotation)
                                                                          * static_cast<Tensor<2, DIM>>(new_stress)
                                                                          * rotation);
            local_quadrature_points_history[q].old_stress = rotated_new_stress;
        }
    }
}
void TopLevel::move_mesh()
{
    std::cout << "    Moving mesh..." << std::flush;

    std::vector<bool> vertex_touched(triangulation->n_vertices(), false);
    for (auto &&cell : dof_handler.active_cell_iterators())
    {
        for (size_t v = 0; v < GeometryInfo<DIM>::vertices_per_cell; ++v)
        {
            if (vertex_touched[cell->vertex_index(v)] == false)
            {
                vertex_touched[cell->vertex_index(v)] = true;

                Point<DIM> vertex_displacement;
                for (size_t d = 0; d < DIM; ++d)
                {
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

    for (size_t i = 0; i < dofs_per_cell; ++i)
    {
        for (size_t j = 0; j < dofs_per_cell; ++j)
        {
            for (size_t q = 0; q < n_q_points; ++q)
            {
                const SymmetricTensor<2, DIM>
                    eps_phi_i = get_strain(scratch_data.fe_values, i, q),
                    eps_phi_j = get_strain(scratch_data.fe_values, j, q);

                copy_data.cell_matrix(i, j) += (eps_phi_i * stress_strain_tensor * eps_phi_j
                                                * scratch_data.fe_values.JxW(q));
            }
        }
    }

    std::vector<Vector<double>> body_force_values(n_q_points, Vector<double>(DIM));
    body_force->vector_value_list(scratch_data.fe_values.get_quadrature_points(), body_force_values);

    const PointHistory *local_quadrature_points_data = reinterpret_cast<PointHistory *>(cell->user_pointer());

    for (size_t i = 0; i < dofs_per_cell; ++i)
    {
        size_t component_i = fe->system_to_component_index(i).first;

        for (size_t q = 0; q < n_q_points; ++q)
        {
            const SymmetricTensor<2, DIM> &old_stress = local_quadrature_points_data[q].old_stress;
            copy_data.cell_rhs(i) += (body_force_values[q](component_i)
                                      * scratch_data.fe_values.shape_value(i, q)
                                      -
                                      old_stress *
                                      get_strain(scratch_data.fe_values, i, q))
                                     * scratch_data.fe_values.JxW(q);
        }
    }

    cell->get_dof_indices(copy_data.local_dofs_indices);
}

void
TopLevel::copy_local_to_global(const TopLevel::AssemblyCopyData &copy_data, TopLevel::LinearSystem &linear_system) const
{
    linear_system.hanging_node_constraints.distribute_local_to_global(copy_data.cell_matrix,
                                                                      copy_data.cell_rhs,
                                                                      copy_data.local_dofs_indices,
                                                                      linear_system.matrix,
                                                                      linear_system.rhs);
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
{}

TopLevel::AssemblyScratchData::AssemblyScratchData(const AssemblyScratchData &scratch)
    :
    fe_values(scratch.fe_values.get_fe(),
              scratch.fe_values.get_quadrature(),
              scratch.fe_values.get_update_flags())
{}

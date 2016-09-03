//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_ELASTICITYSOLVER_HPP
#define SOLVERS_ELASTICITYSOLVER_HPP

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/base/work_stream.h>

#include <iostream>
#include <fstream>

#include "src/libsolvers/global.hpp"
#include "src/libsolvers/BoundaryConditions/BoundaryConditions.hpp"
#include "QuadraturePointsHistory.hpp"
#include "Utils.hpp"

namespace bc = BoundaryConditions;

namespace ElasticityEquation
{

class ElasticitySolver
{
public:
    ElasticitySolver(dealii::Triangulation<DIM> &triangulation,
                     const dealii::FESystem<DIM> &fe,
                     const dealii::QGauss<DIM> &quadrature,
                     const dealii::Function<DIM> &body_force,
                     bc::FunctionBoundaryConditions &boundary_conditions);
    ~ElasticitySolver();
    void run(double timestep, double end_time);

private:
    const dealii::SmartPointer<dealii::Triangulation<DIM>> triangulation;
    const dealii::SmartPointer<const dealii::FESystem<DIM>> fe;
    const dealii::SmartPointer<const dealii::QGauss<DIM>> quadrature;

    dealii::DoFHandler<DIM> dof_handler;

    const dealii::SmartPointer<const dealii::Function<DIM>> body_force;
    const dealii::SmartPointer<bc::FunctionBoundaryConditions> boundary_conditions;

    dealii::Vector<double> incremental_displacement;

    QuadraturePointsHistory quadrature_points_history;
    static const dealii::SymmetricTensor<4, DIM> stress_strain_tensor;

    double present_time;
    double present_timestep;
    double end_time;

    size_t timestep_no;

    void do_initial_timestep();
    void do_timestep();

    void solve_timestep();
    void move_mesh();
    void output_results() const;

    class LinearSystem
    {
    public:
        LinearSystem(const dealii::DoFHandler<DIM> &dof_handler);

        void solve(dealii::Vector<double> &solution) const;

        dealii::ConstraintMatrix hanging_node_constraints;
        dealii::SparsityPattern sparsity_pattern;
        dealii::SparseMatrix<double> matrix;
        dealii::Vector<double> rhs;
    };

    struct AssemblyScratchData
    {
        AssemblyScratchData(const dealii::FiniteElement<DIM> &fe,
                            const dealii::Quadrature<DIM> &quadrature);
        AssemblyScratchData(const AssemblyScratchData &scratch);
        dealii::FEValues<DIM> fe_values;
    };

    struct AssemblyCopyData
    {
        dealii::FullMatrix<double> cell_matrix;
        dealii::Vector<double> cell_rhs;

        std::vector<dealii::types::global_dof_index> local_dofs_indices;
    };

    void assemble_linear_system(LinearSystem &linear_system);

    void local_assemble_system(const typename dealii::DoFHandler<DIM>::active_cell_iterator &cell,
                               AssemblyScratchData &scratch_data,
                               AssemblyCopyData &copy_data) const;

    void copy_local_to_global(const AssemblyCopyData &copy_data,
                              LinearSystem &linear_system) const;
};
}


#endif //SOLVERS_ELASTICITYSOLVER_HPP

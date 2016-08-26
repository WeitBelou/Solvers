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

//#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>

#include <iostream>
#include <fstream>

#include "global.hpp"
#include "BoundaryConditions.hpp"
#include "QuadraturePointsHistory.hpp"
#include "Utils.hpp"

namespace ElasticitySolver {
//begin namespace ElasticitySolver
using namespace dealii;
namespace bc = BoundaryConditions;
namespace ut = Utils;
namespace ph = PointsHistory;

class TopLevel
{
public:
    TopLevel(Triangulation <DIM> &triangulation,
             const FESystem <DIM> &fe,
             const QGauss <DIM> &quadrature,
             const Function <DIM> &body_force,
             BoundaryConditions::FunctionTimeBoundaryConditions &boundary_conditions);
    ~TopLevel();
    void run();

private:
    const SmartPointer<Triangulation<DIM>> triangulation;
    const SmartPointer<const FESystem<DIM>> fe;
    const SmartPointer<const QGauss<DIM>> quadrature;

    DoFHandler<DIM> dof_handler;

    const SmartPointer<const Function<DIM>> body_force;
    const SmartPointer<bc::FunctionTimeBoundaryConditions> boundary_conditions;

    Vector<double> incremental_displacement;

    ph::QuadraturePointsHistory quadrature_points_history;
    static const SymmetricTensor<4, DIM> stress_strain_tensor;

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
        LinearSystem(const DoFHandler<DIM> &dof_handler);

        void solve(Vector<double> &solution) const;

        ConstraintMatrix hanging_node_constraints;
        SparsityPattern sparsity_pattern;
        SparseMatrix<double> matrix;
        Vector<double> rhs;
    };

    struct AssemblyScratchData {
        AssemblyScratchData(const FiniteElement <DIM> &fe,
                            const Quadrature <DIM> &quadrature);
        AssemblyScratchData(const AssemblyScratchData &scratch);
        FEValues<DIM> fe_values;
    };

    struct AssemblyCopyData {
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;

        std::vector<types::global_dof_index> local_dofs_indices;
    };

    void assemble_linear_system(LinearSystem &linear_system);

    void local_assemble_system(const typename DoFHandler<DIM>::active_cell_iterator &cell,
                               AssemblyScratchData &scratch_data,
                               AssemblyCopyData &copy_data) const;

    void copy_local_to_global(const AssemblyCopyData &copy_data,
                              LinearSystem &linear_system) const;
};

//end namespace ElasticitySolver
}


#endif //SOLVERS_ELASTICITYSOLVER_HPP

//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_LAPLACESOLVER_HPP
#define SOLVERS_LAPLACESOLVER_HPP

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <sstream>

#include "src/libsolvers/Postprocessor/EvaluationBase.hpp"

namespace LaplaceSolver
{
//begin namespace LaplaceSolver
using namespace dealii;

/**
 * @class Solver.
 * General solver class.
 */
template <int dim>
class GeneralSolver
{
public:
    GeneralSolver (Triangulation<dim> &triangulation,
                   const FiniteElement<dim> &fe,
                   const Quadrature<dim> &quadrature,
                   const Function<dim> &boundaryValues);
    virtual ~GeneralSolver ();

    virtual void solve_problem ();
    virtual void refine_grid () = 0;

    virtual size_t n_dofs () const;

    virtual void
    postprocess (const Postprocessor::EvaluationBase<dim> &postprocessor) const;

protected:
    const SmartPointer<Triangulation<dim>> triangulation;
    const SmartPointer<const FiniteElement<dim>> fe;
    const SmartPointer<const Quadrature<dim>> quadrature;
    DoFHandler<dim> dofHandler;

    Vector<double> solution;

    const SmartPointer<const Function<dim>> boundaryValues;

    virtual void assemble_rhs (Vector<double> &rhs) = 0;

private:
    struct LinearSystem
    {
        LinearSystem (const DoFHandler<dim> &dofHandler);

        void solve (Vector<double> &solution) const;

        ConstraintMatrix hangingNodeConstraints;
        SparsityPattern sparsityPattern;
        SparseMatrix<double> matrix;
        Vector<double> rhs;
    };

    struct AssemblyScratchData
    {
        AssemblyScratchData (const FiniteElement<dim> &fe,
                             const Quadrature<dim> &quadrature);
        AssemblyScratchData (const AssemblyScratchData &scratchData);

        FEValues<dim> feValues;
    };

    struct AssemblyCopyData
    {
        FullMatrix<double> cellMatrix;
        std::vector<types::global_dof_index> localDofIndices;
    };

    void assemble_linear_system (LinearSystem &linearSystem);

    void
    local_assemble_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                           AssemblyScratchData &scratchData,
                           AssemblyCopyData &copyData) const;

    void copy_local_to_global (const AssemblyCopyData &copyData,
                               LinearSystem &linearSystem) const;
};

/**
 * GeneralSolver constructor.
 * @param triangulation
 * @param fe
 * @param quadrature
 * @param boundaryValues
 */
template <int dim>
GeneralSolver<dim>::GeneralSolver (Triangulation<dim> &triangulation,
                                   const FiniteElement<dim> &fe,
                                   const Quadrature<dim> &quadrature,
                                   const Function<dim> &boundaryValues)
    :
    triangulation (&triangulation),
    fe (&fe),
    quadrature (&quadrature),
    dofHandler (triangulation),
    boundaryValues (&boundaryValues)
{

}
/**
 * GeneralSolver destructor.
 */
template <int dim>
GeneralSolver<dim>::~GeneralSolver ()
{
    dofHandler.clear ();
}

/**
 * Solve LaplaceProblem.
 */
template <int dim>
void GeneralSolver<dim>::solve_problem ()
{
    dofHandler.distribute_dofs (*fe);
    solution.reinit (dofHandler.n_dofs ());

    LinearSystem linearSystem (dofHandler);
    assemble_linear_system (linearSystem);
    linearSystem.solve (solution);
}

/**
 * Assemble linear system for Laplace problem.
 * @param linearSystem Linear system to assemble.
 * @details Uses Deal.II threads to run in parallel.
 */
template <int dim>
void
GeneralSolver<dim>::assemble_linear_system (GeneralSolver::LinearSystem &linearSystem)
{
    auto assembleRhsFunc = [this] (Vector<double> &rhs)
    {
        this->assemble_rhs (rhs);
    };

    Threads::Task<>
        rhsTask = Threads::new_task (assembleRhsFunc,
                                     linearSystem.rhs);

    auto assembleMatrixFunc = [this] (auto cell,
                                      auto scratchData,
                                      auto copyData)
    {
        this->local_assemble_matrix (cell,
                                     scratchData,
                                     copyData);
    };

    auto copyMatrixFunc = [this] (auto copyData,
                                  auto linearSystem)
    {
        this->copy_local_to_global (copyData,
                                    linearSystem);
    };


    WorkStream::run (dofHandler.begin_active (),
                     dofHandler.end (),
                     assembleMatrixFunc,
                     copyMatrixFunc,
                     AssemblyScratchData (*fe, *quadrature),
                     AssemblyCopyData ());

    linearSystem.hangingNodeConstraints.condense (linearSystem.matrix);

    std::map<types::global_dof_index, double> boundaryValueMap;
    VectorTools::interpolate_boundary_values (dofHandler,
                                              0,
                                              *boundaryValues,
                                              boundaryValueMap);
    rhsTask.join ();
    linearSystem.hangingNodeConstraints.condense (linearSystem.rhs);

    MatrixTools::apply_boundary_values (boundaryValueMap,
                                        linearSystem.matrix,
                                        solution,
                                        linearSystem.rhs);

}

/**
 * AssemblyScratchData constructor. Create data in which we store finite
 * element and quadrature data
 * @param fe
 * @param quadrature
 */
template <int dim>
GeneralSolver<dim>::AssemblyScratchData::AssemblyScratchData (const FiniteElement<dim> &fe,
                                                              const Quadrature<dim> &quadrature)
    :
    feValues (fe,
              quadrature,
              update_gradients | update_JxW_values)
{

}

/**
 * AssemblyScratchData copy constructor.
 * @param scratchData
 */
template <int dim>
GeneralSolver<dim>::AssemblyScratchData::AssemblyScratchData (const AssemblyScratchData &scratchData)
    :
    feValues (scratchData.feValues.get_fe (),
              scratchData.feValues.get_quadrature (),
              scratchData.feValues.get_update_flags ())
{

}

/**
 * Assemble local cell matrix.
 * @param cell Iterator to cell
 * @param scratchData Object which stores FEValue data
 * @param copyData Object in which we will write local matrix.
 */
template <int dim>
void
GeneralSolver<dim>::local_assemble_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           GeneralSolver::AssemblyScratchData &scratchData,
                                           GeneralSolver::AssemblyCopyData &copyData) const
{
    const size_t dofsPerCell = fe->dofs_per_cell;
    const size_t nQPoints = quadrature->size ();

    copyData.cellMatrix.reinit (dofsPerCell, dofsPerCell);
    copyData.localDofIndices.resize (dofsPerCell);
    scratchData.feValues.reinit (cell);

    for (size_t q = 0; q < nQPoints; ++q)
    {
        for (size_t i = 0; i < dofsPerCell; ++i)
        {
            for (size_t j = 0; j < dofsPerCell; ++j)
            {
                copyData.cellMatrix (i, j)
                    += scratchData.feValues.shape_grad (i, q) *
                       scratchData.feValues.shape_grad (j, q) *
                       scratchData.feValues.JxW (q);
            }
        }
    }

    cell->get_dof_indices (copyData.localDofIndices);
}

/**
 * Copy data from cell matrix to global matrix
 * @param copyData Data to copy
 * @param linearSystem System in which matrix we copy data.
 */
template <int dim>
void
GeneralSolver<dim>::copy_local_to_global (const GeneralSolver::AssemblyCopyData &copyData,
                                          GeneralSolver::LinearSystem &linearSystem) const
{
    const size_t nLocalDofIndices = copyData.localDofIndices.size ();
    for (size_t i = 0; i < nLocalDofIndices; ++i)
    {
        for (size_t j = 0; j < nLocalDofIndices; ++j)
        {
            linearSystem.matrix.add (copyData.localDofIndices[i],
                                     copyData.localDofIndices[j],
                                     copyData.cellMatrix (i, j));
        }
    }
}

/**
 * LinearSystem constructor.
 * Create linear system, matrix and right hand side, from degrees of freedom.
 * @param dofHandler
 */
template <int dim>
GeneralSolver<dim>::LinearSystem::LinearSystem (const DoFHandler<dim> &dofHandler)
{
    hangingNodeConstraints.clear ();

    auto mhnc = [] (const DoFHandler<dim> &d,
                    ConstraintMatrix &c)
    {
        DoFTools::make_hanging_node_constraints (d, c);
    };

    Threads::Task<> sideTask = Threads::new_task (mhnc,
                                                  dofHandler,
                                                  hangingNodeConstraints);

    DynamicSparsityPattern dsp (dofHandler.n_dofs (), dofHandler.n_dofs ());
    DoFTools::make_sparsity_pattern (dofHandler, dsp);

    sideTask.join ();

    hangingNodeConstraints.close ();
    hangingNodeConstraints.condense (dsp);
    sparsityPattern.copy_from (dsp);

    matrix.reinit (sparsityPattern);
    rhs.reinit (dofHandler.n_dofs ());
}

/**
 * Solve linear system.
 * @param solution Vector in which we write answer.
 */
template <int dim>
void GeneralSolver<dim>::LinearSystem::solve (Vector<double> &solution) const
{
    SolverControl solverControl (1000, 1e-12);
    SolverCG<> cg (solverControl);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize (matrix, 1.2);

    cg.solve (matrix, solution, rhs, preconditioner);

    hangingNodeConstraints.distribute (solution);
}

/**
 * Do postproccessing (e.g. output results, get value in some points)
 * @param postprocessor
 */
template <int dim>
void
GeneralSolver<dim>::postprocess (const Postprocessor::EvaluationBase<dim> &postprocessor) const
{
    postprocessor (dofHandler, solution);
}

/**
 * @return Number of degrees of freedom
 */
template <int dim>
size_t GeneralSolver<dim>::n_dofs () const
{
    return dofHandler.n_dofs ();
}

//end namespace LaplaceSolver
}


#endif //SOLVERS_LAPLACESOLVER_HPP

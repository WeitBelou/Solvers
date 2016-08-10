//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_PRIMALSOLVER_HPP
#define SOLVERS_PRIMALSOLVER_HPP
#include "GeneralSolver.hpp"
#include "src/libsolvers/Refinement/Refiner.hpp"

namespace LaplaceSolver
{
//begin namespace LaplaceSolver

/**
 * @class PrimalSolver.
 * Solve Laplace problem assuming that right hand side is
 * of problem given as functional object.
 */
template <int dim>
class PrimalSolver: public GeneralSolver<dim>
{
public:
    virtual void refine_grid () override;
public:
    PrimalSolver (Triangulation<dim> &triangulation,
                  const FiniteElement<dim> &fe,
                  const Quadrature<dim> &quadrature,
                  const Function<dim> &boundaryValues,
                  const Function<dim> &rhsFunction);

    virtual ~PrimalSolver ();

protected:
    const SmartPointer<const Function<dim>> rhsFunction;
    virtual void assemble_rhs (Vector<double> &rhs) override;
};

template <int dim>
PrimalSolver<dim>::PrimalSolver (Triangulation<dim> &triangulation,
                            const FiniteElement<dim> &fe,
                            const Quadrature<dim> &quadrature,
                            const Function<dim> &boundaryValues,
                            const Function<dim> &rhsFunction)
    :
    GeneralSolver<dim> (triangulation, fe, quadrature, boundaryValues),
    rhsFunction (&rhsFunction)
{

}

template <int dim>
PrimalSolver<dim>::~PrimalSolver ()
{

}

template <int dim>
void PrimalSolver<dim>::assemble_rhs (Vector<double> &rhs)
{
    FEValues<dim> feValues (*this->fe, *this->quadrature,
                            update_values | update_quadrature_points |
                            update_JxW_values);

    const size_t dofsPerCell = this->fe->dofs_per_cell;
    const size_t nQPoints = this->quadrature->size ();

    Vector<double> cellRhs (dofsPerCell);
    std::vector<double> rhsValues(nQPoints);
    std::vector<types::global_dof_index> localDofIndices (dofsPerCell);

    for (auto cell : this->dofHandler.active_cell_iterators ())
    {
        cellRhs = 0;
        feValues.reinit (cell);

        rhsFunction->value_list (feValues.get_quadrature_points (), rhsValues);

        for (size_t q = 0; q < nQPoints; ++q)
        {
            for (size_t i = 0; i < dofsPerCell; ++i)
            {
                cellRhs (i) += (feValues.shape_value (i, q)
                                * rhsValues[q]
                                * feValues.JxW (q));
            }
        }

        cell->get_dof_indices (localDofIndices);
        for (size_t i = 0; i < dofsPerCell; ++i)
        {
            rhs(localDofIndices[i]) += cellRhs(i);
        }
    }
}

template <int dim>
void PrimalSolver<dim>::refine_grid ()
{
    Refinement::GlobalRefiner<dim> refiner(*this->triangulation, 1);
    refiner.do_refinement ();
}

//end namespace LaplaceSolver
}
#endif //SOLVERS_PRIMALSOLVER_HPP

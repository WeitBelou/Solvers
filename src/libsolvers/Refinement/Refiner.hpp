//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_REFINER_HPP
#define SOLVERS_REFINER_HPP

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

namespace Refinement
{
//begin namespace Refinement
using namespace dealii;

/**
 * @class Refiner.
 * Abstract basic class for all refiners.
 */
template <int dim>
class Refiner
{
public:
    Refiner (Triangulation<dim> & triangulation);

    virtual ~Refiner() {}

    virtual void do_refinement () const = 0;

protected:
    const SmartPointer<Triangulation<dim>> triangulation;
};

template <int dim>
Refiner<dim>::Refiner (Triangulation<dim> & triangulation)
    :
    triangulation(&triangulation)
{

}

/**
 * @class GlobalRefiner
 * Refine grid "triangulation" "n" times.
 */
template <int dim>
class GlobalRefiner : public Refiner<dim>
{
public:
    GlobalRefiner (Triangulation<dim> & triangulation, int n);
    virtual void do_refinement () const override;
private:
    const int n;
};

template <int dim>
GlobalRefiner<dim>::GlobalRefiner (Triangulation<dim> &triangulation, int n)
    :
    Refiner<dim> (triangulation),
    n(n)
{

}

template <int dim>
void GlobalRefiner<dim>::do_refinement () const
{
    this->triangulation->refine_global (n);
}

//end namespace Refinement
}


#endif //SOLVERS_REFINER_HPP
